#!/usr/bin/env python3

import sys
import math
from collections import defaultdict

def parse_sam_line(line):
    """Parse a SAM line and extract relevant fields"""
    fields = line.strip().split('\t')
    
    pair_status = None
    alignment_score = None
    
    # Extract YT (pair status) and AS (alignment score) tags
    for field in fields[11:]:
        if field.startswith('YT:Z:'):
            pair_status = field.split(':')[2]
        elif field.startswith('AS:i:'):
            alignment_score = int(field.split(':')[2])
    
    return {
        'qname': fields[0],
        'flag': int(fields[1]),
        'rname': fields[2],
        'pos': int(fields[3]),
        'mapq': int(fields[4]),
        'cigar': fields[5],
        'rnext': fields[6],
        'pnext': int(fields[7]),
        'tlen': int(fields[8]),
        'seq': fields[9],
        'qual': fields[10],
        'pair_status': pair_status,
        'alignment_score': alignment_score,
        'line': line,
        'fields': fields
    }

def calculate_normalized_score(alignment_score, query_length):
    """Calculate normalized alignment score"""
    if alignment_score is None or query_length == 0:
        return 0
    return alignment_score / math.log(query_length)

def create_unmapped_mate(mapped_alignment):
    """Create unmapped mate for a secondary alignment missing its pair"""
    fields = mapped_alignment['fields'][:]
    old_flag = mapped_alignment['flag']
    
    # Flip read1/read2 bits and set unmapped
    new_flag = (old_flag ^ 192) | 4  # XOR flips bits 6&7 (64|128), OR adds unmapped
    new_flag &= ~8  # Ensure mate mapped flag is clear
    
    # Update fields for unmapped mate
    fields[1] = str(new_flag)
    fields[5] = '*'      # CIGAR = * (unmapped)
    fields[6] = '='      # RNEXT = same chromosome
    fields[7] = mapped_alignment['fields'][3]  # PNEXT = mate's position
    
    return '\t'.join(fields)

def apply_score_filter(alignments, threshold, group_by_ref=False):
    """Apply score threshold filtering to alignments"""
    if group_by_ref:
        # Group by reference and apply threshold per group
        by_ref = defaultdict(list)
        for alignment in alignments:
            by_ref[alignment['rname']].append(alignment)
        
        kept = []
        for ref_alignments in by_ref.values():
            if len(ref_alignments) == 2:
                max_score = max(a['normalized_score'] for a in ref_alignments)
                if max_score >= threshold:
                    kept.extend(ref_alignments)
            elif len(ref_alignments) == 1:
                if ref_alignments[0]['normalized_score'] >= threshold:
                    kept.extend(ref_alignments)
        return kept
    else:
        # Apply threshold to group as whole
        if len(alignments) == 2:
            max_score = max(a['normalized_score'] for a in alignments)
            return alignments if max_score >= threshold else []
        elif len(alignments) == 1:
            return alignments if alignments[0]['normalized_score'] >= threshold else []
        return []

def add_missing_mates(alignments, group_by_ref=False):
    """Add missing mates for UP reads"""
    if group_by_ref:
        by_ref = defaultdict(list)
        for alignment in alignments:
            by_ref[alignment['rname']].append(alignment)
        
        for ref_alignments in by_ref.values():
            _add_mates_to_group(ref_alignments, lambda other, mate_flag, rname: 
                              other['flag'] & mate_flag and other['rname'] == rname)
    else:
        _add_mates_to_group(alignments, lambda other, mate_flag, rname: 
                           other['flag'] & mate_flag and other['flag'] < 256)

def _add_mates_to_group(alignments, mate_check_fn):
    """Helper to add mates to a group of alignments"""
    for alignment in alignments[:]:  # Use slice copy
        if alignment['pair_status'] == 'UP':
            mate_flag = 128 if (alignment['flag'] & 64) else 64
            has_mate = any(mate_check_fn(other, mate_flag, alignment['rname']) 
                          for other in alignments)
            
            if not has_mate:
                unmapped_mate_line = create_unmapped_mate(alignment)
                unmapped_mate = parse_sam_line(unmapped_mate_line)
                unmapped_mate['normalized_score'] = 0
                alignments.append(unmapped_mate)

def process_alignment_group(alignments, score_threshold):
    """Process all alignments for a single read name"""
    # Separate primary and secondary alignments
    primary = [a for a in alignments if a['flag'] < 256]
    secondary = [a for a in alignments if a['flag'] >= 256]
    
    # Calculate normalized scores for all alignments
    for alignment in primary + secondary:
        query_length = len(alignment['seq'])
        alignment['normalized_score'] = calculate_normalized_score(
            alignment['alignment_score'], query_length)
    
    # Apply score filtering
    primary_kept = apply_score_filter(primary, score_threshold)
    secondary_kept = apply_score_filter(secondary, score_threshold, group_by_ref=True)
    
    # Add missing mates
    add_missing_mates(primary_kept)
    add_missing_mates(secondary_kept, group_by_ref=True)
    
    # Sort and combine results
    primary_sorted = sorted(primary_kept, key=lambda x: x['flag'])
    secondary_sorted = sorted(secondary_kept, key=lambda x: (x['rname'], x['flag']))
    
    return primary_sorted + secondary_sorted

def filter_viral_sam(input_sam, contaminant_ids_file, output_sam, score_threshold):
    """
    Filter viral SAM file by:
    1. Removing contaminant reads
    2. Applying score threshold (keep pair if either read exceeds threshold)
    3. Adding missing mates for UP reads
    4. Sorting output
    """
    
    # Load contaminant read IDs
    contaminant_ids = set()
    with open(contaminant_ids_file, 'r') as f:
        for line in f:
            contaminant_ids.add(line.strip())
    
    # Group alignments by read name
    alignments_by_qname = defaultdict(list)
    
    with open(input_sam, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
                
            parsed = parse_sam_line(line)
            alignments_by_qname[parsed['qname']].append(parsed)
    
    # Process each read group
    filtered_alignments = []
    
    for qname in sorted(alignments_by_qname.keys()):
        if qname in contaminant_ids:
            continue
        
        alignments = alignments_by_qname[qname]
        kept_alignments = process_alignment_group(alignments, score_threshold)
        filtered_alignments.extend(kept_alignments)
    
    # Write output
    with open(output_sam, 'w') as f:
        for alignment in filtered_alignments:
            f.write(alignment['line'])
            if not alignment['line'].endswith('\n'):
                f.write('\n')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: filter_viral_sam.py <input_sam> <contaminant_ids> <output_sam> <score_threshold>", file=sys.stderr)
        sys.exit(1)
    
    input_sam = sys.argv[1]
    contaminant_ids_file = sys.argv[2]
    output_sam = sys.argv[3]
    score_threshold = float(sys.argv[4])
    
    filter_viral_sam(input_sam, contaminant_ids_file, output_sam, score_threshold)
