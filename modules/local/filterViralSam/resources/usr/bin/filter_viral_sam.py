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
    new_flag = old_flag
    
    # Flip read1/read2 bits
    if old_flag & 128:  # This is read2, create read1
        new_flag &= ~128    # Remove read2 flag
        new_flag |= 64      # Add read1 flag
    elif old_flag & 64:   # This is read1, create read2
        new_flag &= ~64     # Remove read1 flag
        new_flag |= 128     # Add read2 flag
    
    # Set this read as unmapped, keep mate mapped
    new_flag |= 4       # Set unmapped
    new_flag &= ~8      # Ensure mate mapped flag is clear
    
    # Update fields - keep same RNAME but unmapped characteristics
    fields[1] = str(new_flag)
    fields[2] = mapped_alignment['fields'][2]  # RNAME = same as mapped mate
    fields[3] = mapped_alignment['fields'][3]  # POS = same as mate (unmapped convention)
    fields[5] = '*'      # CIGAR = * (unmapped)
    fields[6] = '='      # RNEXT = same chromosome
    fields[7] = mapped_alignment['fields'][3]  # PNEXT = mate's position
    
    return '\t'.join(fields)

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
            # Skip header lines (should not be present based on our workflow)
            if line.startswith('@'):
                continue
                
            parsed = parse_sam_line(line)
            alignments_by_qname[parsed['qname']].append(parsed)
    
    # Process each read group
    filtered_alignments = []
    
    for qname in sorted(alignments_by_qname.keys()):
        alignments = alignments_by_qname[qname]
        
        # Skip if any read in this group is a contaminant
        if qname in contaminant_ids:
            continue
        
        # Separate primary and secondary alignments
        primary = [a for a in alignments if a['flag'] < 256]
        secondary = [a for a in alignments if a['flag'] >= 256]
        
        # Process primary alignments
        primary_kept = []
        for alignment in primary:
            query_length = len(alignment['seq'])
            normalized_score = calculate_normalized_score(alignment['alignment_score'], query_length)
            alignment['normalized_score'] = normalized_score
            primary_kept.append(alignment)
        
        # Apply pair-based score filtering for primary alignments
        if len(primary_kept) >= 2:
            # Check if either read in the pair exceeds threshold
            max_score = max(a['normalized_score'] for a in primary_kept)
            if max_score < score_threshold:
                primary_kept = []
        elif len(primary_kept) == 1:
            # Single read - check if it exceeds threshold
            if primary_kept[0]['normalized_score'] < score_threshold:
                primary_kept = []
        
        # Add missing mates for UP reads in primary
        for alignment in primary_kept[:]:  # Use slice copy to allow modification
            if alignment['pair_status'] == 'UP':
                # Check if mate exists
                mate_flag = 128 if (alignment['flag'] & 64) else 64
                has_mate = any(other['flag'] & mate_flag and other['flag'] < 256 
                              for other in primary_kept)
                
                if not has_mate:
                    unmapped_mate_line = create_unmapped_mate(alignment)
                    unmapped_mate = parse_sam_line(unmapped_mate_line)
                    unmapped_mate['normalized_score'] = 0  # Unmapped has no score
                    primary_kept.append(unmapped_mate)
        
        # Process secondary alignments similarly
        secondary_kept = []
        for alignment in secondary:
            query_length = len(alignment['seq'])
            normalized_score = calculate_normalized_score(alignment['alignment_score'], query_length)
            alignment['normalized_score'] = normalized_score
            secondary_kept.append(alignment)
        
        # Apply pair-based score filtering for secondary alignments by reference
        secondary_by_ref = defaultdict(list)
        for alignment in secondary_kept:
            secondary_by_ref[alignment['rname']].append(alignment)
        
        secondary_final = []
        for rname, ref_alignments in secondary_by_ref.items():
            # Apply score threshold check for this reference
            if len(ref_alignments) >= 2:
                max_score = max(a['normalized_score'] for a in ref_alignments)
                if max_score >= score_threshold:
                    secondary_final.extend(ref_alignments)
            elif len(ref_alignments) == 1:
                if ref_alignments[0]['normalized_score'] >= score_threshold:
                    secondary_final.extend(ref_alignments)
        
        # Add missing mates for UP reads in secondary
        secondary_by_ref_final = defaultdict(list)
        for alignment in secondary_final:
            secondary_by_ref_final[alignment['rname']].append(alignment)
        
        for rname, ref_alignments in secondary_by_ref_final.items():
            for alignment in ref_alignments[:]:  # Use slice copy
                if alignment['pair_status'] == 'UP':
                    # Check if mate exists for this reference
                    mate_flag = 128 if (alignment['flag'] & 64) else 64
                    has_mate = any(other['flag'] & mate_flag and other['rname'] == rname 
                                  for other in ref_alignments)
                    
                    if not has_mate:
                        unmapped_mate_line = create_unmapped_mate(alignment)
                        unmapped_mate = parse_sam_line(unmapped_mate_line)
                        unmapped_mate['normalized_score'] = 0
                        ref_alignments.append(unmapped_mate)
        
        # Combine and sort all kept alignments for this read
        all_kept = primary_kept + secondary_final
        
        # Sort: primary first (by flag), then secondary by reference name then flag
        primary_sorted = sorted([a for a in all_kept if a['flag'] < 256], 
                               key=lambda x: x['flag'])
        secondary_sorted = sorted([a for a in all_kept if a['flag'] >= 256], 
                                 key=lambda x: (x['rname'], x['flag']))
        
        filtered_alignments.extend(primary_sorted + secondary_sorted)
    
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