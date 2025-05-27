#!/usr/bin/env python3
import sys
from collections import defaultdict

def parse_sam_line(line):
    fields = line.strip().split('\t')
    
    pair_status = None
    for field in fields[11:]:
        if field.startswith('YT:Z:'):
            pair_status = field.split(':')[2]
            break
    
    return {
        'qname': fields[0],
        'flag': int(fields[1]),
        'rname': fields[2],
        'pos': int(fields[3]),
        'pnext': int(fields[7]),
        'pair_status': pair_status,
        'line': line,
        'fields': fields
    }

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
    fields[3] = mapped_alignment['fields'][3]      # POS = 0 (unmapped)
    fields[5] = '*'      # CIGAR = * (unmapped)
    fields[6] = '='      # RNEXT = same chromosome
    fields[7] = mapped_alignment['fields'][3]  # PNEXT = mate's position
    
    return '\t'.join(fields) + '\n'

def sort_sam_alignments(input_file, output_file):
    with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
        alignments_by_qname = defaultdict(list)
        
        for line in inf:
            if line.startswith('@'):
                outf.write(line)
            else:
                parsed = parse_sam_line(line)
                alignments_by_qname[parsed['qname']].append(parsed)
        
        for qname in sorted(alignments_by_qname.keys()):
            alignments = alignments_by_qname[qname]
            
            primary = [a for a in alignments if a['flag'] < 256]
            secondary = [a for a in alignments if a['flag'] >= 256]
            
            primary.sort(key=lambda x: x['flag'])
            
            # For secondary alignments, check if mate exists
            secondary_by_ref = defaultdict(list)
            for a in secondary:
                secondary_by_ref[a['rname']].append(a)
                
                # Check if this secondary alignment has a mate
                mate_flag = 128 if (a['flag'] & 64) else 64  # Expected mate flag
                has_mate = any(other['flag'] & mate_flag and 
                             other['rname'] == a['rname'] and
                             other['flag'] >= 256  # Also secondary
                             for other in alignments)
                
                if not has_mate and a['pair_status'] == 'UP':
                    unmapped_mate = create_unmapped_mate(a)
                    unmapped_parsed = parse_sam_line(unmapped_mate)
                    secondary_by_ref[a['rname']].append(unmapped_parsed)
            
            # Write output
            for alignment in primary:
                outf.write(alignment['line'])
            
            for rname in sorted(secondary_by_ref.keys()):
                ref_alignments = secondary_by_ref[rname]
                ref_alignments.sort(key=lambda x: x['flag'])
                for alignment in ref_alignments:
                    outf.write(alignment['line'])

if __name__ == '__main__':
    sort_sam_alignments(sys.argv[1], sys.argv[2])
