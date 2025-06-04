#!/usr/bin/env python3

import math
import logging
import gzip
import argparse
import datetime
from datetime import datetime, timezone
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Iterator
from collections import defaultdict
from Bio import SeqIO

# Configure logging
class UTCFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime('%Y-%m-%d %H:%M:%S UTC')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter('[%(asctime)s] %(message)s')
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)


def parse_args() -> argparse.Namespace:

    desc = "Given a sorted (by read ID) fastq file of filtered reads, a sorted (by read ID) SAM file with no header, " \
           "and a minimum normalized alignment score threshold, filter the SAM file by the following steps:" \
           "Keep only the filtered reads from the SAM file, then apply score thresholding on those reads, and finally" \
           "add missing mates for unpaired reads."

    parser = argparse.ArgumentParser(
        description=desc
    )
    parser.add_argument(
        'input_sam',
        help='Input SAM file to filter'
    )
    parser.add_argument(
        'filtered_fastq',
        help='FASTQ file containing filtered reads to keep'
    )
    parser.add_argument(
        'output_sam',
        help='Output filtered SAM file'
    )
    parser.add_argument(
        'score_threshold',
        type=float,
        help='Minimum normalized alignment score threshold'
    )
    
    return parser.parse_args()

def open_by_suffix(filename, mode="r", debug=False):
    """Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


@dataclass
class SamAlignment:
    """
    Represents a single SAM alignment record with parsed fields and metadata.
    """
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    pair_status: Optional[str]
    alignment_score: Optional[int]
    normalized_score: float
    line: str
    fields: List[str]

    @classmethod
    def from_sam_line(cls, line: str) -> 'SamAlignment':
        """Parse a SAM line and create a SamAlignment object."""
        fields = line.strip().split('\t')
        
        pair_status = None
        alignment_score = None
        
        # Extract YT (pair status) and AS (alignment score) tags
        for field in fields[11:]:
            if field.startswith('YT:Z:'):
                pair_status = field.split(':')[2]
            elif field.startswith('AS:i:'):
                alignment_score = int(field.split(':')[2])
        
        return cls(
            qname=fields[0],
            flag=int(fields[1]),
            rname=fields[2],
            pos=int(fields[3]),
            mapq=int(fields[4]),
            cigar=fields[5],
            rnext=fields[6],
            pnext=int(fields[7]),
            tlen=int(fields[8]),
            seq=fields[9],
            qual=fields[10],
            pair_status=pair_status,
            alignment_score=alignment_score,
            normalized_score=0.0,  # Will be calculated later
            line=line,
            fields=fields
        )

    def calculate_normalized_score(self) -> None:
        """Calculate and set the normalized alignment score."""
        if self.alignment_score is None or len(self.seq) == 0:
            self.normalized_score = 0.0
        else:
            self.normalized_score = self.alignment_score / math.log(len(self.seq))

    def create_unmapped_mate(self) -> 'SamAlignment':
        """
        Create a synthetic unmapped mate for a read in a SAM file.

        When Bowtie2 is run with secondary alignments enabled, its SAM output can cause issues for our downstream script, `process_bowtie_viral_sam.py`. This script processes the SAM file line by line, expecting that consecutive entries with the same `read_id` contain consistent metadata in key columns (such as alignment status and FLAG).

        The problem arises specifically when secondary alignments are present. In some cases, Bowtie2 will output an unmapped primary alignment for one mate, but then also output secondary alignments, all sharing the same `read_id`. This causes `process_bowtie_viral_sam.py` to encounter consecutive entries with the same `read_id`, but conflicting metadata—e.g., the first entry is marked as unmapped, however the second entry has the same `read_id`. This inconsistency leads the script to throw an error.

        To address this, we generate a synthetic unmapped mate entry whenever its absence would create metadata conflicts for entries with the same `read_id`. This synthesized entry uses information from the mapped mate, but key fields—most importantly the FLAG—are updated to accurately indicate an unmapped status. This ensures that every set of consecutive reads with the same `read_id` has consistent metadata, allowing the downstream script to operate without errors.

        Those interested in learning more about SAM format can refer to this, https://samtools.github.io/hts-specs/SAMv1.pdf.
        
        This is used for UP (unpaired) reads that need synthetic unmapped mates.
        The SAM flag manipulation:
        - XOR with 192 (0b11000000) flips read1/read2 bits (64|128)
        - OR with 4 (0b100) sets the unmapped flag
        - AND with ~8 (NOT 0b1000) clears the mate mapped flag
        """
        fields = self.fields[:]
        old_flag = self.flag
        
        # SAM flag manipulation to create synthetic unmapped mate records
        # Note: The FLAG column in the SAM file uses bits to indicate information about the alignment.
        # SAM flag bits 0x40 (64) and 0x80 (128) indicate read position within the template:
        # 0x40: first segment in template (typically read1 in paired-end)  
        # 0x80: last segment in template (typically read2 in paired-end)
        # To create a synthetic unmapped mate, we flip these bits so the mate appears
        # as the complementary read in the pair (read1 ↔ read2).
        # We set the unmapped flag (0x4) and clear the mate unmapped flag (0x8)
        # to indicate this read is unmapped while its mate has a valid alignment.
        new_flag = (old_flag ^ 192) | 4  # XOR flips bits 6&7 (64|128), OR adds unmapped
        new_flag &= ~8  # Ensure mate unmapped flag is clear
        
        # Update fields for unmapped mate
        fields[1] = str(new_flag)
        fields[5] = '*'      # CIGAR = * (unmapped)
        fields[6] = '='      # RNEXT = same chromosome as mate
        fields[7] = self.fields[3]  # PNEXT = mate's position
        
        mate_line = '\t'.join(fields)
        mate = SamAlignment.from_sam_line(mate_line)
        mate.normalized_score = 0.0
        
        return mate

def stream_filtered_ids(fastq_file: str) -> Iterator[str]:
    """
    Stream filtered IDs one at a time from sorted FASTQ file.
    Assumes FASTQ is sorted by read ID and handles paired reads correctly.
    
    For paired reads (read1/read2), yields the base read ID only once to avoid
    duplicates when comparing with SAM query names.
    
    Args:
        fastq_file: Path to FASTQ file (can be gzipped)
        
    Yields:
        Read IDs (without @ prefix and /1, /2 suffixes) in sorted order
    """
    previous_read_id = None
    
    try:
        with open_by_suffix(fastq_file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Remove /1, /2 suffixes and space-separated parts if present
                read_id = record.id.split()[0].split('/')[0]
                
                # Only yield if this is a new read ID (handles paired reads)
                if read_id != previous_read_id:
                    yield read_id
                    previous_read_id = read_id
    except Exception as e:
        logger.error(f"Error reading FASTQ file {fastq_file}: {e}")
        raise

def apply_score_filter(alignments: List[SamAlignment], threshold: float, group_by_ref: bool = False) -> List[SamAlignment]:
    """
    Apply score threshold filtering to alignments.
    
    Args:
        alignments: List of alignments to filter
        threshold: Minimum normalized score threshold
        group_by_ref: If True, group alignments by reference sequence (rname) and apply 
                     threshold per group. If False, apply threshold to the entire group.
                     
    Reference grouping is used for secondary alignments where each reference 
    represents a different viral genome that the read mapped to.
    """
    if group_by_ref:
        # Group by reference sequence and apply threshold per group
        by_ref: Dict[str, List[SamAlignment]] = defaultdict(list)
        for alignment in alignments:
            by_ref[alignment.rname].append(alignment)
        
        kept = []
        for ref_alignments in by_ref.values():
            if len(ref_alignments) == 2:
                max_score = max(a.normalized_score for a in ref_alignments)
                if max_score >= threshold:
                    kept.extend(ref_alignments)
            elif len(ref_alignments) == 1:
                if ref_alignments[0].normalized_score >= threshold:
                    kept.extend(ref_alignments)
        return kept
    else:
        # Apply threshold to group as whole
        if len(alignments) == 2:
            max_score = max(a.normalized_score for a in alignments)
            return alignments if max_score >= threshold else []
        elif len(alignments) == 1:
            return alignments if alignments[0].normalized_score >= threshold else []
        return []

def add_missing_mates(alignments: List[SamAlignment], group_by_ref: bool = False) -> None:
    """
    Add missing mates for UP (unpaired) reads.
    
    Args:
        alignments: List of alignments to process (modified in place)
        group_by_ref: If True, only consider mates within the same reference group
    """
    if group_by_ref:
        by_ref: Dict[str, List[SamAlignment]] = defaultdict(list)
        for alignment in alignments:
            by_ref[alignment.rname].append(alignment)
        
        for ref_alignments in by_ref.values():
            _add_mates_to_group(ref_alignments, 
                              lambda other, mate_flag, rname: other.flag & mate_flag and other.rname == rname)
        
        # Propagate changes back to main alignments list
        alignments.clear()
        for ref_alignments in by_ref.values():
            alignments.extend(ref_alignments)
    else:
        _add_mates_to_group(alignments, 
                           lambda other, mate_flag, rname: other.flag & mate_flag and other.flag < 256)

def _add_mates_to_group(alignments: List[SamAlignment], mate_check_fn) -> None:
    """Helper to add mates to a group of alignments."""
    for alignment in alignments[:]:  # Use slice copy to avoid modifying list during iteration
        if alignment.pair_status == 'UP':
            # Determine which mate flag to look for (read1 vs read2)
            mate_flag = 128 if (alignment.flag & 64) else 64
            has_mate = any(mate_check_fn(other, mate_flag, alignment.rname) 
                          for other in alignments)
            
            if not has_mate:
                unmapped_mate = alignment.create_unmapped_mate()
                alignments.append(unmapped_mate)

def process_alignment_group(alignments: List[SamAlignment], score_threshold: float) -> List[SamAlignment]:
    """
    Process all alignments for a single read name.
    
    Args:
        alignments: All alignments for one read name
        score_threshold: Minimum normalized score threshold
        
    Returns:
        Filtered and processed alignments
    """
    # Separate primary and secondary alignments (flag < 256 vs >= 256)
    primary = [a for a in alignments if a.flag < 256]
    secondary = [a for a in alignments if a.flag >= 256]
    
    # Calculate normalized scores for all alignments
    for alignment in primary + secondary:
        alignment.calculate_normalized_score()
    
    # Apply score filtering
    primary_kept = apply_score_filter(primary, score_threshold)
    secondary_kept = apply_score_filter(secondary, score_threshold, group_by_ref=True)
    
    # Add missing mates
    add_missing_mates(primary_kept)
    add_missing_mates(secondary_kept, group_by_ref=True)
    
    # Sort and combine results
    primary_sorted = sorted(primary_kept, key=lambda x: x.flag)
    secondary_sorted = sorted(secondary_kept, key=lambda x: (x.rname, x.flag))
    
    return primary_sorted + secondary_sorted

def stream_sam_by_qname(sam_file: str) -> Iterator[Tuple[str, List[SamAlignment]]]:
    """
    Stream SAM file and yield groups of alignments by query name.
    
    This memory-efficient approach assumes the SAM file is sorted by query name,
    which is typically the case for Bowtie2 output.
    
    Args:
        sam_file: Path to SAM file
        
    Yields:
        Tuples of (qname, list_of_alignments)
    """
    current_qname = None
    current_alignments = []
    
    with open_by_suffix(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header lines
                
            alignment = SamAlignment.from_sam_line(line)
            
            if current_qname is None:
                current_qname = alignment.qname
                current_alignments = [alignment]
            elif alignment.qname == current_qname:
                current_alignments.append(alignment)
            else:
                # New query name, yield previous group
                yield current_qname, current_alignments
                current_qname = alignment.qname
                current_alignments = [alignment]
        
        # Yield final group
        if current_qname is not None:
            yield current_qname, current_alignments

def filter_viral_sam_memory_efficient(input_sam: str, filtered_fastq: str, output_sam: str, score_threshold: float) -> None:
    """
    Filter viral SAM file using a memory-efficient streaming approach.
    
    Processing steps:
    1. Stream filtered read IDs and SAM alignments simultaneously
    2. Use two-pointer merge to keep only filtered reads
    3. Apply score threshold filtering (keep pair if either read exceeds threshold)
    4. Add missing mates for UP reads
    5. Write output immediately without accumulation
    
    Assumes both filtered FASTQ and SAM files are sorted by query ID.
    
    Args:
        input_sam: Input SAM file path (sorted by query name)
        filtered_fastq: FASTQ file containing filtered reads to keep (sorted by read ID)
        output_sam: Output filtered SAM file path
        score_threshold: Minimum normalized alignment score threshold
    """
    logger.info("Starting memory-efficient viral SAM filtering")
    
    # Create iterators for streaming processing
    filtered_iter = stream_filtered_ids(filtered_fastq)
    sam_iter = stream_sam_by_qname(input_sam)
    
    # Two-pointer merge algorithm
    current_filtered = next(filtered_iter, None)
    alignments_processed = 0
    alignments_kept = 0
    groups_processed = 0
    non_filtered_skipped = 0
    
    logger.info(f"Processing SAM file with score threshold {score_threshold}")
    
    with open_by_suffix(output_sam, 'w') as outf:
        for qname, alignments in sam_iter:
            alignments_processed += len(alignments)
            groups_processed += 1
            
            # Advance filtered iterator until >= qname
            while current_filtered is not None and current_filtered < qname:
                logger.debug(f"Advancing past filtered ID: {current_filtered}")
                current_filtered = next(filtered_iter, None)
            
            # Skip if this qname is not in filtered reads
            if current_filtered != qname:
                logger.debug(f"Skipping non-filtered read: {qname}")
                non_filtered_skipped += 1
                continue
                
            # Process and write immediately
            kept_alignments = process_alignment_group(alignments, score_threshold)
            alignments_kept += len(kept_alignments)
            
            if len(kept_alignments) != len(alignments):
                logger.debug(f"Score filtering: {qname} kept {len(kept_alignments)}/{len(alignments)} alignments")
            
            # Write alignments sorted by flag
            for alignment in sorted(kept_alignments, key=lambda x: (x.rname, x.flag)):
                outf.write(alignment.line)
                if not alignment.line.endswith('\n'):
                    outf.write('\n')
    
    logger.info(f"Processed {groups_processed} read groups, {alignments_processed} alignments")
    logger.info(f"Skipped {non_filtered_skipped} non-filtered groups, kept {alignments_kept} alignments")

def main() -> None:
    logger.info("Initializing script.")
    logger.info("Parsing arguments.")
    args = parse_args()
    filter_viral_sam_memory_efficient(
        args.input_sam,
        args.filtered_fastq,
        args.output_sam,
        args.score_threshold
    )

if __name__ == "__main__":
    main()
