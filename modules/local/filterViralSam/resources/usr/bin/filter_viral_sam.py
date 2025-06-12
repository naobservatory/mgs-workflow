#!/usr/bin/env python3

import math
import logging
import gzip
import bz2
import argparse
import datetime
from datetime import datetime, timezone
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Iterator
from collections import defaultdict
from Bio import SeqIO


# Configure logging
class UTCFormatter(logging.Formatter):
    """
    Custom logging formatter that displays timestamps in UTC.
    
    Returns:
        Formatted log timestamps in UTC timezone
    """
    def formatTime(self, record, datefmt=None):
        """
        Format log timestamps in UTC timezone.
        
        Args:
            record: LogRecord object containing timestamp data
            datefmt: Optional date format string (unused)
            
        Returns:
            Formatted timestamp string in UTC
        """
        dt = datetime.fromtimestamp(record.created, timezone.utc)
        return dt.strftime("%Y-%m-%d %H:%M:%S UTC")


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = UTCFormatter("[%(asctime)s] %(message)s")
handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(handler)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for viral SAM filtering.
    
    Returns:
        Parsed arguments containing input_sam, filtered_fastq, output_sam, and score_threshold
    """
    desc = (
        "Given a sorted (by read ID) fastq file of filtered reads, a sorted (by read ID) SAM file with no header, "
        "and a minimum normalized alignment score threshold, filter the SAM file by the following steps:"
        "Keep only the filtered reads from the SAM file, then apply score thresholding on those reads, and finally"
        "add missing mates for unpaired reads."
    )
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("input_sam", help="Input SAM file to filter")
    parser.add_argument(
        "filtered_fastq", help="FASTQ file containing filtered reads to keep"
    )
    parser.add_argument("output_sam", help="Output filtered SAM file")
    parser.add_argument(
        "score_threshold",
        type=float,
        help="Minimum normalized alignment score threshold",
    )
    return parser.parse_args()


def open_by_suffix(filename: str, mode: str = "r"):
    """
    Parse the suffix of a filename to determine the right open method
    to use, then open the file. Can handle .gz, .bz2, and uncompressed files.
    
    Args:
        filename: Path to file to open
        mode: File open mode (default "r")
        
    Returns:
        File handle appropriate for the file compression type
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode + "t")
    elif filename.endswith(".bz2"):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


@dataclass
class SamAlignment:
    """
    Represents a single SAM alignment record with parsed fields and metadata.
    """

    qname: str  # query name/read id
    flag: int  # bitwise FLAG
    rname: str  # reference genome name
    pos: int  # alignment position
    mapq: int  # mapping quality
    cigar: str  # CIGAR string
    rnext: str  # reference genome name of the mate
    pnext: int  # alignment position of the mate
    tlen: int  # insert size
    seq: str  # sequence
    qual: str  # quality score for sequence
    pair_status: Optional[str]  # pair status
    alignment_score: Optional[int]  # alignment score
    normalized_score: float | None  # normalized alignment score
    line: str  # original SAM line
    fields: list[str]

    @classmethod
    def from_sam_line(cls, line: str) -> "SamAlignment":
        """
        Parse a SAM line and create a SamAlignment object.
        
        Args:
            line: Raw SAM format line
            
        Returns:
            SamAlignment object with parsed fields and metadata
        """
        fields = line.strip().split("\t")
        pair_status = None
        alignment_score = None
        # Extract YT (pair status) and AS (alignment score) tags
        for field in fields[11:]:
            if field.startswith("YT:Z:"):
                pair_status = field.split(":")[2]
            elif field.startswith("AS:i:"):
                alignment_score = int(field.split(":")[2])
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
            normalized_score=None,  # Will be calculated later
            line=line,
            fields=fields,
        )

    def calculate_normalized_score(self) -> None:
        """
        Calculate and set the normalized alignment score.
        
        Normalizes the alignment score by dividing by the log of sequence length.
        Sets normalized_score to 0.0 if alignment_score is None or sequence is empty.
        
        Returns:
            None (modifies self.normalized_score in place)
        """
        if self.alignment_score is None or len(self.seq) == 0:
            self.normalized_score = 0.0
        else:
            self.normalized_score = self.alignment_score / math.log(len(self.seq))

    def create_unmapped_mate(self) -> "SamAlignment":
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
        - AND with ~8 (NOT 0b1000) clears the mate unmapped flag

        Returns:
            SamAlignment object representing the synthetic unmapped mate

        """
        fields = self.fields[:]
        old_flag = self.flag
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
        fields[5] = "*"  # CIGAR = * (unmapped)
        fields[6] = "="  # RNEXT = same chromosome as mate
        fields[7] = self.fields[3]  # PNEXT = mate's position
        mate_line = "\t".join(fields)
        mate = SamAlignment.from_sam_line(mate_line)
        return mate


def group_alignments_by_mates(
    alignments: list[SamAlignment],
) -> Dict[int, list[SamAlignment]]:
    """
    Group alignments by mate pairs.
    
    For UP (unpaired) reads, groups alignments by position/mate information.
    For CP (concordant pair) reads, groups by position and reference name since both pairs must be present.
    
    Args:
        alignments: List of SamAlignment objects to group
        
    Returns:
        Dictionary mapping group IDs to lists of alignments in each group
    """
    if not alignments:
        return {}
    
    # Check pair status - all alignments in a group should have the same pair_status
    pair_status = alignments[0].pair_status
    
    by_group: Dict[int, list[SamAlignment]] = defaultdict(list)
    group_idx = 0
    
    if pair_status == "UP":
        # Group by unique key: (rname, rnext, min(pos, pnext), max(pos, pnext))
        pair_key_to_group = {}
        
        for alignment in alignments:
            # Skip unmapped or self-referential
            if alignment.pnext == 0 or alignment.pos == alignment.pnext:
                by_group[group_idx].append(alignment)
                group_idx += 1
                continue
            
            # Skip cross-reference genome mappings
            if alignment.rnext != "=" and alignment.rname != alignment.rnext:
                by_group[group_idx].append(alignment)
                group_idx += 1
                continue
            
            # Create key for same-reference mappings
            pos_min = min(alignment.pos, alignment.pnext)
            pos_max = max(alignment.pos, alignment.pnext)
            
            pair_key = (
                alignment.qname,
                alignment.rname,
                pos_min,
                pos_max
            )
            
            if pair_key not in pair_key_to_group:
                pair_key_to_group[pair_key] = group_idx
                group_idx += 1
            
            group_id = pair_key_to_group[pair_key]
            by_group[group_id].append(alignment)
    else:
        # For non-UP reads, map tuple keys to int keys
        pair_key_to_group: Dict[tuple, int] = {}
        
        for alignment in alignments:
            pos_min = min(alignment.pos, alignment.pnext)
            pos_max = max(alignment.pos, alignment.pnext)
            abs_tlen = abs(alignment.tlen)
            
            pair_key = (
                alignment.qname,
                alignment.rname,
                alignment.rnext,
                pos_min,
                pos_max,
                abs_tlen
            )
            
            # Get or create group_idx for this pair_key
            if pair_key not in pair_key_to_group:
                pair_key_to_group[pair_key] = group_idx
                group_idx += 1
            
            group_id = pair_key_to_group[pair_key]
            by_group[group_id].append(alignment)

    return dict(by_group)


def apply_score_filter(
    alignments: list[SamAlignment], threshold: float, secondary: bool = False
) -> Dict[int, list[SamAlignment]]:
    """
    Filter alignment groups based on normalized score threshold.
    
    For secondary alignments, groups alignments by mates and keeps groups where
    at least one alignment exceeds the threshold. For primary alignments,
    keeps all alignments if any exceeds the threshold.
    
    Args:
        alignments: List of SamAlignment objects to filter
        threshold: Minimum normalized score threshold
        secondary: Whether these are secondary alignments (affects grouping logic)
        
    Returns:
        Dictionary of group IDs to alignment lists that pass the threshold
    """
    if secondary:
        groups = group_alignments_by_mates(alignments)
        kept_groups = {}

        for group_id, group_alignments in groups.items():
            max_score = max(
                a.normalized_score
                for a in group_alignments
                if a.normalized_score is not None
            )
            if max_score >= threshold:
                kept_groups[group_id] = group_alignments
            else:
                logger.info("Filtering group %s", group_id)

        return kept_groups
    else:
        max_score = max(
            a.normalized_score for a in alignments if a.normalized_score is not None
        )
        if max_score >= threshold:
            return {0: alignments}  # Single group for primary
        else:
            return {}


def add_missing_mates_by_groups(
    groups: Dict[int, list[SamAlignment]],
) -> Dict[int, list[SamAlignment]]:
    """
    Add missing mate reads based on provided alignment groups.
    
    For each group, checks if both read1 and read2 are present. If either is missing
    and the existing alignment has pair_status "UP" (unpaired), creates a synthetic
    unmapped mate to complete the pair.
    
    Args:
        groups: Dictionary mapping group IDs to lists of alignments
        
    Returns:
        Dictionary mapping group IDs to lists of alignments including synthetic mates
    """
    result = {}

    for group_id, group_alignments in groups.items():
        # Check if we have both read1 and read2 in this group
        has_read1 = any(a.flag & 64 for a in group_alignments)
        has_read2 = any(a.flag & 128 for a in group_alignments)

        # Start with existing alignments
        group_result = group_alignments[:]

        # Only create synthetic mates if we're missing either read1 or read2
        if not (has_read1 and has_read2):
            for alignment in group_alignments:
                if alignment.pair_status == "UP":
                    unmapped_mate = alignment.create_unmapped_mate()
                    group_result.append(unmapped_mate)

        # Sort this group by flag
        result[group_id] = sorted(group_result, key=lambda x: x.flag)

    return result


def process_alignment_group(
    alignments: list[SamAlignment], score_threshold: float
) -> list[SamAlignment]:
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
    assert len(primary) <= 2
    secondary = [a for a in alignments if a.flag >= 256]

    # Calculate normalized scores for all alignments
    for alignment in primary + secondary:
        alignment.calculate_normalized_score()

    # Get groups that pass score threshold
    primary_groups = apply_score_filter(primary, score_threshold)

    if primary_groups == {}:
        logger.info("Filtering primary group")
        return []

    secondary_groups = apply_score_filter(secondary, score_threshold, secondary=True)

    # Add missing mates using the groups (returns sorted groups)
    primary_with_mates = add_missing_mates_by_groups(primary_groups)
    secondary_with_mates = add_missing_mates_by_groups(secondary_groups)

    # Combine all groups into final result
    result = []
    for group_alignments in primary_with_mates.values():
        result.extend(group_alignments)
    for group_alignments in secondary_with_mates.values():
        result.extend(group_alignments)
    
    return result


def stream_filtered_fastq(fastq_file: str) -> Iterator[str]:
    """
    Stream filtered IDs one at a time from sorted FASTQ file.
    Assumes FASTQ is sorted by read ID and handles paired reads correctly.

    For interleaved reads (read1/read2), yields the base read ID only once to avoid
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
                read_id = record.id.split()[0].split("/")[0]

                # Only yield if this is a new read ID (handles paired reads)
                if read_id != previous_read_id:
                    yield read_id
                    previous_read_id = read_id
    except Exception as e:
        logger.error(f"Error reading FASTQ file {fastq_file}: {e}")
        raise


def stream_sam_by_qname(sam_file: str) -> Iterator[Tuple[str, list[SamAlignment]]]:
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
            if isinstance(line, bytes):
                line = line.decode()
            if line.startswith("@"):
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


def filter_viral_sam(
    input_sam: str, filtered_fastq: str, output_sam: str, score_threshold: float
) -> None:
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
    filtered_read_ids_iter = stream_filtered_fastq(filtered_fastq)
    sam_iter = stream_sam_by_qname(input_sam)
    # Two-pointer merge algorithm
    curr_read_id = next(filtered_read_ids_iter, None)
    alignments_processed = 0
    alignments_kept = 0
    groups_processed = 0
    non_filtered_skipped = 0
    filt_bool = False
    logger.info(f"Processing SAM file with score threshold {score_threshold}")
    with open_by_suffix(output_sam, "w") as outf:
        last_curr_align_read_id = None
        for curr_align_read_id, alignments in sam_iter:
            # Check if SAM file is sorted
            if (
                last_curr_align_read_id is not None
                and curr_align_read_id < last_curr_align_read_id
            ):
                logger.error(
                    f"SAM file not sorted! {curr_align_read_id} appeared after {last_curr_align_read_id}"
                )
                raise ValueError(
                    f"SAM file not sorted by curr_align_read_id at {curr_align_read_id}, after {last_curr_align_read_id}"
                )
            last_curr_align_read_id = curr_align_read_id
            alignments_processed += len(alignments)
            groups_processed += 1
            # Advance filtered iterator until >= curr_align_read_id
            while curr_read_id is not None and curr_read_id < curr_align_read_id:
                if filt_bool:
                  logger.debug(f"Advancing past filtered ID: {curr_read_id}")
                  curr_read_id = next(filtered_read_ids_iter, None)
                  filt_bool = False
                else:
                  raise ValueError(
                        f"SAM file not sorted by curr_align_read_id at {curr_align_read_id}, after {last_curr_align_read_id}"
                    )

            # Skip if this curr_align_read_id is not in filtered reads
            if curr_read_id is not None and curr_read_id > curr_align_read_id:
                  logger.debug(
                      f"Skipping alignment that is not found in fastq file: {curr_align_read_id}"
                  )
                  non_filtered_skipped += 1
                  filt_bool = False
                  continue

            filt_bool = True
            # Process and write immediately
            logger.info(f"Starting {curr_align_read_id}")
            kept_alignments = process_alignment_group(alignments, score_threshold)
            alignments_kept += len(kept_alignments)
            # Write alignments sorted by flag
            for alignment in kept_alignments:
                outf.write(alignment.line)
                if not alignment.line.endswith("\n"):
                    outf.write("\n")
    logger.info(
        f"Processed {groups_processed} read groups, {alignments_processed} alignments"
    )
    logger.info(
        f"Skipped {non_filtered_skipped} non-filtered groups, kept {alignments_kept} alignments"
    )


def main() -> None:
    """
    Main entry point for viral SAM filtering script.
    
    Parses command-line arguments and executes the filtering pipeline.
    
    Returns:
        None
    """
    logger.info("Initializing script.")
    logger.info("Parsing arguments.")
    args = parse_args()
    filter_viral_sam(
        args.input_sam, args.filtered_fastq, args.output_sam, args.score_threshold
    )


if __name__ == "__main__":
    main()
