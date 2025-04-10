/*
Take in reads in either single-end or interleaved format, then handle
them based on endedness:

- Single-end reads are just passed out again as-is
- Interleaved reads are merged with BBMERGE, then unmerged reads are
    concatenated (with an intervening "N"), to produce a single output
    sequence per input read pair.

The convoluted way of handling the single-end vs interleaved cases
is necessitated by the fact that `single_end` is a channel, not a bare boolean.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ } from "../../../modules/local/joinFastq"
include { SUMMARIZE_BBMERGE } from "../../../modules/local/summarizeBBMerge"

/***********
| WORKFLOW |
***********/

workflow MERGE_JOIN_READS {
    take:
        reads_ch // Single-end or interleaved FASTQ sequences
        single_end // Boolean channel: true if input reads are single-ended, false if interleaved
    main:
        // Split single-end value channel into two branches, one of which will be empty
        single_end_check = single_end.branch{
            single: it
            paired: !it
        }
        // Forward reads into one of two channels based on endedness (the other will be empty)
        reads_ch_single = single_end_check.single.combine(reads_ch).map{it -> [it[1], it[2]] }
        reads_ch_paired = single_end_check.paired.combine(reads_ch).map{it -> [it[1], it[2]] }
        // In paired-end case, merge and join
        merged_ch = BBMERGE(reads_ch_paired)
        single_read_ch_paired = JOIN_FASTQ(merged_ch.reads, false).reads
        summarize_bbmerge_ch_paired = SUMMARIZE_BBMERGE(merged_ch.reads).summary
        // In single-end case, take unmodified reads
        single_read_ch_single = reads_ch_single
        summarize_bbmerge_ch_single = Channel.empty()
        single_read_ch = single_read_ch_paired.mix(single_read_ch_single)
        summarize_bbmerge_ch = summarize_bbmerge_ch_paired.mix(summarize_bbmerge_ch_single)
    emit:
        input_reads = reads_ch
        single_reads = single_read_ch
        bbmerge_summary = summarize_bbmerge_ch
}
