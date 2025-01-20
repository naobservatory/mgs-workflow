/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_SINGLE_TARGET as SUBSET_SINGLE } from "../../../modules/local/subsetReads"
include { SUBSET_READS_PAIRED_TARGET as SUBSET_PAIRED } from "../../../modules/local/subsetReads"
include { FASTP_STREAMED as FASTP } from "../../../modules/local/fastp"
include { INTERLEAVE_FASTQ } from "../../../modules/local/interleaveFastq"

/***********
| WORKFLOW |
***********/

workflow SUBSET_TRIM_STREAMED {
    take:
      reads_ch
      n_reads
      adapter_path
      single_end
    main:
        if (single_end) {
            // TODO: Consider using output of COUNT_READS rather than re-counting during subsetting (could save a few minutes of clock-time)
            subset_ch = SUBSET_SINGLE(reads_ch, n_reads, "fastq")
            inter_ch  = subset_ch
            fastp_ch  = FASTP(subset_ch, adapter_path, false)
        } else {
            subset_ch = SUBSET_PAIRED(reads_ch, n_reads, "fastq")
            inter_ch  = INTERLEAVE_FASTQ(subset_ch).output
            fastp_ch  = FASTP(inter_ch, adapter_path, true)
        }
    emit:
        subset_reads = inter_ch
        trimmed_subset_reads = fastp_ch.reads
        test_failed = fastp_ch.failed
}
