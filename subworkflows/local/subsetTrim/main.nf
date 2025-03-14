/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_SINGLE_TARGET as SUBSET_SINGLE } from "../../../modules/local/subsetReads"
include { SUBSET_READS_PAIRED_TARGET as SUBSET_PAIRED } from "../../../modules/local/subsetReads"
include { FASTP } from "../../../modules/local/fastp"
include { FILTLONG } from "../../../modules/local/filtlong"
include { INTERLEAVE_FASTQ } from "../../../modules/local/interleaveFastq"

/***********
| WORKFLOW |
***********/

workflow SUBSET_TRIM {
    take:
      reads_ch
      n_reads
      adapter_path
      single_end
      ont
      random_seed
    main:
        // Split single-end value channel into two branches, one of which will be empty
        single_end_check = single_end.branch{
            single: it
            paired: !it
        }
        // Forward reads into one of two channels based on endedness (the other will be empty)
        reads_ch_single = single_end_check.single.combine(reads_ch).map{it -> [it[1], it[2]] }
        reads_ch_paired = single_end_check.paired.combine(reads_ch).map{it -> [it[1], it[2]] }
        // Subset reads according to endedness (other channel will be empty)
        subset_ch_single = SUBSET_SINGLE(reads_ch_single, n_reads, "fastq", random_seed)
        subset_ch_paired = SUBSET_PAIRED(reads_ch_paired, n_reads, "fastq", random_seed)
        // Interleave reads based on endedness (other channel will be empty)
        inter_ch_single = subset_ch_single
        inter_ch_paired = INTERLEAVE_FASTQ(subset_ch_paired).output
        inter_ch = inter_ch_single.mix(inter_ch_paired)
        // Read cleaning
        if (ont) {
            cleaned_ch = FILTLONG(inter_ch, 100, 15000, 90)
        } else {
            cleaned_ch = FASTP(inter_ch, adapter_path, single_end.map{!it})
        }
    emit:
        subset_reads = inter_ch
        trimmed_subset_reads = cleaned_ch.reads
        test_failed = ont ? Channel.empty() : cleaned_ch.failed // TODO: Capture rejected ONT reads somehow
}
