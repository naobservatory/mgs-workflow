/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_SINGLE_TARGET as SUBSET_SINGLE } from "../../../modules/local/subsetReads"
include { SUBSET_READS_PAIRED_TARGET as SUBSET_PAIRED } from "../../../modules/local/subsetReads"
include { FASTP } from "../../../modules/local/fastp"
include { FILTLONG } from "../../../modules/local/filtlong"
include { MINIMAP2_ONT as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { SAMTOOLS_FILTER } from "../../../modules/local/samtools"
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
      human_read_filtering
      random_seed
    main:
        if (single_end) {
            subset_ch = SUBSET_SINGLE(reads_ch, n_reads, "fastq", random_seed)
            inter_ch  = subset_ch
        } else {
            subset_ch = SUBSET_PAIRED(reads_ch, n_reads, "fastq", random_seed)
            inter_ch  = INTERLEAVE_FASTQ(subset_ch).output
        }
        if (ont) {
            cleaned_ch = FILTLONG(inter_ch)
            if (human_read_filtering) {
                minimap2_human_index = "s3://nao-mgs-simon/ont-indices/2024-12-14/minimap2-human-index/chm13v2.0.mmi"
                minimap2_ch = MINIMAP2_HUMAN(cleaned_ch, minimap2_human_index, "human")
                cleaned_ch = SAMTOOLS_FILTER(minimap2_ch.sam, "no-human")
            }
        } else {
            cleaned_ch = FASTP(inter_ch, adapter_path, !single_end)
        }
    emit:
        subset_reads = inter_ch
        trimmed_subset_reads = cleaned_ch.reads
        test_failed = ont ? Channel.empty() : cleaned_ch.failed // TODO: Capture rejected ONT reads somehow
}
