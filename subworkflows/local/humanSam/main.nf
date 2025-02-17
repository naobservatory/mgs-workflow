include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2_ONT as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { SAMTOOLS_KEEP_AS_SAM } from "../../../modules/local/samtools"
include { MERGE_SAM } from "../../../modules/local/samtools"

workflow HUMAN_SAM {
    take:
        reads_ch
        minimap2_human_index
    main:
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(reads_ch)

        // Identify HV reads
        minimap2_human_sam_ch = MINIMAP2_HUMAN(masked_ch, minimap2_human_index, "human")

        // Drop reads that didn't align to human
        human_sam_ch = SAMTOOLS_KEEP_AS_SAM(minimap2_human_sam_ch, "human")

        sam_files_ch = human_sam_ch.sam.collect()

        // Merge SAM files
        merged_sam_ch = MERGE_SAM(sam_files_ch, "human")
    emit:
        sam = merged_sam_ch.merged_sam
}
