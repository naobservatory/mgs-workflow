include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MERGE_SAM } from "../../../modules/local/samtools"

workflow HUMAN_SAM {
    take:
        reads_ch
        minimap2_human_index
    main:
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(reads_ch)

        // Identify HV reads
        minimap2_human_sam_ch = MINIMAP2_HUMAN(masked_ch, minimap2_human_index, "human")

        // Merging doesn't yet work, to fix.
        merged_sam_ch = MERGE_SAM(minimap2_human_sam_ch.sam.collect(), "human")
    emit:
        human_sam = merged_sam_ch.merged_sam
}
