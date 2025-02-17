include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2_ONT as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MERGE_SAM } from "../../../modules/local/samtools"

workflow HUMAN_SAM {
    take:
        reads_ch
    main:
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(reads_ch)

        // Identify HV reads
        minimap2_human_index = "s3://nao-mgs-simon/ont-indices/2024-12-14/minimap2-human-index/chm13v2.0.mmi"
        minimap2_human_sam_ch = MINIMAP2_HUMAN(masked_ch, minimap2_human_index, "human")

        // Extract just the SAM file path from the tuple
        sam_files_ch = minimap2_human_sam_ch.sam.map { it[1] }

        // Merging doesn't yet work, to fix.
        merged_sam_ch = MERGE_SAM(sam_files_ch.collect(), "human")
    emit:
        sam = merged_sam_ch.merged_sam
}