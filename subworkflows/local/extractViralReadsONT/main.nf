/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { MERGE_SAM } from "../../../modules/local/samtools"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_ONT {
    take:
        reads_ch
        ref_dir
        host_taxon
    main:
        // Get reference_paths
        // TODO: RENAME HV INDEX TO VIRUS INDEX
        minimap2_hv_index = "${projectDir}/test-index/mm2-virus-index"
        minimap2_human_index = "${projectDir}/test-index/mm2-human-index"
        minimap2_contam_index = "${projectDir}/test-index/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"

        // Drop non-complex reads
        // masked_ch = DUSTMASKER_FASTQ_GZIPPED(reads_ch)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(reads_ch, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
        contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
        no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_contam_ch, minimap2_hv_index, "virus", false)

        virus_sam_ch = virus_minimap2_ch.sam.map { it[1] }

        // Merge SAM files
        merged_sam_ch = MERGE_SAM(virus_sam_ch, "hv")

        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(merged_sam_ch, genome_meta_path, virus_db_path, host_taxon)

    emit:
        hv_tsv = hv_tsv_ch.output
}