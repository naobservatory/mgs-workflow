/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { MERGE_SAM as MERGE_SAM_VIRUS } from "../../../modules/local/samtools"
include { MERGE_SAM as MERGE_SAM_HUMAN } from "../../../modules/local/samtools"
include { FILTLONG } from "../../../modules/local/filtlong"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
include { BLAST_VIRAL } from "../../../subworkflows/local/blastViral"
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
        minimap2_virus_index = "${projectDir}/test-index/mm2-virus-index"
        minimap2_human_index = "${projectDir}/test-index/mm2-human-index"
        minimap2_contam_index = "${projectDir}/test-index/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"

        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 100, 10000, 90)

        // Mark non-complex read sections
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(filtered_ch)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
        contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
        no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_contam_ch, minimap2_virus_index, "virus_1", false)
        // Pull out files only
        virus_sam_ch = virus_minimap2_ch.sam.map { it[1] }

        // Merge SAM files
        merged_virus_sam_ch = MERGE_SAM_VIRUS(virus_sam_ch.collect(), "hv")

        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(merged_virus_sam_ch, genome_meta_path, virus_db_path, host_taxon)

    emit:
        hv_tsv = hv_tsv_ch.output
}