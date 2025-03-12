/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2_NON_STREAMED as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { CONCATENATE_TSVS as CONCATENATE_HV_TSVS } from "../../../modules/local/concatenateTsvs"
include { ADD_SAMPLE_COLUMN as LABEL_HV_TSVS } from "../../../modules/local/addSampleColumn"
include { FILTLONG } from "../../../modules/local/filtlong"
include { MASK_FASTQ_READS } from "../../../modules/local/maskRead"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
include { BLAST_VIRAL } from "../../../subworkflows/local/blastViral"
include { EXTRACT_SHARED_FASTQ_READS } from "../../../modules/local/extractSharedFastq"

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
        minimap2_virus_index = "${ref_dir}/results/mm2-virus-index"
        minimap2_human_index = "${ref_dir}/results/mm2-human-index"
        minimap2_contam_index = "${ref_dir}/results/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"


        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 50, 90)

        // Mask non-complex read sections
        masked_ch = MASK_FASTQ_READS(filtered_ch, 25, 0.55)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
        contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
        no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_human_ch, minimap2_virus_index, "hv", false)

        virus_sam_ch = virus_minimap2_ch.sam
        virus_fastq_ch = virus_minimap2_ch.reads_mapped

        // Pull out clean reads from mapped reads
        clean_matched_subset_ch = EXTRACT_SHARED_FASTQ_READS(virus_fastq_ch.join(filtered_ch.reads))

        // Create common channel for HV SAM and clean HV reads
        sam_and_reads_ch = virus_sam_ch.join(clean_matched_subset_ch.output)

        // Make clean reads ready for BLAST
        clean_reads_ch = clean_matched_subset_ch.output.map { it[1] }
        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(sam_and_reads_ch, genome_meta_path, virus_db_path, host_taxon)
        hv_tsv_labeled_ch = LABEL_HV_TSVS(hv_tsv_ch.output, "sample", "hv_tsv")

        // Concatenate HV TSVs
        hv_tsvs = hv_tsv_ch.output.map { it[1] }.collect()
        merged_tsv_ch = CONCATENATE_HV_TSVS(hv_tsvs, "hv")

    emit:
        hv_tsv = merged_tsv_ch.output
        hv_fastq = clean_matched_subset_ch.output
        test_minimap2_virus = virus_sam_ch
        test_fastq_filtered_human = human_minimap2_ch.reads_unmapped
        test_fastq_filtered_contam = contam_minimap2_ch.reads_unmapped
}
