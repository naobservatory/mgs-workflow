/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { FILTLONG } from "../../../modules/local/filtlong"
include { MASK_FASTQ_READS } from "../../../modules/local/maskRead"
include { PROCESS_VIRAL_MINIMAP2_SAM_LCA as PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2SamLca"
include { EXTRACT_SHARED_FASTQ_READS } from "../../../modules/local/extractSharedFastq"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { LCA_TSV } from "../../../modules/local/lcaTsv"
include { SORT_TSV as SORT_MINIMAP2_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_LCA } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { FILTER_TSV_COLUMN_BY_VALUE } from "../../../modules/local/filterTsvColumnByValue"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_ONT_LCA {
    take:
        reads_ch
        ref_dir
        taxid_artificial
    main:
        // Get reference_paths
        minimap2_virus_index = "${ref_dir}/results/mm2-virus-index"
        minimap2_human_index = "${ref_dir}/results/mm2-human-index"
        minimap2_contam_index = "${ref_dir}/results/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        nodes_db = "${ref_dir}/results/taxonomy-nodes.dmp"
        names_db = "${ref_dir}/results/taxonomy-names.dmp"
        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 50, 15000, 90)
        // Mask non-complex read sections
        masked_ch = MASK_FASTQ_READS(filtered_ch, 25, 0.55)
        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false, "")
        no_human_ch = human_minimap2_ch.reads_unmapped
        // Identify other contaminants
        contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false, "")
        no_contam_ch = contam_minimap2_ch.reads_unmapped
        // Identify virus reads with multiple alignments for LCA analysis
        virus_minimap2_ch = MINIMAP2_VIRUS(no_contam_ch, minimap2_virus_index, "virus", false, "-N 10")
        virus_sam_ch = virus_minimap2_ch.sam
        // Group cleaned reads and sam files by sample
        sam_fastq_ch = virus_sam_ch.join(filtered_ch)
        // Generate TSV of viral hits, and sort
        processed_minimap2_ch = PROCESS_VIRAL_MINIMAP2_SAM(sam_fastq_ch, genome_meta_path, virus_db_path)
        processed_minimap2_sorted_ch = SORT_MINIMAP2_VIRAL(processed_minimap2_ch.output, "seq_id")
        // Run LCA on viral hits TSV
        lca_ch = LCA_TSV(processed_minimap2_sorted_ch.sorted, nodes_db, names_db, "seq_id", 
            "aligner_taxid", "aligner_length_normalized_score", taxid_artificial, "aligner")
        // Sort the LCA TSV files by seq_id for joining
        lca_sorted_ch = SORT_LCA(lca_ch.output, "seq_id")
        // Join LCA and Minimap2 TSV on seq_id to combine taxonomic and alignment data
        joined_input_ch = lca_sorted_ch.sorted.join(processed_minimap2_sorted_ch.sorted, by: 0)
        joined_ch = JOIN_TSVS(joined_input_ch, "seq_id", "inner", "lca_minimap2")
        // Filter to keep only primary alignments (aligner_secondary_status=False)
        primary_lca_ch = FILTER_TSV_COLUMN_BY_VALUE(joined_ch.output, "aligner_secondary_status", false, true)
        tsv_labeled_ch = ADD_SAMPLE_COLUMN(primary_lca_ch.output, "sample", "viral_minimap2") 
        // Concatenate TSVs of viral hits
        viral_tsvs = tsv_labeled_ch.output.map { it[1] }.collect()
        merged_tsv_ch = CONCATENATE_TSVS(viral_tsvs, "virus_hits_final") 
        // Pull out clean reads from mapped reads to feed into BLAST
        virus_fastq_ch = virus_minimap2_ch.reads_mapped
        clean_virus_fastq_ch = EXTRACT_SHARED_FASTQ_READS(virus_fastq_ch.join(filtered_ch.reads))
        fastq_ch = CONCATENATE_FILES(clean_virus_fastq_ch.output.map{ it[1] }.collect(), "virus_hits_final", "fastq.gz")
    emit:
        hits_final = merged_tsv_ch.output
        hits_fastq = fastq_ch.output
        test_minimap2_virus = virus_sam_ch
        test_fastq_filtered_human = human_minimap2_ch.reads_unmapped
        test_fastq_filtered_contam = contam_minimap2_ch.reads_unmapped
}
