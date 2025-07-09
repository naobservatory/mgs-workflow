// Process LCA and aligner columns by sorting, joining, filtering, selecting, and renaming columns
/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_LCA } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_ALIGNER } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { FILTER_TSV_COLUMN_BY_VALUE } from "../../../modules/local/filterTsvColumnByValue"
include { SELECT_TSV_COLUMNS } from "../../../modules/local/selectTsvColumns"
include { REHEAD_TSV } from "../../../modules/local/reheadTsv"

/***********
| WORKFLOW |
***********/

workflow PROCESS_LCA_ALIGNER_OUTPUT {
    take:
        lca_tsv        // Output from LCA_TSV process
        aligner_tsv    // Output from PROCESS_VIRAL_{MINIMAP2,BOWTIE2}_SAM
        column_prefix  // Prefix to add to specified columns
    main:
       // Define columns to keep, separating by ones to prefix and ones to not
        col_keep_no_prefix = ["seq_id", "aligner_taxid_lca", "aligner_taxid_top", 
                              "aligner_length_normalized_score_mean", "aligner_taxid_lca_natural",
                              "aligner_n_assignments_natural", "aligner_length_normalized_score_mean_natural",
                              "aligner_taxid_lca_artificial", "aligner_n_assignments_artificial", 
                              "aligner_length_normalized_score_mean_artificial"]
        col_keep_add_prefix = ["genome_id_all", "taxid_all", "fragment_length", 
                               "best_alignment_score", "best_alignment_score_rev",
                               "edit_distance", "edit_distance_rev", "ref_start", 
                               "ref_start_rev", "query_len", "query_len_rev",
                               "query_seq", "query_seq_rev", "query_rc", 
                               "query_rc_rev", "query_qual", "query_qual_rev", 
                               "pair_status"]
        // Step 1: Sort both TSV files by seq_id for joining
        lca_sorted_ch = SORT_LCA(lca_tsv, "seq_id")
        aligner_sorted_ch = SORT_ALIGNER(aligner_tsv, "seq_id")
        // Step 2: Join LCA and aligner TSV on seq_id
        joined_input_ch = lca_sorted_ch.sorted.join(aligner_sorted_ch.sorted, by: 0)
        joined_ch = JOIN_TSVS(joined_input_ch, "seq_id", "inner", "lca_aligner")
        // Step 3: Filter to keep only primary alignments
        filtered_ch = FILTER_TSV_COLUMN_BY_VALUE(joined_ch.output, "classification", "primary", true)
        // Step 4: Select specific columns
        col_keep = (col_keep_no_prefix + col_keep_add_prefix).join(",")
        selected_ch = SELECT_TSV_COLUMNS(filtered_ch.output, col_keep, "keep")
        // Step 5: Rename columns with prefix using REHEAD_TSV
        old_cols = col_keep_add_prefix.join(",")
        new_cols = col_keep_add_prefix.collect { "${column_prefix}${it}" }.join(",")
        renamed_ch = REHEAD_TSV(selected_ch.output, old_cols, new_cols)
    emit:
        output = renamed_ch.output
        test_lca = lca_tsv
        test_aligner = aligner_tsv
        test_col_add_prefix = col_keep_add_prefix
        test_col_no_prefix = col_keep_no_prefix
}
