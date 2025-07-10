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
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS as CONCATENATE_ALIGNER } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_LCA } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_VIRAL_HITS } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow PROCESS_LCA_ALIGNER_OUTPUT {
    take:
        lca_tsv        // Output from LCA_TSV process
        aligner_tsv    // Output from PROCESS_VIRAL_{MINIMAP2,BOWTIE2}_SAM
        col_keep_no_prefix // Columns to keep without prefix
        col_keep_add_prefix // Columns to keep with prefix
        column_prefix  // Prefix to add to specified columns
    main:
        // Step 1: Sort LCA tsv by seq_id (aligner_tsv is already sorted)
        lca_sorted_ch = SORT_LCA(lca_tsv, "seq_id")
        // Step 2: Join LCA and aligner TSV on seq_id
        joined_input_ch = lca_sorted_ch.sorted.join(aligner_tsv, by: 0)
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

        // Combine bowtie2 output
        aligner_combined_ch = aligner_tsv.map{ _sample, file -> file }.collect().ifEmpty([])
        aligner_concat_ch = CONCATENATE_ALIGNER(aligner_combined_ch, "aligner_final")

        // Combine LCA output
        lca_labeled_ch = ADD_SAMPLE_COLUMN(lca_tsv, "sample", "viral_lca")
        lca_combined_ch = lca_labeled_ch.output.map{ _sample, file -> file }.collect().ifEmpty([])
        lca_concat_ch = CONCATENATE_LCA(lca_combined_ch, "lca_final")

        // Combine final output
        label_combined_ch = renamed_ch.output.map{ _sample, file -> file }.collect().ifEmpty([])
        concat_ch = CONCATENATE_VIRAL_HITS(label_combined_ch, "virus_hits_final")

    emit:
        viral_hits_tsv = concat_ch.output
        aligner_tsv = aligner_concat_ch.output
        lca_tsv = lca_concat_ch.output
}
