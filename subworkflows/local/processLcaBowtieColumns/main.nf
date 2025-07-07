// Process LCA and Bowtie2 columns by sorting, joining, filtering, selecting, and renaming columns
/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_LCA } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BOWTIE2 } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { FILTER_TSV_COLUMN_BY_VALUE } from "../../../modules/local/filterTsvColumnByValue"
include { SELECT_TSV_COLUMNS } from "../../../modules/local/selectTsvColumns"
include { REHEAD_TSV } from "../../../modules/local/reheadTsv"

/***********
| WORKFLOW |
***********/

workflow PROCESS_LCA_BOWTIE_COLUMNS {
    take:
        lca_tsv        // Output from LCA_TSV process
        bowtie2_tsv    // Output from PROCESS_VIRAL_BOWTIE2_SAM
        col_keep_no_prefix   // List of columns to keep without prefix
        col_keep_add_prefix  // List of columns to keep with prefix
    main:
        // Step 1: Sort both TSV files by seq_id for joining
        lca_sorted_ch = SORT_LCA(lca_tsv, "seq_id")
        bowtie2_sorted_ch = SORT_BOWTIE2(bowtie2_tsv, "seq_id")
        // Step 2: Join LCA and Bowtie2 TSV on seq_id
        joined_input_ch = lca_sorted_ch.sorted.join(bowtie2_sorted_ch.sorted, by: 0)
        joined_ch = JOIN_TSVS(joined_input_ch, "seq_id", "inner", "lca_bowtie2")
        // Step 3: Filter to keep only primary alignments
        filtered_ch = FILTER_TSV_COLUMN_BY_VALUE(joined_ch.output, "classification", "primary", true)
        // Step 4: Select specific columns
        col_keep = (col_keep_no_prefix + col_keep_add_prefix).join(",")
        selected_ch = SELECT_TSV_COLUMNS(filtered_ch.output, col_keep, "keep")
        // Step 5: Rename columns with prefix using REHEAD_TSV
        old_cols = col_keep_add_prefix.join(",")
        new_cols = col_keep_add_prefix.collect { "prim_align_${it}" }.join(",")
        renamed_ch = REHEAD_TSV(selected_ch.output, old_cols, new_cols)
    emit:
        output = renamed_ch.output
}
