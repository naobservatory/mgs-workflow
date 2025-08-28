/*
Split a virus hit TSV by assigned taxonomic species, inferred from assigned
taxid by traversing the viral taxonomy tree.

Returns split hits in TSV and FASTQ formats, each as a flattened channel of
2-tuples ["group_species", FILE]
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_SEQ_ID } from "../../../modules/local/sortTsv"
include { CHECK_TSV_DUPLICATES } from "../../../modules/local/checkTsvDuplicates"
include { SELECT_TSV_COLUMNS } from "../../../modules/local/selectTsvColumns"
include { REHEAD_TSV } from "../../../modules/local/reheadTsv"
include { SORT_TSV as SORT_DB_TAXID } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { PARTITION_TSV } from "../../../modules/local/partitionTsv"
include { SORT_TSV as SORT_GROUP_TAXID } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_JOINED_SPECIES } from "../../../modules/local/sortTsv"
include { EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF_LABELED as EXTRACT_FASTQ } from "../../../modules/local/extractViralHitsToFastqNoref"
include { ADD_CONDITIONAL_TSV_COLUMN } from "../../../modules/local/addConditionalTsvColumn"

/***********
| WORKFLOW |
***********/

workflow SPLIT_VIRAL_TSV_BY_SELECTED_TAXID {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        db // Viral taxonomy DB
    main:
       // 0. Check for duplicate read IDs (throw error if found)
        sort_seq_id_ch = SORT_SEQ_ID(groups, "seq_id").sorted
        check_ch = CHECK_TSV_DUPLICATES(sort_seq_id_ch, "seq_id").output
        // 1. Prepare taxonomy DB for joining
        db_ch = Channel.of("db").combine(db)
        db_select_ch = SELECT_TSV_COLUMNS(db_ch, "taxid,taxid_species", "keep").output
        db_rehead_ch = REHEAD_TSV(db_select_ch, "taxid", "aligner_taxid_lca").output
        db_sorted_ch = SORT_DB_TAXID(db_rehead_ch, "aligner_taxid_lca").sorted
        db_raw_ch = db_sorted_ch.map{ _db_label, db_contents -> db_contents }
        // 2. Join hits TSVs to prepared taxonomy DB
        sorted_ch = SORT_GROUP_TAXID(check_ch, "aligner_taxid_lca").sorted
        join_prep_ch = sorted_ch.combine(db_raw_ch)
        join_ch = JOIN_TSVS(join_prep_ch, "aligner_taxid_lca", "left", "taxonomy").output
        // 3. Partition on species taxid        
        updated_col = ADD_CONDITIONAL_TSV_COLUMN(join_ch, [
            chk_col: "taxid_species",
            match_val: "NA",
            if_col: "aligner_taxid_lca",
            else_col: "taxid_species",
            new_hdr: "selected_taxid"
        ]).tsv
        join_sorted_ch = SORT_JOINED_SPECIES(updated_col, "selected_taxid").sorted
        part_ch = PARTITION_TSV(join_sorted_ch, "selected_taxid").output
        // 4. Restructure channel to separate species
        partitioned_flattened_ch = part_ch.flatMap{
            group, filepaths ->
                // Make sure paths are a list even if there's just one
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    // Define regex pattern for extracting species taxid
                    def filename = path.last()
                    def pattern = "^partition_(.*?)_sorted_selected_taxid_added_selected_taxid_${group}_taxonomy_left_joined_aligner_taxid_lca\\.tsv\\.gz\$"
                    def matcher = (filename =~ pattern)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${group}, ${path}, ${path.last()}"
                        throw new IllegalArgumentException(msg)
                    }
                    // Convert to 2-tuples with key combining group and species
                    ["${group}_${matcher[0][1]}", path]
                }
            }
        // 5. Extract into interleaved FASTQ format
        fastq_ch = EXTRACT_FASTQ(partitioned_flattened_ch, false).output
    emit:
        tsv = partitioned_flattened_ch
        fastq = fastq_ch
        test_in   = groups
        test_db   = db
        test_sort = sorted_ch
        test_join = join_ch
        test_part = part_ch
}
