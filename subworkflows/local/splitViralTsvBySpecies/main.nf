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

/***********
| WORKFLOW |
***********/

workflow SPLIT_VIRAL_TSV_BY_SPECIES {
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
        db_rehead_ch = REHEAD_TSV(db_select_ch, "taxid", "bowtie2_taxid_best").output
        db_sorted_ch = SORT_DB_TAXID(db_rehead_ch, "bowtie2_taxid_best").sorted
        db_raw_ch = db_sorted_ch.map{ group, db -> db }
        // 2. Join hits TSVs to prepared taxonomy DB
        sorted_ch = SORT_GROUP_TAXID(check_ch, "bowtie2_taxid_best").sorted
        join_prep_ch = sorted_ch.combine(db_raw_ch)
        join_ch = JOIN_TSVS(join_prep_ch, "bowtie2_taxid_best", "left", "taxonomy").output
        // 3. Partition on species taxid
        join_sorted_ch = SORT_JOINED_SPECIES(join_ch, "taxid_species").sorted
        part_ch = PARTITION_TSV(join_sorted_ch, "taxid_species").output
        // 4. Restructure channel to separate species
        partitioned_flattened_ch = part_ch.flatMap{
            group, filepaths ->
                // Make sure paths are a list even if there's just one
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    // Define regex pattern for extracting species taxid
                    def filename = path.last()
                    def pattern = /^partition_(.*?)_sorted_taxid_species_${group}_taxonomy_left_joined_bowtie2_taxid_best\.tsv\.gz$/
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
