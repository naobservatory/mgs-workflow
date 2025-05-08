/*
Split a virus hit TSV by assigned taxonomic species, inferred from assigned
taxid by traversing the viral taxonomy tree.

Returns split hits in TSV and FASTQ formats, each as a flattened channel of
2-tuples ["group_species", FILE]
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { REHEAD_TSV } from "../../../modules/local/reheadTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { PARTITION_TSV } from "../../../modules/local/partitionTsv"
include { SORT_TSV as SORT_GROUP_TAXID } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_DB_TAXID } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_JOINED_SPECIES } from "../../../modules/local/sortTsv"
include { EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF_LABELED as EXTRACT_FASTQ } from "../../../modules/local/extractViralHitsToFastqNoref"

/***********
| WORKFLOW |
***********/

workflow SPLIT_VIRAL_TSV_BY_SPECIES {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        db // Viral taxonomy DB
        single_end
    main:
        // 1. Join to viral taxonomy DB
        rehead_ch = REHEAD_TSV(groups, "aligner_taxid", "taxid").output
        rehead_sorted_ch = SORT_GROUP_TAXID(rehead_ch, "taxid").sorted
        db_sorted_ch = SORT_DB_TAXID(
            Channel.of("db").combine(db),
            "taxid"
        ).sorted.map{ group, db -> db }
        join_prep_ch = rehead_sorted_ch.combine(db_sorted_ch)
        join_ch = JOIN_TSVS(join_prep_ch, "taxid", "left", "taxonomy").output
        // 2. Partition on species taxid
        join_sorted_ch = SORT_JOINED_SPECIES(join_ch, "taxid_species").sorted
        part_ch = PARTITION_TSV(join_sorted_ch, "taxid_species").output
        // 3. Restructure channel to separate species
        // First rearrange each element from [group, [paths]] to [[group_species1, path1], [group_species2, path2], ...]
        partitioned_flattened_ch = part_ch.flatMap{
            group, filepaths ->
                // Make sure paths are a list even if there's just one
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    // Define regex pattern for extracting species taxid
                    def filename = path.last()
                    def pattern = /^partition_(.*?)_sorted_taxid_species_${group}_taxonomy_left_joined_taxid\.tsv\.gz$/
                    def matcher = (filename =~ pattern)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${group}, ${path}, ${path.last()}"
                        throw new IllegalArgumentException(msg)
                    }
                    // Convert to 2-tuples with key combining group and species
                    ["${group}_${matcher[0][1]}", path]
                }
            }
        // 4. Extract into interleaved FASTQ format
        fastq_ch = EXTRACT_FASTQ(partitioned_flattened_ch, false, single_end).output
    emit:
        tsv = partitioned_flattened_ch
        fastq = fastq_ch
        test_in   = groups
        test_db   = db
        test_sort = rehead_sorted_ch
        test_join = join_ch
        test_part = part_ch
}
