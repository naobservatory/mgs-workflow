/*
Perform efficient post-hoc validation of putative viral reads identified by the RUN workflow.

A. Partition putative hits by assigned species
B. Extract into FASTQ and merge read pairs into single sequences
C. [TODO] Cluster sequences from each species and identify representative sequences
D. [TODO] Align representative sequences against a large reference DB
E. [TODO] Compare taxids assigned in (4) to those assigned by RUN workflow
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
include { BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ } from "../../../modules/local/joinFastq"
include { VSEARCH_CLUSTER } from "../../../modules/local/vsearch"

/***********
| WORKFLOW |
***********/

workflow VALIDATE_VIRAL_ASSIGNMENTS {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        db // Viral taxonomy DB
        cluster_identity // Identity threshold for VSEARCH clustering
        cluster_min_len // Minimum sequence length for VSEARCH clustering
    main:
        // 1. Join to viral taxonomy DB
        rehead_ch = REHEAD_TSV(groups, "bowtie2_taxid_best", "taxid").output
        rehead_sorted_ch = SORT_GROUP_TAXID(rehead_ch, "taxid").sorted
        db_sorted_ch = SORT_DB_TAXID(
            Channel.of("db").combine(db),
            "taxid"
        ).sorted.map{ sample, db -> db }
        join_prep_ch = rehead_sorted_ch.combine(db_sorted_ch)
        join_ch = JOIN_TSVS(join_prep_ch, "taxid", "left", "taxonomy").output
        // 2. Partition on species taxid
        join_sorted_ch = SORT_JOINED_SPECIES(join_ch, "taxid_species").sorted
        part_ch = PARTITION_TSV(join_sorted_ch, "taxid_species").output
        // 3. Restructure channel to separate species
        // First rearrange each element from [sample, [paths]] to [[sample_species1, path1], [sample_species2, path2], ...]
        partitioned_flattened_ch = part_ch.flatMap{
            sample, filepaths ->
                // Make sure paths are a list even if there's just one
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    // Define regex pattern for extracting species taxid
                    def filename = path.last()
                    def pattern = /^partition_(.*?)_sorted_taxid_species_${sample}_taxonomy_left_joined_taxid\.tsv\.gz$/
                    def matcher = (filename =~ pattern)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${sample}, ${path}, ${path.last()}"
                        throw new IllegalArgumentException(msg)
                    }
                    // Convert to 2-tuples with key combining sample and species
                    ["${sample}_${matcher[0][1]}", path]
                }
            }
        // 4. Extract into interleaved FASTQ format
        fastq_ch = EXTRACT_FASTQ(partitioned_flattened_ch, false).output
        // 5. Merge and concatenate pairs into single sequences for clustering
        bbmerge_ch = BBMERGE(fastq_ch).reads
        concat_ch = JOIN_FASTQ(bbmerge_ch, false).reads
        // 6. Cluster merged reads
        cluster_ch = VSEARCH_CLUSTER(concat_ch, cluster_identity, 0, cluster_min_len)
        // TODO: Implement C,D,E from docstring
    emit:
        test_in   = groups
        test_db   = db
        test_sort = rehead_sorted_ch
        test_join = join_ch
        test_part = part_ch
        test_part_2 = partitioned_flattened_ch
        test_fastq = fastq_ch
        test_concat = concat_ch
        test_cluster_reps = cluster_ch.reps
        test_cluster_summ = cluster_ch.summary
}
