/*
Given three TSVs...
1. A viral hits TSV containing preliminary taxid assignments for each sequence
2. A cluster TSV mapping each sequence to its cluster representative
3. A validation TSV containing validation statuses for cluster representatives
...merge them together and propagate validation information from cluster representatives
to individual hits.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_CLUSTER_TSV } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_VALIDATION_TSV } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_JOIN_TSV } from "../../../modules/local/sortTsv"
include { JOIN_TSVS as JOIN_CLUSTER_VALIDATION } from "../../../modules/local/joinTsvs"
include { JOIN_TSVS as JOIN_HITS_CLUSTER } from "../../../modules/local/joinTsvs"
include { SELECT_TSV_COLUMNS as DROP_TAXID } from "../../../modules/local/selectTsvColumns"
include { SELECT_TSV_COLUMNS as DROP_GROUP } from "../../../modules/local/selectTsvColumns"

/***********
| WORKFLOW |
***********/

workflow PROPAGATE_VALIDATION_INFORMATION {
    take:
        hits_tsv // Viral hits TSV, including seq_id and original taxid assignment
        cluster_tsv // Cluster TSV, including seq_id and cluster_rep_id
        validation_tsv // Validation TSV, including cluster_rep_id and validation status
        taxid_column // Column header for taxid in hits TSV
    main:
        // 0. Drop redundant columns from validation TSV
        drop_taxid_ch = DROP_TAXID(validation_tsv, taxid_column, "drop").output
        // 1. Sort cluster and validation TSVs by cluster_rep_id
        sort_cluster_tsv = SORT_CLUSTER_TSV(cluster_tsv, "vsearch_cluster_rep_id").sorted
        sort_validation_tsv = SORT_VALIDATION_TSV(drop_taxid_ch, "vsearch_cluster_rep_id").sorted
        // 2. Left-join cluster and validation TSVs by cluster_rep_id
        combine_1_ch = sort_cluster_tsv.combine(sort_validation_tsv, by: 0)
        join_1_ch = JOIN_CLUSTER_VALIDATION(combine_1_ch, "vsearch_cluster_rep_id", "left", "representative").output
        // 3. Sort intermediate TSV by seq_id (hits TSV should already be sorted)
        sort_join_1_tsv = SORT_JOIN_TSV(join_1_ch, "seq_id").sorted
        // 4. Strict-join intermediate and hits TSVs by seq_id
        // Both tables should contain the same sequences in the same order
        drop_group_ch = DROP_GROUP(sort_join_1_tsv, "group", "drop").output
        combine_2_ch = hits_tsv.combine(drop_group_ch, by: 0)
        join_2_ch = JOIN_HITS_CLUSTER(combine_2_ch, "seq_id", "strict", "validation").output
    emit:
        output = join_2_ch
        test_intermediate = join_1_ch
        test_drop_taxid = drop_taxid_ch
        test_drop_group = drop_group_ch
}
