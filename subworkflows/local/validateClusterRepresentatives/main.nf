/*
Given a viral hits TSV containing representative sequences,
and tabular validation information for each representative,
merge the two tables together on sequence ID, then compare
the two taxids to evaluate validation status.

The input map distance_params should specify the following fields
for computing taxonomic distance between the original and validated taxids:
- taxid_field_1: Column header for original taxid
- taxid_field_2: Column header for validated taxid
- distance_field_1: Column header for original taxid distance
- distance_field_2: Column header for validated taxid distance
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SELECT_TSV_COLUMNS } from "../../../modules/local/selectTsvColumns"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { REHEAD_TSV as REHEAD_QSEQID} from "../../../modules/local/reheadTsv"
include { COMPUTE_TAXID_DISTANCE } from "../../../modules/local/computeTaxidDistance"
include { REHEAD_TSV as REHEAD_SEQ_ID } from "../../../modules/local/reheadTsv"

/***********
| WORKFLOW |
***********/

workflow VALIDATE_CLUSTER_REPRESENTATIVES {
    take:
        hits_tsv // Viral hits TSV, including seq_id and original taxid assignment
        lca_tsv // LCA TSV from validation against core_nt
        ref_dir // Path to reference directory containing taxonomy information
        distance_params // Map specifying input taxid fields and output distance fields
    main:
        // 0. Get reference paths
        nodes_db = "${ref_dir}/results/taxonomy-nodes.dmp"
        // 1. Prepare inputs for joining
        // Subset hits TSV to only seq_id and taxid columns
        def hits_taxid_column_str = "seq_id,${distance_params.taxid_field_1}"
        select_ch = SELECT_TSV_COLUMNS(hits_tsv, hits_taxid_column_str, "keep").output
        // Rename qseqid to seq_id in LCA TSV
        rehead_ch = REHEAD_QSEQID(lca_tsv, "qseqid", "seq_id").output
        // Combine channels for joining
        combine_ch = select_ch.combine(rehead_ch, by: 0)
        // 2. Inner-join hits and LCA tables by seq_id
        // Inner join subsets hits TSV to sequences in LCA TSV (i.e. cluster representatives)
        join_ch = JOIN_TSVS(combine_ch, "seq_id", "inner", "representative").output
        // 3. Compute taxonomic distance between original and validated taxids
        dist_ch = COMPUTE_TAXID_DISTANCE(join_ch, distance_params, nodes_db).output
        // NB: As implemented, this will produce a negative distance if the original
        // taxid is too high (ancestor of LCA taxid), and a positive distance if the
        // original taxid is too low (descendant of LCA taxid).
        // 4. Rename seq_id to cluster_rep_id for downstream processing
        rename_ch = REHEAD_SEQ_ID(dist_ch, "seq_id", "vsearch_cluster_rep_id").output
    emit:
        output = rename_ch
        test_dist = dist_ch
        test_select = select_ch
        test_rehead = rehead_ch
        test_join = join_ch
}
