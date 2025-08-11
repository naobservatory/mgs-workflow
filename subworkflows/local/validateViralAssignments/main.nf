/*
Perform efficient post-hoc validation of putative viral reads identified by the RUN workflow.

A. Partition putative hits by assigned species
B. Extract into FASTQ and merge read pairs into single sequences
C. Cluster sequences from each species and identify representative sequences
D. Align representative sequences against a large reference DB
E. Compare taxids assigned in (4) to those assigned by RUN workflow
F. [TODO] Propagate validation information from cluster representatives to other hits
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SPLIT_VIRAL_TSV_BY_SELECTED_TAXID } from "../../../subworkflows/local/splitViralTsvBySelectedTaxid"
include { CLUSTER_VIRAL_ASSIGNMENTS } from "../../../subworkflows/local/clusterViralAssignments"
include { CONCATENATE_FILES_ACROSS_SELECTED_TAXID } from "../../../subworkflows/local/concatenateFilesAcrossSelectedTaxid"
include { CONCATENATE_TSVS_ACROSS_SELECTED_TAXID } from "../../../subworkflows/local/concatenateTsvsAcrossSelectedTaxid"
include { BLAST_FASTA } from "../../../subworkflows/local/blastFasta"
include { VALIDATE_CLUSTER_REPRESENTATIVES } from "../../../subworkflows/local/validateClusterRepresentatives"
include { PROPAGATE_VALIDATION_INFORMATION } from "../../../subworkflows/local/propagateValidationInformation"
include { SELECT_TSV_COLUMNS } from "../../../modules/local/selectTsvColumns"
include { COPY_FILE as COPY_HITS } from "../../../modules/local/copyFile"
include { COPY_FILE as COPY_BLAST } from "../../../modules/local/copyFile"

/***********
| WORKFLOW |
***********/

workflow VALIDATE_VIRAL_ASSIGNMENTS {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        db // Viral taxonomy DB
        params_map // Map containing all parameters
    main:
        // Extract parameters from map
        cluster_identity = params_map.cluster_identity
        cluster_min_len = params_map.cluster_min_len
        n_clusters = params_map.n_clusters
        ref_dir = params_map.ref_dir
        blast_db_prefix = params_map.blast_db_prefix
        perc_id = params_map.perc_id
        qcov_hsp_perc = params_map.qcov_hsp_perc
        blast_max_rank = params_map.blast_max_rank
        blast_min_frac = params_map.blast_min_frac
        taxid_artificial = params_map.taxid_artificial
        
        // 1. Split viral hits TSV by species
        split_ch = SPLIT_VIRAL_TSV_BY_SELECTED_TAXID(groups, db)
        // 2. Cluster sequences within species and obtain representatives of largest clusters
        cluster_params = [
            cluster_identity: cluster_identity,
            cluster_min_len: cluster_min_len,
            n_clusters: n_clusters
        ]
        cluster_ch = CLUSTER_VIRAL_ASSIGNMENTS(split_ch.fastq, cluster_params, Channel.of(false))
        // 3. Concatenate data across species (prepare for group-level BLAST)
        concat_fasta_ch = CONCATENATE_FILES_ACROSS_SELECTED_TAXID(cluster_ch.fasta, "cluster_reps")
        concat_cluster_ch = CONCATENATE_TSVS_ACROSS_SELECTED_TAXID(cluster_ch.tsv, "cluster_info")
        // 4. Run BLAST on concatenated cluster representatives (single job per group)
        blast_fasta_params = [
            ref_dir: ref_dir,
            blast_db_prefix: blast_db_prefix,
            perc_id: perc_id,
            qcov_hsp_perc: qcov_hsp_perc,
            blast_max_rank: blast_max_rank,
            blast_min_frac: blast_min_frac,
            taxid_artificial: taxid_artificial,
            lca_prefix: "validation"
        ]
        blast_ch = BLAST_FASTA(concat_fasta_ch.output, blast_fasta_params)
        // 5. Validate original group hits against concatenated BLAST results
        distance_params = [
            taxid_field_1: "aligner_taxid_lca",
            taxid_field_2: "validation_staxid_lca",
            distance_field_1: "validation_distance_aligner",
            distance_field_2: "validation_distance_validation"
        ]
        validate_params = [
            ref_dir: ref_dir,
            distance_params: distance_params
        ]
        validate_ch = VALIDATE_CLUSTER_REPRESENTATIVES(groups, blast_ch.lca, validate_params)
        // 6. Propagate validation information back to individual hits
        propagate_ch = PROPAGATE_VALIDATION_INFORMATION(groups, concat_cluster_ch.output,
            validate_ch.output, "aligner_taxid_lca")
        // 7. Cleanup and generate final outputs
        regrouped_drop_ch = SELECT_TSV_COLUMNS(propagate_ch.output, "taxid_species,selected_taxid", "drop").output 
        output_hits_ch = COPY_HITS(regrouped_drop_ch, "validation_hits.tsv.gz")
        output_blast_ch = COPY_BLAST(blast_ch.blast, "validation_blast.tsv.gz")
    emit:
        // Main output
        annotated_hits = output_hits_ch
        // Intermediate output
        blast_results = output_blast_ch
        // Extra outputs for testing
        test_in   = groups
        test_split_tsv = split_ch.tsv
        test_cluster_tab = cluster_ch.tsv
        test_reps_fasta = cluster_ch.fasta
        test_concat_fasta = concat_fasta_ch.output
        test_concat_cluster = concat_cluster_ch.output
        test_blast_db = blast_ch.blast
        test_blast_query = blast_ch.query
        test_blast_lca = blast_ch.lca
        test_validate = validate_ch.output
        test_propagate = propagate_ch.output
}
