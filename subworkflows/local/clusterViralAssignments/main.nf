/*
Cluster a channel of viral sequences with VSEARCH, process the output,
and extract the representative sequences of the top N largest clusters.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MERGE_JOIN_READS } from "../../../subworkflows/local/mergeJoinReads"
include { VSEARCH_CLUSTER } from "../../../modules/local/vsearch"
include { PROCESS_VSEARCH_CLUSTER_OUTPUT } from "../../../modules/local/processVsearchClusterOutput"
include { DOWNSAMPLE_FASTN_BY_ID } from "../../../modules/local/downsampleFastnById"
include { CONVERT_FASTQ_FASTA } from "../../../modules/local/convertFastqFasta"

/***********
| WORKFLOW |
***********/

workflow CLUSTER_VIRAL_ASSIGNMENTS {
    take:
        reads_ch // Single-end or interleaved FASTQ sequences
        params_map // Map containing clustering parameters
        single_end // Is the input read data single-ended (true) or interleaved (false)?
    main:
        // Extract parameters from map
        cluster_identity = params_map.cluster_identity
        cluster_min_len = params_map.cluster_min_len
        n_clusters = params_map.n_clusters
        
        // 1. Merge and join interleaved sequences to produce a single sequence per input pair
        merge_ch = MERGE_JOIN_READS(reads_ch, single_end)
        // 2. Cluster merged reads
        vsearch_params = [
            identity_threshold: cluster_identity,
            identity_method: 0,
            min_seq_length: cluster_min_len
        ]
        cluster_ch = VSEARCH_CLUSTER(merge_ch.single_reads, vsearch_params)
        // 3. Extract clustering information and representative IDs
        cluster_info_ch = PROCESS_VSEARCH_CLUSTER_OUTPUT(cluster_ch.summary, n_clusters, "vsearch")
        // 4. Extract representative sequences for the N largest clusters for each species
        id_prep_ch = merge_ch.single_reads.combine(cluster_info_ch.ids, by: 0)
        rep_fastq_ch = DOWNSAMPLE_FASTN_BY_ID(id_prep_ch).output
        rep_fasta_ch = CONVERT_FASTQ_FASTA(rep_fastq_ch).output
    emit:
        tsv = cluster_info_ch.output
        ids = cluster_info_ch.ids
        fastq = rep_fastq_ch
        fasta = rep_fasta_ch
        test_merged = merge_ch.single_reads
        test_cluster_reps = cluster_ch.reps
        test_cluster_summ = cluster_ch.summary
}
