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

include { SPLIT_VIRAL_TSV_BY_SPECIES } from "../../../subworkflows/local/splitViralTsvBySpecies"
include { CLUSTER_VIRAL_ASSIGNMENTS } from "../../../subworkflows/local/clusterViralAssignments"
include { BLAST_FASTA } from "../../../subworkflows/local/blastFasta"
include { VALIDATE_CLUSTER_REPRESENTATIVES } from "../../../subworkflows/local/validateClusterRepresentatives"
include { PROPAGATE_VALIDATION_INFORMATION } from "../../../subworkflows/local/propagateValidationInformation"
include { CONCATENATE_TSVS_ACROSS_SPECIES as CONCATENATE_CLUSTERING_INFO } from "../../../subworkflows/local/concatenateTsvsAcrossSpecies"
include { CONCATENATE_TSVS_ACROSS_SPECIES as CONCATENATE_BLAST_RESULTS } from "../../../subworkflows/local/concatenateTsvsAcrossSpecies"

/***********
| WORKFLOW |
***********/

workflow VALIDATE_VIRAL_ASSIGNMENTS {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        db // Viral taxonomy DB
        cluster_identity // Identity threshold for VSEARCH clustering
        cluster_min_len // Minimum sequence length for VSEARCH clustering
        n_clusters // Number of cluster representatives to validate for each specie
        ref_dir // Path to reference directory containing BLAST DB
        blast_db_prefix // Prefix for BLAST reference DB files (e.g. "nt")
        perc_id // Minimum %ID required for BLAST to return an alignment
        qcov_hsp_perc // Minimum query coverage required for BLAST to return an alignment
        blast_max_rank // Only keep alignments that are in the top-N for that query by bitscore
        blast_min_frac // Only keep alignments that have at least this fraction of the best bitscore for that query
        taxid_artificial // Parent taxid for artificial sequences in NCBI taxonomy
    main:
        // 1. Split viral hits TSV by species
        split_ch = SPLIT_VIRAL_TSV_BY_SPECIES(groups, db)
        // 2. Cluster sequences within species and obtain representatives of largest clusters
        cluster_ch = CLUSTER_VIRAL_ASSIGNMENTS(split_ch.fastq, cluster_identity,
            cluster_min_len, n_clusters, Channel.of(false))
        // 3. BLAST cluster representatives
        blast_ch = BLAST_FASTA(cluster_ch.fasta, ref_dir, blast_db_prefix,
            perc_id, qcov_hsp_perc, blast_max_rank, blast_min_frac, taxid_artificial,
            "validation")
        // 4. Validate hit TSV taxids against BLAST results
        validate_ch = VALIDATE_CLUSTER_REPRESENTATIVES(split_ch.tsv, blast_ch.lca, 
            "bowtie2_taxid_best", // Column header for original taxid in hits TSV
            "validation_staxid_lca_natural", // LCA taxid computed from BLAST results, excluding artificial sequences
            "validation_distance", ref_dir)
        // 5. Propagate validation information back to individual hits
        propagate_ch = PROPAGATE_VALIDATION_INFORMATION(split_ch.tsv, cluster_ch.tsv,
            validate_ch.output, "bowtie2_taxid_best")
        // 6. Concatenate validation info and BLAST results across species (to regenerate per-group information)
        regrouped_ch = CONCATENATE_CLUSTERING_INFO(propagate_ch.output, "validation")
        regrouped_blast_ch = CONCATENATE_BLAST_RESULTS(blast_ch.blast, "validation")
    emit:
        // Main output
        annotated_hits = regrouped_ch.output
        // Intermediate output
        blast_results = regrouped_blast_ch.output
        // Extra outputs for testing
        test_in   = groups
        test_db   = db
        test_split_tsv = split_ch.tsv
        test_split_fastq = split_ch.fastq
        test_cluster_tab = cluster_ch.tsv
        test_cluster_ids = cluster_ch.ids
        test_reps_fastq = cluster_ch.fastq
        test_reps_fasta = cluster_ch.fasta
        test_blast_db = blast_ch.blast
        test_blast_query = blast_ch.query
        test_blast_lca = blast_ch.lca
        test_validate = validate_ch.output
        test_propagate = propagate_ch.output
}
