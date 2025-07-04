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
include { CONCATENATE_FILES_ACROSS_SPECIES } from "../../../subworkflows/local/concatenateFilesAcrossSpecies"
include { CONCATENATE_TSVS_ACROSS_SPECIES } from "../../../subworkflows/local/concatenateTsvsAcrossSpecies"
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
        // 3. Concatenate data across species (prepare for group-level BLAST)
        concat_fasta_ch = CONCATENATE_FILES_ACROSS_SPECIES(cluster_ch.fasta, "cluster_reps")
        concat_cluster_ch = CONCATENATE_TSVS_ACROSS_SPECIES(cluster_ch.tsv, "cluster_info")
        // 4. Run BLAST on concatenated cluster representatives (single job per group)
        blast_ch = BLAST_FASTA(concat_fasta_ch.output, ref_dir, blast_db_prefix,
            perc_id, qcov_hsp_perc, blast_max_rank, blast_min_frac, taxid_artificial,
            "validation")
        // 5. Validate original group hits against concatenated BLAST results
        distance_params = [
            taxid_field_1: "aligner_taxid",
            taxid_field_2: "validation_staxid_lca_natural",
            distance_field_1: "validation_distance_aligner",
            distance_field_2: "validation_distance_validation"
        ]
        validate_ch = VALIDATE_CLUSTER_REPRESENTATIVES(groups, blast_ch.lca,
            ref_dir, distance_params)
        // 6. Propagate validation information back to individual hits
        propagate_ch = PROPAGATE_VALIDATION_INFORMATION(groups, concat_cluster_ch.output,
            validate_ch.output, "aligner_taxid")
        // 7. Cleanup and generate final outputs
        regrouped_drop_ch = SELECT_TSV_COLUMNS(propagate_ch.output, "taxid_species", "drop").output
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
