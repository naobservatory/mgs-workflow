/*
Perform efficient post-hoc validation of putative viral reads identified by the RUN workflow.

A. Partition putative hits by assigned species
B. Extract into FASTQ and merge read pairs into single sequences
C. Cluster sequences from each species and identify representative sequences
D. [TODO] Align representative sequences against a large reference DB
E. [TODO] Compare taxids assigned in (4) to those assigned by RUN workflow
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SPLIT_VIRAL_TSV_BY_SPECIES } from "../../../subworkflows/local/splitViralTsvBySpecies"
include { CLUSTER_VIRAL_ASSIGNMENTS } from "../../../subworkflows/local/clusterViralAssignments"
include { BLAST_FASTA } from "../../../subworkflows/local/blastFasta"
include { ADD_SAMPLE_COLUMN as LABEL_GROUP_SPECIES } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS_LABELED } from "../../../modules/local/concatenateTsvs"
include { ADD_SAMPLE_COLUMN as LABEL_GROUP } from "../../../modules/local/addSampleColumn"

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
        blast_db_dir // Path to BLAST reference DB
        blast_db_prefix // Prefix for BLAST reference DB files (e.g. "nt")
        perc_id // Minimum %ID required for BLAST to return an alignment
        qcov_hsp_perc // Minimum query coverage required for BLAST to return an alignment
        blast_max_rank // Only keep alignments that are in the top-N for that query by bitscore
        blast_min_frac // Only keep alignments that have at least this fraction of the best bitscore for that query
    main:
        // 1. Split viral hits TSV by species
        split_ch = SPLIT_VIRAL_TSV_BY_SPECIES(groups, db)
        // 2. Cluster sequences within species and obtain representatives of largest clusters
        cluster_ch = CLUSTER_VIRAL_ASSIGNMENTS(split_ch.fastq, cluster_identity,
            cluster_min_len, n_clusters, Channel.of(false))
        // 3. BLAST cluster representatives
        blast_ch = BLAST_FASTA(cluster_ch.fasta, blast_db_dir, blast_db_prefix,
            perc_id, qcov_hsp_perc, blast_max_rank, blast_min_frac)
        // 4. Concatenate clustering information across species to regenerate per-group information
        // NB: This concatenation stage will move down as more steps are added, but will need to happen eventually
        // and is useful for testing, so I'm implementing it now. It should probably get moved into its own
        // workflow eventually.
        // First label each element with species tag
        to_concat_ch = cluster_ch.tsv
        to_concat_labeled_ch = LABEL_GROUP_SPECIES(to_concat_ch, "group_species", "group_species").output
        // Then change each element from [group_species, path] to [group, path]
        split_label_ch = to_concat_labeled_ch.map{
            label, path ->
                def pattern = /^(.*?)_(\d+)$/
                def matcher = (label =~ pattern)
                if (!matcher) {
                    def msg = "Group label doesn't match required pattern: ${label}, ${path}, ${pattern}"
                    throw new IllegalArgumentException(msg)
                }
                [matcher[0][1], path]
        }
        // Then group elements by group label and concatenate
        regrouped_label_ch = split_label_ch.groupTuple()
        regrouped_concat_ch = CONCATENATE_TSVS_LABELED(regrouped_label_ch, "clusters").output
        regrouped_ch = LABEL_GROUP(regrouped_concat_ch, "group", "clusters").output
        // TODO: Implement D,E from docstring
    emit:
        test_in   = groups
        test_db   = db
        test_split_tsv = split_ch.tsv
        test_split_fastq = split_ch.fastq
        test_cluster_tab = cluster_ch.tsv
        test_cluster_ids = cluster_ch.ids
        test_reps_fastq = cluster_ch.fastq
        test_reps_fasta = cluster_ch.fasta
        test_regrouped = regrouped_ch
}
