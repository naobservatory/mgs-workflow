/*
Given a FASTA file of query sequences and a BLAST DB,
align the former against the latter, then sort and filter
the tabular BLAST output based on bitscore.
*/


/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BLASTN as BLAST } from "../../../modules/local/blast"
include { SORT_FILE as SORT_BLAST_1 } from "../../../modules/local/sortFile"
include { SORT_FILE as SORT_BLAST_2 } from "../../../modules/local/sortFile"
include { FILTER_TSV } from "../../../modules/local/filterTsv"
include { FILTER_BLAST } from "../../../modules/local/filterBlast"

/***********
| WORKFLOW |
***********/

workflow BLAST_FASTA {
    take:
        query_fasta // Fasta of sequences to align
        blast_db_dir // Path to BLAST reference DB
        blast_db_prefix // Prefix for BLAST reference DB files (e.g. "nt")
        perc_id // Minimum %ID required for BLAST to return an alignment
        qcov_hsp_perc // Minimum query coverage required for BLAST to return an alignment
        blast_max_rank // Only keep alignments that are in the top-N for that query by bitscore
        blast_min_frac // Only keep alignments that have at least this fraction of the best bitscore for that query
    main:
        // 1. Run BLAST
        blast_ch = BLAST(query_fasta, blast_db_dir, blast_db_prefix,
            perc_id, qcov_hsp_perc)
        // 2. Filter BLAST output to only one row per query/subject combination
        sort_str_1 = "-t\$\'\\t\' -k1,1 -k2,2 -k7,7nr -k9,9nr" // Sort by query, subject, bitscore, length
        sort_ch_1 = SORT_BLAST_1(blast_ch.output, sort_str_1, "blast")
        filter_ch_1 = FILTER_TSV(sort_ch_1.output, "1,2", "blast") // Take first (best) row for each query/subject combo
        // 3. Filter BLAST output to only high-scoring alignments within each query
        sort_str_2 = "-t\$\'\\t\' -k1,1 -k7,7nr" // Sort by query and bitscore only
        sort_ch_2 = SORT_BLAST_2(filter_ch_1.output, sort_str_2, "blast")
        filter_ch_2 = FILTER_BLAST(sort_ch_2.output, blast_max_rank, blast_min_frac).output // Filter on relative bitscores
    emit:
        blast = filter_ch_2
        query = query_fasta
}
