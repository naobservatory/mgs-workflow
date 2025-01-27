/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_FASTQ } from "../../../modules/local/subsetFastq"
include { CONVERT_FASTQ_FASTA } from "../../../modules/local/convertFastqFasta"
include { BLASTN_STREAMED as BLAST } from "../../../modules/local/blast"
include { SORT_FILE as SORT_BLAST_1 } from "../../../modules/local/sortFile"
include { SORT_FILE as SORT_BLAST_2 } from "../../../modules/local/sortFile"
include { FILTER_TSV } from "../../../modules/local/filterTsv"
include { FILTER_BLAST_STREAMED as FILTER_BLAST } from "../../../modules/local/filterBlast"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL_STREAMED {
    take:
        viral_fastq // Interleaved or single-end
        blast_db_dir
        blast_db_prefix
        read_fraction
        blast_max_rank
        blast_min_frac
    main:
        // 1. Subset viral reads for BLAST
        reads_in = Channel.of("merged")
            | combine(viral_fastq)
        subset_ch = (read_fraction < 1)
            ? SUBSET_FASTQ(reads_in, read_fraction)
            : reads_in
        // 2. Convert to FASTA
        fasta_ch = CONVERT_FASTQ_FASTA(subset_ch.output)
        // 3. Run BLAST
        blast_ch = BLAST(fasta_ch.output, blast_db_dir, blast_db_prefix)
        // 4. Sort and filter BLAST output
        // First sort by query, subject, bitscore and length
        sort_str_1 = "-t\$\'\\t\' -k1,1 -k2,2 -k7,7nr -k9,9nr"
        sort_ch_1 = SORT_BLAST_1(blast_ch.output, sort_str_1, "blast")
        // Then take the first (best) row for each query/subject combination
        filter_ch_1 = FILTER_TSV(sort_ch_1.output, "1,2", "blast")
        // Next sort by query and bitscore only
        sort_str_2 = "-t\$\'\\t\' -k1,1 -k7,7nr"
        sort_ch_2 = SORT_BLAST_2(filter_ch_1.output, sort_str_2, "blast")
        // Then filter by bitscore within each query
        filter_ch_2 = FILTER_BLAST(sort_ch_2.output, blast_max_rank, blast_min_frac)
        // TODO: Incomplete from here
        // pair_ch = PAIR_BLAST(filter_out_1, filter_out_2)
    emit:
        blast_subset = filter_ch_2.output
        // blast_paired = pair_ch
}
