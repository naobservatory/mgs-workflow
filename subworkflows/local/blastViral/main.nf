/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_FASTQ } from "../../../modules/local/subsetFastq"
include { CONVERT_FASTQ_FASTA } from "../../../modules/local/convertFastqFasta"
include { BLASTN as BLAST } from "../../../modules/local/blast"
include { SORT_FILE as SORT_BLAST_1 } from "../../../modules/local/sortFile"
include { SORT_FILE as SORT_BLAST_2 } from "../../../modules/local/sortFile"
include { FILTER_TSV } from "../../../modules/local/filterTsv"
include { FILTER_BLAST } from "../../../modules/local/filterBlast"
include { SORT_FASTQ } from "../../../modules/local/sortFastq"
include { COPY_FILE } from "../../../modules/local/copyFile"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL {
    take:
        viral_fastq // Interleaved or single-end
        blast_db_dir
        blast_db_prefix
        read_fraction
        blast_max_rank
        blast_min_frac
        random_seed
        perc_id
        qcov_hsp_perc
    main:
        // 1. Subset viral reads for BLAST
        reads_in = Channel.of("merged")
            | combine(viral_fastq)
        subset_ch = (read_fraction < 1)
            ? SUBSET_FASTQ(reads_in, read_fraction, random_seed).output
            : reads_in
        subset_sorted_ch = SORT_FASTQ(subset_ch).output
        // 2. Convert to FASTA
        fasta_ch = CONVERT_FASTQ_FASTA(subset_sorted_ch).output
        // 3. Run BLAST
        blast_ch = BLAST(fasta_ch, blast_db_dir, blast_db_prefix, perc_id, qcov_hsp_perc)
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
        filter_ch_2 = FILTER_BLAST(sort_ch_2.output, blast_max_rank, blast_min_frac).output
        // Rename subset FASTA file for output
        copy_ch = COPY_FILE(fasta_ch, "blast_input_subset.fasta.gz")
    emit:
        blast_subset = filter_ch_2
        subset_reads = copy_ch
        test_input = reads_in
}
