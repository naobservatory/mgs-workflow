/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_FASTN } from "../../../modules/local/subsetFastn"
include { SORT_FASTQ } from "../../../modules/local/sortFastq"
include { CONVERT_FASTQ_FASTA } from "../../../modules/local/convertFastqFasta"
include { BLAST_FASTA } from "../../../subworkflows/local/blastFasta"
include { COPY_FILE } from "../../../modules/local/copyFile"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL {
    take:
        viral_fastq // Interleaved or single-end
        ref_dir // Path to reference directory containing BLAST DB
        blast_db_prefix // Prefix for BLAST reference DB files (e.g. "nt")
        read_fraction
        blast_max_rank
        blast_min_frac
        random_seed
        perc_id
        qcov_hsp_perc
        taxid_artificial // Parent taxid for artificial sequences in NCBI taxonomy
    main:
        // 1. Subset viral reads for BLAST
        reads_in = Channel.of("merged")
            | combine(viral_fastq)
        subset_ch = (read_fraction < 1)
            ? SUBSET_FASTN(reads_in, read_fraction, random_seed).output
            : reads_in
        subset_sorted_ch = SORT_FASTQ(subset_ch).output
        // 2. Convert to FASTA
        fasta_ch = CONVERT_FASTQ_FASTA(subset_sorted_ch).output
        // 3. Run BLAST and process output
        blast_ch = BLAST_FASTA(fasta_ch, ref_dir, blast_db_prefix,
            perc_id, qcov_hsp_perc, blast_max_rank, blast_min_frac,
            taxid_artificial, "validation")
        // 4. Rename subset FASTA file for output
        copy_ch = COPY_FILE(fasta_ch, "blast_input_subset.fasta.gz")
    emit:
        blast_subset = blast_ch.blast
        subset_reads = copy_ch
        test_input = reads_in
}
