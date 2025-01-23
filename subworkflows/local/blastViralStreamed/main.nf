/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_FASTQ } from "../../../modules/local/subsetFastq"
include { CONVERT_FASTQ_FASTA } from "../../../modules/local/convertFastqFasta"
include { BLASTN_STREAMED as BLAST } from "../../../modules/local/blast"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL {
    take:
        viral_fastq // Interleaved or single-end
        blast_db_dir
        blast_db_prefix
        read_fraction
    main:
        // 1. Subset viral reads for BLAST
        reads_in = Channel.of("merged")
            | combine(viral_fastq)
        subset_ch = (read_fraction < 1)
            ? SUBSET_FASTQ(reads_in, read_fraction)
            : reads_in
        // 2. Convert to FASTA
        fasta_ch = CONVERT_FASTQ_FASTA(subset_ch)
        // 3. Run BLAST
        blast_ch = BLAST(fasta_ch, blast_db_dir, blast_db_prefix)
        // Filter BLAST output (keeping forward and reverse separate)
        filter_ch_1 = FILTER_BLAST_FWD(blast_ch_1, 5)
        filter_ch_2 = FILTER_BLAST_REV(blast_ch_2, 5)
        // Combine alignment information across pairs
        filter_out_1 = filter_ch_1.map{ file ->
            def newFile = file.parent / "fwd_${file.name}"
            file.copyTo(newFile)
            return newFile
            }
        filter_out_2 = filter_ch_2.map{ file ->
            def newFile = file.parent / "rev_${file.name}"
            file.copyTo(newFile)
            return newFile
            }
        pair_ch = PAIR_BLAST(filter_out_1, filter_out_2)
    emit:
        blast_subset = subset_ch
        blast_paired = pair_ch
}
