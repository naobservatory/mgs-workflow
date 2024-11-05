/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED_MERGED } from "../../../modules/local/subsetReads" addParams(suffix: "fasta")
include { FASTQ_TO_FASTA } from "../../../modules/local/fastqToFasta"
include { BLAST_PAIRED_NT } from "../../../modules/local/blast" addParams(cpus: params.blast_cpus, mem: params.blast_mem)
include { FILTER_BLAST } from "../../../modules/local/filterBlast" addParams(mem: params.blast_filter_mem)
include { PAIR_BLAST } from "../../../modules/local/pairBlast"

/***********
| WORKFLOW |
***********/

workflow BLAST_HV {
    take:
        hv_fasta
        blast_nt_dir
        read_fraction
    main:
        // Subset HV reads for BLAST
        subset_ch = FASTQ_TO_FASTA(hv_fasta)
        // BLAST putative HV hits against nt
        blast_ch = BLAST_PAIRED_NT(subset_ch.fastas, blast_nt_dir)
        // Process BLAST output
        filter_ch = FILTER_BLAST(blast_ch)
        pair_ch = PAIR_BLAST(filter_ch)
    emit:
        blast_raw = blast_ch
        blast_filtered = filter_ch
        blast_paired = pair_ch
    publish:
        subset_ch >> "results/hv"
        pair_ch >> "results/hv"
}
