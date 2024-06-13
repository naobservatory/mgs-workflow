/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED } from "../modules/local/subsetReads"
include { BLAST_PAIRED_NT } from "../modules/local/blast" addParams(cpus: params.blast_cpus, mem: params.blast_mem)
include { FILTER_BLAST } from "../modules/local/filterBlast" addParams(mem: params.blast_filter_mem)
include { PAIR_BLAST } from "../modules/local/pairBlast"

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
        subset_ch = SUBSET_READS_PAIRED(hv_fasta, read_fraction)
        // BLAST putative HV hits against nt
        blast_ch = BLAST_PAIRED_NT(subset_ch, blast_nt_dir)
        // Process BLAST output
        filter_ch = FILTER_BLAST(blast_ch)
        pair_ch = PAIR_BLAST(filter_ch)
    emit:
        blast_raw = blast_ch
        blast_filtered = filter_ch
        blast_paired = pair_ch
}
