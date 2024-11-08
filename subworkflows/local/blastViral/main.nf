/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED_MERGED } from "../../../modules/local/subsetReads" addParams(suffix: "fasta")
include { BLAST_PAIRED_LOCAL } from "../../../modules/local/blast" addParams(cpus: params.blast_cpus, mem: params.blast_mem)
include { FILTER_BLAST } from "../../../modules/local/filterBlast" addParams(mem: params.blast_filter_mem)
include { PAIR_BLAST } from "../../../modules/local/pairBlast"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL {
    take:
        viral_fasta
        blast_db_dir
        blast_db_prefix
        read_fraction
    main:
        // Subset viral reads for BLAST
        if ( read_fraction < 1 ) {
            subset_ch = SUBSET_READS_PAIRED_MERGED(viral_fasta, read_fraction)
        } else {
            subset_ch = viral_fasta
        }
        // BLAST putative viral hits against prepared DB
        blast_ch = BLAST_PAIRED_LOCAL(subset_ch, blast_db_dir, blast_db_prefix)
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