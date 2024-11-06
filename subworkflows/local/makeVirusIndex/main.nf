/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTA_GZIPPED } from "../../../modules/local/dustmasker"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2" addParams(outdir: "bt2-virus-index")

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_INDEX {
    take:
        virus_genome_fasta
    main:
        mask_ch = DUSTMASKER_FASTA_GZIPPED(virus_genome_fasta)
        bowtie2_ch = BOWTIE2_INDEX(mask_ch)
    emit:
        bt2 = bowtie2_ch
}
// TODO: Add k-mer index generation
// TODO: Consider changing masking tool
