/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTA_GZIPPED } from "../../../modules/local/dustmasker"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2"
include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_INDEX {
    take:
        virus_genome_fasta
    main:
        mask_ch = DUSTMASKER_FASTA_GZIPPED(virus_genome_fasta)
        bowtie2_ch = BOWTIE2_INDEX(mask_ch, "bt2-virus-index")
        minimap2_ch = MINIMAP2_INDEX(mask_ch, "mm2-virus-index")
    emit:
        bt2 = bowtie2_ch
        mm2 = minimap2_ch.output
}
// TODO: Add k-mer index generation
// TODO: Consider changing masking tool
