/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME } from "../../../modules/local/downloadGenome"
include { BBMAP_INDEX } from "../../../modules/local/bbmap"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2"

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_INDEX {
    take:
        human_genome_url
    main:
        genome_ch = DOWNLOAD_GENOME(human_genome_url, "human_genome")
        bbmap_ch = BBMAP_INDEX(genome_ch, "bbm-human-index")
        bowtie2_ch = BOWTIE2_INDEX(genome_ch, "bt2-human-index")
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
