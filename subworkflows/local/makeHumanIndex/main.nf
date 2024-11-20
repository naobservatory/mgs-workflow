/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_SINGLE_GENOME } from "../../../modules/local/downloadSingleGenome"
include { BBMAP_INDEX } from "../../../modules/local/bbmap"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2"

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_INDEX {
    take:
        human_genome_url
        bbmap_mem
    main:
        genome_ch = DOWNLOAD_SINGLE_GENOME(human_genome_url, "human_genome")
        bbmap_ch = BBMAP_INDEX(genome_ch, "bbm-human-index", bbmap_mem)
        bowtie2_ch = BOWTIE2_INDEX(genome_ch, "bt2-human-index")
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
