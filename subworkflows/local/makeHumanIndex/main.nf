/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME } from "../../../modules/local/downloadGenome" addParams(name: "human_genome")
include { BBMAP_INDEX } from "../../../modules/local/bbmap" addParams(outdir: "bbm-human-index")
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2" addParams(outdir: "bt2-human-index")

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_INDEX {
    take:
        human_genome_url
    main:
        genome_ch = DOWNLOAD_GENOME(human_genome_url)
        bbmap_ch = BBMAP_INDEX(genome_ch)
        bowtie2_ch = BOWTIE2_INDEX(genome_ch)
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
