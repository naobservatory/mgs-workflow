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
        ref_ch = channel
            .fromList(human_genome_url.entrySet())
            .map { entry -> 
                tuple(entry.value, entry.key)  // (url, name)
            }
        downloaded_ch = DOWNLOAD_GENOME(ref_ch)
        genome_ch = downloaded_ch.collect()
        bbmap_ch = BBMAP_INDEX(genome_ch, "bbm-human-index")
        bowtie2_ch = BOWTIE2_INDEX(genome_ch, "bt2-human-index")
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
