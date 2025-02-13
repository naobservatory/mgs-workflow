/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME } from "../../../modules/local/downloadGenome"
include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_MINIMAP2_INDEX {
    take:
        human_genome_url
    main:
        ref_ch = channel
            .of(tuple(human_genome_url, "human"))  // (url, name)
        genome_ch = DOWNLOAD_GENOME(ref_ch)
        mm2_ch = MINIMAP2_INDEX(genome_ch, "minimap2-human-index")
    emit:
        mm2 = mm2_ch
}
