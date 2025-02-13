/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTA_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_MINIMAP2_INDEX {
    take:
        virus_genome_fasta
    main:
        // Need to check if masking works with minimap2 indexing.
        mask_ch = DUSTMASKER_FASTA_GZIPPED(virus_genome_fasta)
        mm2_ch = MINIMAP2_INDEX(mask_ch, "minimap2-virus-index")
    emit:
        mm2 = mm2_ch
}
