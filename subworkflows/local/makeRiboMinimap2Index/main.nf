/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_RIBO_MINIMAP2_INDEX {
    take:
        ribo_ref
    main:
        mm2_ch = MINIMAP2_INDEX(ribo_ref, "minimap2-ribo-index")
    emit:
        mm2 = mm2_ch
}
