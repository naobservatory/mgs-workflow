/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_RIBO_INDEX {
    take:
        ribo_ref
    main:
        minimap2_ch = MINIMAP2_INDEX(ribo_ref_uniq_headers_ch,"mm2-ribo-index")
    emit:
        mm2 = minimap2_ch
}
