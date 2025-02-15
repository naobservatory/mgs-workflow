/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"
include { MARK_DUPLICATE_FASTA_HEADERS } from "../../../modules/local/markDuplicateHeaders"

/***********
| WORKFLOW |
***********/

workflow MAKE_RIBO_INDEX {
    take:
        ribo_ref
    main:
        ribo_ref_uniq_headers_ch = MARK_DUPLICATE_FASTA_HEADERS(ribo_ref)
        minimap2_ch = MINIMAP2_INDEX(ribo_ref_uniq_headers_ch,"mm2-ribo-index")
    emit:
        mm2 = minimap2_ch
}
