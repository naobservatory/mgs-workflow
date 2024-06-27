/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_NCBI_TAXONOMY } from "../../../modules/local/downloadNcbiTaxonomy"
include { EXTRACT_NCBI_TAXONOMY } from "../../../modules/local/extractNcbiTaxonomy"

/***********
| WORKFLOW |
***********/

workflow GET_NCBI_TAXONOMY {
    take:
        taxonomy_url
    main:
        dl_ch = DOWNLOAD_NCBI_TAXONOMY(taxonomy_url)
        ext_ch = EXTRACT_NCBI_TAXONOMY(dl_ch)
    emit:
        nodes = ext_ch.nodes
        names = ext_ch.names
}
