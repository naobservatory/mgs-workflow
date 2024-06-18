/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_HUMAN_VIRUS_DB } from "../modules/local/downloadHumanVirusDB"
include { EXPAND_HUMAN_VIRUS_DB } from "../modules/local/expandHumanVirusDB"
include { GET_HV_DESCENDENTS } from "../modules/local/getHvDescendents"
include { FINALIZE_HUMAN_VIRUS_DB } from "../modules/local/finalizeHumanVirusDB"

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_VIRUS_DB {
    take:
        virus_host_db_url
        nodes_ch
        names_ch
    main:
        dl_ch = DOWNLOAD_HUMAN_VIRUS_DB(virus_host_db_url)
        exp_ch = EXPAND_HUMAN_VIRUS_DB(dl_ch, nodes_ch, names_ch)
        desc_ch = GET_HV_DESCENDENTS(exp_ch)
        out_ch = FINALIZE_HUMAN_VIRUS_DB(exp_ch, desc_ch.descendents, nodes_ch, names_ch)
    emit:
        hv = out_ch
        taxids = desc_ch.taxids
}
