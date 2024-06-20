import groovy.json.JsonOutput

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { JOIN_RIBO_REF } from "../modules/local/joinRiboRef"
include { DOWNLOAD_BLAST_DB } from "../modules/local/downloadBlastDB" addParams(db: "nt")
include { MAKE_HUMAN_INDEX } from "../subworkflows/local/makeHumanIndex"
include { MAKE_CONTAMINANT_INDEX } from "../subworkflows/local/makeContaminantIndex"
include { MAKE_HUMAN_VIRUS_DB } from "../subworkflows/local/makeHumanVirusDB"
include { MAKE_TOTAL_VIRUS_DB } from "../modules/local/makeTotalVirusDB"
include { GET_NCBI_TAXONOMY } from "../subworkflows/local/getNcbiTaxonomy"
include { MAKE_HUMAN_VIRUS_INDEX } from "../subworkflows/local/makeHumanVirusIndex"
include { COPY_FILE as COPY_KRAKEN } from "../modules/local/copyFile" addParams(outpath: "kraken-db.tar.gz")

/****************
| MAIN WORKFLOW |
****************/

workflow INDEX {
    // Make human-viral and total-viral reference DBs
    GET_NCBI_TAXONOMY(params.taxonomy_url)
    MAKE_HUMAN_VIRUS_DB(params.virus_host_db_url, GET_NCBI_TAXONOMY.out.nodes, GET_NCBI_TAXONOMY.out.names)
    MAKE_TOTAL_VIRUS_DB(MAKE_HUMAN_VIRUS_DB.out.hv, GET_NCBI_TAXONOMY.out.nodes, GET_NCBI_TAXONOMY.out.names)
    // Alignment indexes
    MAKE_HUMAN_INDEX(params.human_url)
    MAKE_CONTAMINANT_INDEX(params.cow_url, params.pig_url, params.mouse_url, params.carp_url, params.ecoli_url, params.contaminants)
    MAKE_HUMAN_VIRUS_INDEX(MAKE_HUMAN_VIRUS_DB.out.taxids, params.hv_patterns_exclude)
    // Other index files
    JOIN_RIBO_REF(params.ssu_url, params.lsu_url)
    DOWNLOAD_BLAST_DB()
    COPY_KRAKEN(params.kraken_db)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "index_params.json")
    publish:
        // Saved inputs
        params_ch >> "input"
        // Taxonomy and virus databases
        GET_NCBI_TAXONOMY.out.nodes >> "results"
        GET_NCBI_TAXONOMY.out.names >> "results"
        MAKE_HUMAN_VIRUS_DB.out.hv >> "results"
        MAKE_TOTAL_VIRUS_DB.out.db >> "results"
        // Alignment indexes
        MAKE_HUMAN_INDEX.out.bbm >> "results"
        MAKE_HUMAN_INDEX.out.bt2 >> "results"
        MAKE_CONTAMINANT_INDEX.out.bbm >> "results"
        MAKE_CONTAMINANT_INDEX.out.bt2 >> "results"
        MAKE_HUMAN_VIRUS_INDEX.out.bt2 >> "results"
        MAKE_HUMAN_VIRUS_INDEX.out.filtered >> "results"
        // Other reference files & directories
        JOIN_RIBO_REF.out.ribo_ref >> "results"
        DOWNLOAD_BLAST_DB.out.db >> "results"
        COPY_KRAKEN.out.file >> "results"
}
