/**************************************************************************************
| WORKFLOW: GENERATE INDEX AND REFERENCE FILES FOR DOWNSTREAM PROCESSING AND ANALYSIS |
**************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

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
include { EXTRACT_TARBALL as EXTRACT_KRAKEN_DB } from "../modules/local/extractTarball"

/****************
| MAIN WORKFLOW |
****************/

workflow INDEX {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
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
    EXTRACT_KRAKEN_DB(params.kraken_db, "kraken_db", true)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "index-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        params_ch >> "input"
        time_ch >> "input"
        version_ch >> "input"
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
        MAKE_HUMAN_VIRUS_INDEX.out.mapping >> "results"
        // Other reference files & directories
        JOIN_RIBO_REF.out.ribo_ref >> "results"
        DOWNLOAD_BLAST_DB.out.db >> "results"
        EXTRACT_KRAKEN_DB.out >> "results"
}
