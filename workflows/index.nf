/**************************************************************************************
| WORKFLOW: GENERATE INDEX AND REFERENCE FILES FOR DOWNSTREAM PROCESSING AND ANALYSIS |
**************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MAKE_VIRUS_TAXONOMY_DB } from "../subworkflows/local/makeVirusTaxonomyDB"
include { MAKE_VIRUS_GENOME_DB } from "../subworkflows/local/makeVirusGenomeDB"
include { JOIN_RIBO_REF } from "../modules/local/joinRiboRef"
include { DOWNLOAD_BLAST_DB } from "../modules/local/downloadBlastDB"
include { MAKE_HUMAN_INDEX } from "../subworkflows/local/makeHumanIndex"
include { MAKE_CONTAMINANT_INDEX } from "../subworkflows/local/makeContaminantIndex"
include { MAKE_VIRUS_INDEX } from "../subworkflows/local/makeVirusIndex"
include { MAKE_RIBO_INDEX } from "../subworkflows/local/makeRiboIndex"
include { GET_TARBALL as GET_KRAKEN_DB } from "../modules/local/getTarball"

/****************
| MAIN WORKFLOW |
****************/

workflow INDEX {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Build viral taxonomy and infection DB
    MAKE_VIRUS_TAXONOMY_DB(params.taxonomy_url, params.virus_host_db_url,
        params.host_taxon_db, params.virus_taxid,
        params.viral_taxids_exclude_hard)
    // Get reference DB of viral genomes of interest
    MAKE_VIRUS_GENOME_DB(params.ncbi_viral_params, MAKE_VIRUS_TAXONOMY_DB.out.db, params.genome_patterns_exclude, params.host_taxa_screen, params.adapters, "20", "3", "0.5", "10")
    // Build alignment indices
    JOIN_RIBO_REF(params.ssu_url, params.lsu_url)
    MAKE_VIRUS_INDEX(MAKE_VIRUS_GENOME_DB.out.fasta)
    MAKE_HUMAN_INDEX(params.human_url)
    MAKE_CONTAMINANT_INDEX(params.genome_urls, params.contaminants)
    MAKE_RIBO_INDEX(JOIN_RIBO_REF.out.ribo_ref)
    // Other index files
    DOWNLOAD_BLAST_DB(params.blast_db_name)
    GET_KRAKEN_DB(params.kraken_db, "kraken_db", true)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "index-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    compatibility_ch = Channel.fromPath("${projectDir}/index-min-pipeline-version.txt")
    publish:
        // Saved inputs
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        compatibility_ch >> "logging"
        // Taxonomy and virus databases
        MAKE_VIRUS_TAXONOMY_DB.out.db >> "results"
        MAKE_VIRUS_TAXONOMY_DB.out.nodes >> "results"
        MAKE_VIRUS_TAXONOMY_DB.out.names >> "results"
        // Virus genome database
        MAKE_VIRUS_GENOME_DB.out.fasta >> "results"
        MAKE_VIRUS_GENOME_DB.out.metadata >> "results"
        // Other reference files & directories
        JOIN_RIBO_REF.out.ribo_ref >> "results"
        DOWNLOAD_BLAST_DB.out.db >> "results"
        GET_KRAKEN_DB.out >> "results"
        // Bowtie2 alignment indexes
        MAKE_HUMAN_INDEX.out.bt2 >> "results"
        MAKE_CONTAMINANT_INDEX.out.bt2 >> "results"
        MAKE_VIRUS_INDEX.out.bt2 >> "results"
        // Minimap2 alignment indices
        MAKE_VIRUS_INDEX.out.mm2 >> "results"
        MAKE_HUMAN_INDEX.out.mm2 >> "results"
        MAKE_RIBO_INDEX.out.mm2 >> "results"
        MAKE_CONTAMINANT_INDEX.out.mm2 >> "results"
}
