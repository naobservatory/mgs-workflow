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
include { MAKE_HUMAN_BOWTIE2_INDEX } from "../subworkflows/local/makeHumanBowtie2Index"
include { MAKE_CONTAMINANT_BOWTIE2_INDEX } from "../subworkflows/local/makeContaminantBowtie2Index"
include { MAKE_VIRUS_BOWTIE2_INDEX } from "../subworkflows/local/makeVirusBowtie2Index"
include { EXTRACT_TARBALL as EXTRACT_KRAKEN_DB } from "../modules/local/extractTarball"
include { MAKE_RIBO_MINIMAP2_INDEX } from "../subworkflows/local/makeRiboMinimap2Index"
include { MAKE_VIRUS_MINIMAP2_INDEX } from "../subworkflows/local/makeVirusMinimap2Index"
include { MAKE_HUMAN_MINIMAP2_INDEX } from "../subworkflows/local/makeHumanMinimap2Index"
include { MAKE_CONTAMINANT_MINIMAP2_INDEX } from "../subworkflows/local/makeContaminantMinimap2Index"

/****************
| MAIN WORKFLOW |
****************/

workflow INDEX {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Build viral taxonomy and infection DB
    MAKE_VIRUS_TAXONOMY_DB(params.taxonomy_url, params.virus_host_db_url, params.host_taxon_db, params.virus_taxid, params.viral_taxids_exclude)
    // Get reference DB of viral genomes of interest
    MAKE_VIRUS_GENOME_DB(params.ncbi_viral_params, MAKE_VIRUS_TAXONOMY_DB.out.db, params.genome_patterns_exclude, params.host_taxa_screen, params.adapters, "20", "3", "0.5", "10")
    // Build bowtie2 alignment indexes
    MAKE_VIRUS_BOWTIE2_INDEX(MAKE_VIRUS_GENOME_DB.out.fasta)
    MAKE_HUMAN_BOWTIE2_INDEX(params.human_url)
    MAKE_CONTAMINANT_BOWTIE2_INDEX(params.genome_urls, params.contaminants)
    // Other index files
    JOIN_RIBO_REF(params.ssu_url, params.lsu_url)
    DOWNLOAD_BLAST_DB(params.blast_db_name)
    EXTRACT_KRAKEN_DB(params.kraken_db, "kraken_db", true)
    // Build minimap2 indices
    MAKE_VIRUS_MINIMAP2_INDEX(MAKE_VIRUS_GENOME_DB.out.fasta)
    MAKE_HUMAN_MINIMAP2_INDEX(params.human_url)
    MAKE_RIBO_MINIMAP2_INDEX(JOIN_RIBO_REF.out.ribo_ref)
    MAKE_CONTAMINANT_MINIMAP2_INDEX(params.genome_urls, params.contaminants)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "index-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // Taxonomy and virus databases
        MAKE_VIRUS_TAXONOMY_DB.out.db >> "results"
        MAKE_VIRUS_TAXONOMY_DB.out.nodes >> "results"
        MAKE_VIRUS_TAXONOMY_DB.out.names >> "results"
        // Virus genome database
        MAKE_VIRUS_GENOME_DB.out.fasta >> "results"
        MAKE_VIRUS_GENOME_DB.out.metadata >> "results"
        // Bowtie2 alignment indexes
        MAKE_HUMAN_BOWTIE2_INDEX.out.bt2 >> "results"
        MAKE_CONTAMINANT_BOWTIE2_INDEX.out.bt2 >> "results"
        MAKE_VIRUS_BOWTIE2_INDEX.out.bt2 >> "results"
        // Other reference files & directories
        JOIN_RIBO_REF.out.ribo_ref >> "results"
        DOWNLOAD_BLAST_DB.out.db >> "results"
        EXTRACT_KRAKEN_DB.out >> "results"
        // Minimap2 alignment indices
        MAKE_VIRUS_MINIMAP2_INDEX.out.mm2 >> "results"
        MAKE_HUMAN_MINIMAP2_INDEX.out.mm2 >> "results"
        MAKE_RIBO_MINIMAP2_INDEX.out.mm2 >> "results"
        MAKE_CONTAMINANT_MINIMAP2_INDEX.out.mm2 >> "results"
}
