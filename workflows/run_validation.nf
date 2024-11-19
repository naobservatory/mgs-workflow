/*****************************************************
| WORKFLOW: POST-HOC VALIDATION OF PUTATIVE viral READS |
*****************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MAKE_VIRUS_READS_FASTA } from "../modules/local/makeVirusReadsFasta"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_VALIDATION {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Define input
    collapsed_ch = params.viral_tsv_collapsed
    // Extract virus reads into FASTA format
    fasta_ch = MAKE_VIRUS_READS_FASTA(collapsed_ch)
    // BLAST validation on host-viral reads
    if ( params.blast_viral_fraction > 0 ) {
        blast_db_path = "${params.ref_dir}/results/core_nt"
        blast_db_prefix = "core_nt"
        BLAST_VIRAL(fasta_ch, blast_nt_path, params.blast_viral_fraction)
        BLAST_VIRAL(fasta_ch, blast_db_path, blast_db_prefix, params.blast_viral_fraction, "32", "256 GB", "32 GB")
    }
    // Publish results (NB: BLAST workflow has its own publish directive)
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        Channel.fromPath("${params.ref_dir}/input/index-params.json") >> "input"
        Channel.fromPath("${params.ref_dir}/input/pipeline-version.txt").collectFile(name: "pipeline-version-index.txt") >> "input"
        params_ch >> "input"
        time_ch >> "input"
        version_ch >> "input"
}
