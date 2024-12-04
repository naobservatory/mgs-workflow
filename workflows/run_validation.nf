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

/****************
| MAIN WORKFLOW |
****************/

// Complete primary workflow
workflow RUN_VALIDATION {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    if ( params.viral_tsv_collapsed == "" ) {
    // Option 1: Directly specify FASTA paths in config file (only used if no RUN output DB given)
        fasta_ch = Channel.value([file(params.viral_fasta_1), file(params.viral_fasta_2)])
    } else {
    // Option 2: Extract read sequences from output DB from RUN workflow (default)
        // Define input
        collapsed_ch = params.viral_tsv_collapsed
        // Extract virus reads into FASTA format
        fasta_ch = MAKE_VIRUS_READS_FASTA(collapsed_ch)
    }
    // BLAST validation on host-viral reads
    if ( params.blast_viral_fraction > 0 ) {
        blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"
        BLAST_VIRAL(fasta_ch, blast_db_path, params.blast_db_prefix, params.blast_viral_fraction)
    }
    // Publish results (NB: BLAST workflow has its own publish directive)
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    index_params_ch = Channel.fromPath("${params.ref_dir}/input/index-params.json")
    .map { file -> file.copyTo("${params.base_dir}/work/params-index.json") }
    index_pipeline_version_ch = Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt")
    .map { file -> file.copyTo("${params.base_dir}/work/pipeline-version-index.txt") }
    publish:
        // Saved inputs
        index_params_ch >> "input"
        index_pipeline_version_ch >> "logging"
        // Saved outputs
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // BLAST outputs
        BLAST_VIRAL.out.blast_subset >> "results"
        BLAST_VIRAL.out.blast_paired >> "results"
}
