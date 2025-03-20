/********************************************************
| WORKFLOW: POST-HOC VALIDATION OF PUTATIVE VIRAL READS |
********************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF as EXTRACT_HITS } from "../modules/local/extractViralHitsToFastqNoref"
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

    // Get input FASTQ
    if ( params.viral_tsv == "" ) {
    // Option 1: Directly specify FASTQ path in config (interleaved/single-end)
        fastq_ch = params.viral_fastq
    } else {
    // Option 2: Extract read sequences from output DB from RUN workflow (default)
        // Define input
        tsv_ch = params.viral_tsv
        fastq_ch = EXTRACT_HITS(tsv_ch, params.drop_unpaired).output
    }

    // BLAST validation on host-viral reads
    blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"
    BLAST_VIRAL(fastq_ch, blast_db_path, params.blast_db_prefix,
        params.blast_viral_fraction, params.blast_max_rank, params.blast_min_frac,
        params.blast_perc_id, params.blast_qcov_hsp_perc, params.random_seed)

    // Publish results (NB: BLAST workflow has its own publish directive)
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    index_params_ch = Channel.fromPath("${params.ref_dir}/input/index-params.json")
        | map { file -> file.copyTo("${params.base_dir}/work/params-index.json") }
    index_pipeline_version_ch = Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt")
        | map { file -> file.copyTo("${params.base_dir}/work/pipeline-version-index.txt") }
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
        BLAST_VIRAL.out.subset_reads >> "results"
}
