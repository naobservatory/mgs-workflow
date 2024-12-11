/***********************************************************************************************
| WORKFLOW: PREPROCESSING ON SHORT-READ MGS DATA (EITHER SINGLE-END OR PAIRED-END) |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { RAW } from "../subworkflows/local/raw"
include { CLEAN } from "../subworkflows/local/clean"
include { PROCESS_OUTPUT } from "../subworkflows/local/processOutput"
include { PROFILE } from "../subworkflows/local/profile"
include { LOAD_SAMPLESHET } from "../subworkflows/local/loadSampleSheet"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_DEV_SE {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

    // Check if grouping column exists in samplesheet
    check_grouping = new File(params.sample_sheet).text.readLines()[0].contains('group') ? true : false
    if (params.grouping != check_grouping) {
        if (params.grouping && !check_grouping) {
            throw new Exception("Grouping enabled in config file, but group column absent from samplesheet.")
        } else if (!params.grouping && check_grouping) {
            throw new Exception("Grouping is not enabled in config file, but group column is present in the samplesheet.")
        }
    // Prepare Kraken DB
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    }

    // Load samplesheet
    LOAD_SAMPLESHET(params.sample_sheet)
    samplesheet_ch = LOAD_SAMPLESHET.out.samplesheet
    group_ch = LOAD_SAMPLESHET.out.group

    // Preprocessing
    RAW(samplesheet_ch, params.n_reads_trunc, "2", "4 GB", "raw_concat", params.single_end)
    CLEAN(RAW.out.reads, params.adapters, "2", "4 GB", "cleaned", params.single_end)

    // Taxonomic profiling
    PROFILE(CLEAN.out.reads, group_ch, kraken_db_path, params.n_reads_profile, params.ref_dir, "0.4", "27", "ribo", params.grouping, params.single_end)

    // Process output
    qc_ch = RAW.out.qc.concat(CLEAN.out.qc)
    PROCESS_OUTPUT(qc_ch)

    // Publish results
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
        Channel.fromPath(params.sample_sheet) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // Intermediate files
        CLEAN.out.reads >> "intermediates/reads/cleaned"
        // QC
        PROCESS_OUTPUT.out.basic >> "results"
        PROCESS_OUTPUT.out.adapt >> "results"
        PROCESS_OUTPUT.out.qbase >> "results"
        PROCESS_OUTPUT.out.qseqs >> "results"
        // Final results
        PROFILE.out.bracken >> "results"
        PROFILE.out.kraken >> "results"
}
