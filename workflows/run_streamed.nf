/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_SAMPLESHEET } from "../subworkflows/local/loadSampleSheet"
include { COUNT_TOTAL_READS } from "../subworkflows/local/countTotalReads"
include { EXTRACT_VIRAL_READS_STREAMED as EXTRACT_VIRAL_READS } from "../subworkflows/local/extractViralReadsStreamed"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_STREAMED {
    // Setting reference paths
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"

    // Load samplesheet
    LOAD_SAMPLESHEET(params.sample_sheet, params.grouping, params.single_end)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    group_ch = LOAD_SAMPLESHEET.out.group
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch)

    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(CLEAN.out.reads, params.ref_dir, kraken_db_path, params.bt2_score_threshold,
        params.adapters, params.host_taxon, "1", "24", "viral")

    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
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
        // Final results
        EXTRACT_VIRAL_READS_STREAMED.out.test_out >> "results"
}
