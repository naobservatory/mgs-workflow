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
include { EXTRACT_VIRAL_READS } from "../subworkflows/local/extractViralReads"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
include { PROFILE } from "../subworkflows/local/profile"
include { PROCESS_OUTPUT } from "../subworkflows/local/processOutput"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_DEV_SE {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

    // Determine read type based on samplesheet header
    single_end = file(params.sample_sheet).readLines()[0].contains('fastq_2') ? false : true

    println "Single end mode: ${single_end}"

    // Prepare samplesheet
    if (single_end) {
        if (params.grouping) {
            samplesheet = Channel
                .fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq), row.group) }
            samplesheet_ch = samplesheet.map { sample, read, group -> tuple(sample, [read]) }
            group_ch = samplesheet.map { sample, read, group -> tuple(sample, group) }
        } else {
            samplesheet = Channel
                .fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq)) }
            samplesheet_ch = samplesheet.map { sample, read -> tuple(sample, [read]) }
            group_ch = Channel.empty()
        }
    } else {
        if (params.grouping) {
            samplesheet = Channel
                .fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2), row.group) }
            samplesheet_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, [read1, read2]) }
            group_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, group) }
        } else {
            samplesheet = Channel
                .fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
            samplesheet_ch = samplesheet.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
            group_ch = Channel.empty()
            }
        }


    // Preprocessing
    RAW(samplesheet_ch, params.n_reads_trunc, "2", "4 GB", "raw_concat", single_end)
    CLEAN(RAW.out.reads, params.adapters, "2", "4 GB", "cleaned", single_end)

    // Process output
    qc_ch = RAW.out.qc.concat(CLEAN.out.qc)
    PROCESS_OUTPUT(qc_ch)

    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")

    publish:
        // Saved inputs
        Channel.fromPath("${params.ref_dir}/input/index-params.json") >> "input"
        Channel.fromPath("${params.ref_dir}/input/pipeline-version.txt").collectFile(name: "pipeline-version-index.txt") >> "logging"
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
}
