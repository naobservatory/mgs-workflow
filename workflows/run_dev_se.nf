/***********************************************************************************************
| WORKFLOW: PREPROCESSING ON SHORT-READ MGS DATA (EITHER SINGLE-END OR PAIRED-END) |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { CLEAN } from "../subworkflows/local/clean" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "cleaned")
include { HV } from "../subworkflows/local/hv" addParams(min_kmer_hits: "3", k: "21", bbduk_suffix: "hv", encoding: "${params.quality_encoding}", fuzzy_match: "${params.fuzzy_match_alignment_duplicates}")
include { BLAST_HV } from "../subworkflows/local/blastHV" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")
include { PROFILE } from "../subworkflows/local/profile" addParams(min_kmer_fraction: "0.4", k: "27", bbduk_suffix: "ribo", kraken_memory: "${params.kraken_memory}")
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
    // Prepare samplesheet
    if (params.read_type == "single_end") {
        samplesheet = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map{row -> tuple(row.sample, file(row.fastq))}
    } else {
        samplesheet = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map{row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
    }
    RAW(samplesheet, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)

    // Taxonomic profiling
    PROFILE(CLEAN.out.reads, kraken_db_path, params.n_reads_profile, params.ref_dir)

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

        PROCESS_OUTPUT.out.basic >> "results/qc"
        PROCESS_OUTPUT.out.adapt >> "results/qc"
        PROCESS_OUTPUT.out.qbase >> "results/qc"
        PROCESS_OUTPUT.out.qseqs >> "results/qc"

        PROFILE.out.bracken >> "results/taxonomy"
        PROFILE.out.kraken >> "results/taxonomy"
}