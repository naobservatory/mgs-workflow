/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { CLEAN } from "../subworkflows/local/clean" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "cleaned")
include { DEDUP } from "../subworkflows/local/dedup" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "dedup")
include { RIBODEPLETION } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribodepletion", min_kmer_fraction: "0.4", k: "27", bbduk_suffix: "ribodepletion")
include { TAXONOMY } from "../subworkflows/local/taxonomy"
include { PROCESS_OUTPUT } from "../subworkflows/local/processOutput"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Prepare libraries
    libraries_ch = Channel
        .fromPath(params.library_tab)
        .splitCsv(header: true)
        .map{row -> [row.sample, row.library]}
        .groupTuple()
    // Preprocessing
    RAW(libraries_ch, params.raw_dir, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)
    DEDUP(CLEAN.out.reads)
    RIBODEPLETION(DEDUP.out.reads, params.ref_dir)
    // merge paired reads, dedup considering RC
    TAXONOMY(RIBODEPLETION.out.reads)
    // Process output
    qc_ch = RAW.out.qc.concat(CLEAN.out.qc, DEDUP.out.qc, RIBODEPLETION.out.qc)
    PROCESS_OUTPUT(qc_ch)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        Channel.fromPath("${params.ref_dir}/input/index-params.json") >> "input"
        Channel.fromPath("${params.ref_dir}/input/pipeline-version.txt").collectFile(name: "pipeline-version-index.txt") >> "input"
        Channel.fromPath(params.sample_tab) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "input"
        version_ch >> "input"
        // Intermediate files
        CLEAN.out.reads >> "intermediates/reads/cleaned"
        RIBODEPLETION.out.reads >> "intermediates/reads/ribodepleted"
        // final joined-deduped
        TAXONOMY.out.joined_reads >> "results/reads"
        // QC
        PROCESS_OUTPUT.out.basic >> "results/qc"
        PROCESS_OUTPUT.out.adapt >> "results/qc"
        PROCESS_OUTPUT.out.qbase >> "results/qc"
        PROCESS_OUTPUT.out.qseqs >> "results/qc"
}
