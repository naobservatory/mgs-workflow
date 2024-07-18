/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { HV_SCREEN_1 } from "../subworkflows/local/hv_screen"
include { HV_SCREEN_2 } from "../subworkflows/local/hv_screen" addParams(adapter_path: params.adapters)
include { HV_SCREEN_3 } from "../subworkflows/local/hv_screen" addParams(min_kmer_fraction: "0.1", k: "21")
include { HV_SCREEN_4 } from "../subworkflows/local/hv_screen" addParams(min_kmer_fraction: "0.1", k: "21", adapter_path: params.adapters)
include { HV_SCREEN_5 } from "../subworkflows/local/hv_screen" addParams(min_kmer_fraction: "0.1", k: "21", adapter_path: params.adapters)
include { BBDUK_MASK } from "../modules/local/bbduk" addParams(label: "human-viral-genomes-filtered")
include { CUTADAPT_MASK } from "../modules/local/cutadapt" addParams(label: "human-viral-genomes-filtered")
include { BBMASK } from "../modules/local/bbmask" addParams(label: "human-viral-genomes-filtered")
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN2 {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Prepare libraries
    libraries_ch = Channel
        .fromPath(params.library_tab)
        .splitCsv(header: true)
        .map{row -> [row.sample, row.library]}
        .groupTuple()
    // Prepare Kraken DB
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    hv_ref_path = "${params.ref_dir}/results/human-viral-genomes-filtered.fasta.gz"
    // Mask HV database for BBDuk
    BBDUK_MASK(hv_ref_path, params.adapters, 11)
    CUTADAPT_MASK(BBDUK_MASK.out.masked, params.adapters)
    BBMASK(CUTADAPT_MASK.out.masked)
    // Prepare raw data
    RAW(libraries_ch, params.raw_dir, params.n_reads_trunc)
    // Compare HV screening methods
    HV_SCREEN_1(RAW.out.reads, params.ref_dir)
    HV_SCREEN_2(RAW.out.reads, params.ref_dir)
    HV_SCREEN_3(RAW.out.reads, params.ref_dir)
    HV_SCREEN_4(RAW.out.reads, params.ref_dir)
    HV_SCREEN_5(RAW.out.reads, BBMASK.out.masked)
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
        // HV Screen Outputs
        HV_SCREEN_1.out.fastq >> "results/hv_screen"
        HV_SCREEN_2.out.fastq >> "results/hv_screen"
        HV_SCREEN_3.out.fastq >> "results/hv_screen"
        HV_SCREEN_3.out.stats >> "results/hv_screen"
        HV_SCREEN_4.out.fastq >> "results/hv_screen"
        HV_SCREEN_4.out.stats >> "results/hv_screen"
        HV_SCREEN_5.out.fastq >> "results/hv_screen"
        HV_SCREEN_5.out.stats >> "results/hv_screen"
}
