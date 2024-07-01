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
include { RIBODEPLETION as RIBO_INITIAL } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_initial", min_kmer_fraction: "0.6", k: "43", bbduk_suffix: "ribo_initial")
include { RIBODEPLETION as RIBO_SECONDARY } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_secondary", min_kmer_fraction: "0.4", k: "27", bbduk_suffix: "ribo_secondary")
include { TAXONOMY as TAXONOMY_FULL } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: 1, kraken_memory: "${params.kraken_memory}")
include { TAXONOMY as TAXONOMY_PRE } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: params.classify_dedup_subset, kraken_memory: "${params.kraken_memory}")
include { TAXONOMY as TAXONOMY_POST } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: params.classify_dedup_subset, kraken_memory: "${params.kraken_memory}")
include { HV } from "../subworkflows/local/hv"
include { BLAST_HV } from "../subworkflows/local/blastHV" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")
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
    // Prepare Kraken DB
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    // Preprocessing
    RAW(libraries_ch, params.raw_dir, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)
    DEDUP(CLEAN.out.reads)
    RIBO_INITIAL(DEDUP.out.reads, params.ref_dir)
    RIBO_SECONDARY(RIBO_INITIAL.out.reads, params.ref_dir)
    // Taxonomic profiling (all ribodepleted)
    TAXONOMY_FULL(RIBO_SECONDARY.out.reads, kraken_db_path)
    // Taxonomic profiling (pre/post dedup)
    TAXONOMY_PRE(CLEAN.out.reads, kraken_db_path)
    TAXONOMY_POST(DEDUP.out.reads, kraken_db_path)
    // Extract and count human-viral reads
    HV(RIBO_INITIAL.out.reads, params.ref_dir, kraken_db_path, params.bt2_score_threshold)
    // BLAST validation on human-viral reads (optional)
    if ( params.blast_hv_fraction > 0 ) {
        blast_nt_path = "${params.ref_dir}/results/nt"
        BLAST_HV(HV.out.fasta, blast_nt_path, params.blast_hv_fraction)
    }
    // Process output
    qc_ch = RAW.out.qc.concat(CLEAN.out.qc, DEDUP.out.qc, RIBO_INITIAL.out.qc, RIBO_SECONDARY.out.qc)
    PROCESS_OUTPUT(qc_ch, TAXONOMY_FULL.out.bracken, TAXONOMY_PRE.out.bracken, TAXONOMY_POST.out.bracken, params.classify_dedup_subset)
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
        // QC
        PROCESS_OUTPUT.out.basic >> "results/qc"
        PROCESS_OUTPUT.out.adapt >> "results/qc"
        PROCESS_OUTPUT.out.qbase >> "results/qc"
        PROCESS_OUTPUT.out.qseqs >> "results/qc"
        // Final results
        HV.out.tsv >> "results/hv"
        HV.out.counts >> "results/hv"
        PROCESS_OUTPUT.out.composition_full >> "results/taxonomy_final"
        PROCESS_OUTPUT.out.composition_pre >> "results/taxonomy_pre_dedup"
        PROCESS_OUTPUT.out.composition_post >> "results/taxonomy_post_dedup"
        TAXONOMY_FULL.out.kraken_reports >> "results/taxonomy_final"
        TAXONOMY_PRE.out.kraken_reports >> "results/taxonomy_pre_dedup"
        TAXONOMY_POST.out.kraken_reports >> "results/taxonomy_post_dedup"
}
