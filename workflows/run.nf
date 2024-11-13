/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
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
workflow RUN {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Prepare samplesheet
    samplesheet = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map{row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
    // Prepare Kraken DB
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    // Preprocessing
    RAW(samplesheet, params.n_reads_trunc, "2", "4 GB", "raw_concat")
    CLEAN(RAW.out.reads, params.adapters, "2", "4 GB", "cleaned")
    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(CLEAN.out.reads, params.ref_dir, kraken_db_path, params.bt2_score_threshold, params.adapters, params.host_taxon, "3", "21", "viral", "${params.quality_encoding}", "${params.fuzzy_match_alignment_duplicates}", "${params.kraken_memory}")
    // BLAST validation on host-viral reads (optional)
    if ( params.blast_viral_fraction > 0 ) {
        blast_db_path = "${params.ref_dir}/results/core_nt"
        blast_db_prefix = "core_nt"
        BLAST_VIRAL(EXTRACT_VIRAL_READS.out.fasta, blast_db_path, blast_db_prefix, params.blast_viral_fraction, "32", "256 GB", "32 GB")
    }
    // Taxonomic profiling
    PROFILE(CLEAN.out.reads, kraken_db_path, params.n_reads_profile, params.ref_dir, "0.4", "27", "ribo", "${params.kraken_memory}")
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
        Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt").collectFile(name: "pipeline-version-index.txt") >> "logging"
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
        EXTRACT_VIRAL_READS.out.tsv >> "results"
        EXTRACT_VIRAL_READS.out.counts >> "results"
        PROFILE.out.bracken >> "results"
        PROFILE.out.kraken >> "results"
}
