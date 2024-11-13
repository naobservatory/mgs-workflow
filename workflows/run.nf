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
include { EXTRACT_VIRAL_READS } from "../subworkflows/local/extractViralReads" addParams(min_kmer_hits: "3", k: "21", bbduk_suffix: "viral", encoding: "${params.quality_encoding}", fuzzy_match: "${params.fuzzy_match_alignment_duplicates}", grouping: params.grouping)
include { BLAST_VIRAL } from "../subworkflows/local/blastViral" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")
include { PROFILE } from "../subworkflows/local/profile" addParams(min_kmer_fraction: "0.4", k: "27", bbduk_suffix: "ribo", kraken_memory: "${params.kraken_memory}", grouping: params.grouping)
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
    if ( params.grouping ) {
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
    // Prepare Kraken DB
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    // Preprocessing
    RAW(samplesheet_ch, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)
    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(CLEAN.out.reads, group_ch, params.ref_dir, kraken_db_path, params.bt2_score_threshold, params.adapters, params.host_taxon)
    // BLAST validation on host-viral reads (optional)
    if ( params.blast_hv_fraction > 0 ) {
        blast_db_path = "${params.ref_dir}/results/core_nt"
        blast_db_prefix = "core_nt"
        BLAST_VIRAL(EXTRACT_VIRAL_READS.out.fasta, blast_db_path, blast_db_prefix, params.blast_hv_fraction)
    }
    // Taxonomic profiling
    PROFILE(CLEAN.out.reads, group_ch, kraken_db_path, params.n_reads_profile, params.ref_dir)
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
