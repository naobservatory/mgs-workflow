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
include { EXTRACT_RAW_READS_FROM_PROCESSED } from "../modules/local/extractRawReadsFromProcessed"
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
    // Preprocessing
    RAW(samplesheet_ch, params.n_reads_trunc, "2", "4 GB", "raw_concat")
    CLEAN(RAW.out.reads, params.adapters, "2", "4 GB", "cleaned")
    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(CLEAN.out.reads, group_ch, params.ref_dir, params.kraken_db_path, params.bt2_score_threshold, params.adapters, params.host_taxon, "3", "21", "viral", "${params.quality_encoding}", "${params.fuzzy_match_alignment_duplicates}", params.grouping)
    // Process intermediate output for chimera detection
    raw_processed_ch = EXTRACT_VIRAL_READS.out.bbduk_match.join(RAW.out.reads, by: 0)
    EXTRACT_RAW_READS_FROM_PROCESSED(raw_processed_ch, "raw_viral_subset")
    // BLAST validation on host-viral reads (optional)
    if ( params.blast_viral_fraction > 0 ) {
        BLAST_VIRAL(EXTRACT_VIRAL_READS.out.fasta, params.blast_db_path, params.blast_db_prefix, params.blast_viral_fraction)
    }
    // Taxonomic profiling
    PROFILE(CLEAN.out.reads, group_ch, params.kraken_db_path, params.n_reads_profile, params.ref_dir, "0.4", "27", "ribo", params.grouping)
    // Process output
    qc_ch = RAW.out.qc.concat(CLEAN.out.qc)
    PROCESS_OUTPUT(qc_ch)
    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
//    index_params_ch = Channel.fromPath("${params.ref_dir}/input/index-params.json")
    index_pipeline_version_ch = Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt")
    .map { file -> file.copyTo("pipeline-version-index.txt") }
    publish:
        // Saved inputs
//        index_params_ch >> "input"
        index_pipeline_version_ch >> "logging"
        Channel.fromPath(params.sample_sheet) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // Intermediate files
        CLEAN.out.reads >> "intermediates/reads/cleaned"
        EXTRACT_RAW_READS_FROM_PROCESSED.out.reads >> "intermediates/reads/raw_viral"
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
