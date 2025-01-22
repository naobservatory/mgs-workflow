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
include { SUBSET_TRIM_STREAMED as SUBSET_TRIM } from "../subworkflows/local/subsetTrimStreamed"
include { RUN_QC_STREAMED as RUN_QC } from "../subworkflows/local/runQcStreamed"
include { PROFILE_STREAMED as PROFILE } from "../subworkflows/local/profileStreamed"
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
    // TODO: Drop grouping and update subworkflow
    LOAD_SAMPLESHEET(params.sample_sheet, params.grouping, params.single_end)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    group_ch = LOAD_SAMPLESHEET.out.group
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch, params.single_end)

    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(samplesheet_ch, params.ref_dir, kraken_db_path,
        params.bt2_score_threshold, params.adapters, params.host_taxon,
        "1", "24", "viral", params.bracken_threshold)

    // TODO: Add BLAST validation

    // Subset reads to target number, and trim adapters
    SUBSET_TRIM(samplesheet_ch, params.n_reads_profile,
        params.adapters, params.single_end)

    // Run QC on subset reads before and after adapter trimming (NB: unchanged in streamed version)
    RUN_QC(SUBSET_TRIM.out.subset_reads, SUBSET_TRIM.out.trimmed_subset_reads, params.single_end)

    // Profile ribosomal and non-ribosomal reads of the subset adapter-trimmed reads
    PROFILE(SUBSET_TRIM.out.trimmed_subset_reads, kraken_db_path, params.ref_dir, "0.4", "27", "ribo",
        params.bracken_threshold, params.single_end)

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
        // Intermediate files
        EXTRACT_VIRAL_READS.out.bbduk_match >> "reads_raw_viral"
        EXTRACT_VIRAL_READS.out.hits_all    >> "intermediates"
        EXTRACT_VIRAL_READS.out.hits_fastq  >> "intermediates"
        // QC
        COUNT_TOTAL_READS.out.read_counts >> "results"
        RUN_QC.out.qc_basic >> "results"
        RUN_QC.out.qc_adapt >> "results"
        RUN_QC.out.qc_qbase >> "results"
        RUN_QC.out.qc_qseqs >> "results"
        // Final results
        EXTRACT_VIRAL_READS.out.hits_filtered >> "results"
        PROFILE.out.bracken >> "results"
        PROFILE.out.kraken >> "results"
}
