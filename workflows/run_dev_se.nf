/***********************************************************************************************
| WORKFLOW: PREPROCESSING ON SHORT-READ MGS DATA (EITHER SINGLE-END OR PAIRED-END) |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_SAMPLESHEET } from "../subworkflows/local/loadSampleSheet"
include { COUNT_TOTAL_READS } from "../subworkflows/local/countTotalReads"
include { SUBSET_TRIM } from "../subworkflows/local/subsetTrim"
include { RUN_QC } from "../subworkflows/local/runQc"
include { PROFILE } from "../subworkflows/local/profile"
include { EXTRACT_VIRAL_READS_ONT } from "../subworkflows/local/extractViralReadsONT"
include { EXTRACT_VIRAL_READS_SHORT } from "../subworkflows/local/extractViralReadsShort"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
include { CHECK_VERSION_COMPATIBILITY } from "../subworkflows/local/checkVersionCompatibility"

nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_DEV_SE {
    // Check index/pipeline version compatibility
    pipeline_version_path = "${projectDir}/pipeline-version.txt"
    index_version_path = "${params.ref_dir}/logging/pipeline-version.txt"
    pipeline_min_index_version_path = "${projectDir}/pipeline-min-index-version.txt"
    index_min_pipeline_version_path = "${params.ref_dir}/logging/index-min-pipeline-version.txt"
    CHECK_VERSION_COMPATIBILITY(pipeline_version_path, index_version_path,
        pipeline_min_index_version_path, index_min_pipeline_version_path)

    // Setting reference paths
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"

    // Load samplesheet
    LOAD_SAMPLESHEET(params.sample_sheet)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str
    single_end_ch = LOAD_SAMPLESHEET.out.single_end

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch, single_end_ch)

    // Extract viral reads
    if ( params.ont ) {
        EXTRACT_VIRAL_READS_ONT(samplesheet_ch, params.ref_dir)
        hv_tsv_ch = EXTRACT_VIRAL_READS_ONT.out.hits_hv
        hv_fastqs = EXTRACT_VIRAL_READS_ONT.out.hits_fastq
    } else {
        hv_tsv_ch = Channel.empty()
        hv_fastqs = Channel.empty()
    }

    // BLAST viral reads
    if ( params.blast_viral_fraction > 0 ) {
        BLAST_VIRAL(hv_fastqs, blast_db_path, params.blast_db_prefix,
            params.blast_viral_fraction, params.blast_max_rank, params.blast_min_frac, params.random_seed,
            params.blast_perc_id, params.blast_qcov_hsp_perc)
        blast_subset_ch = BLAST_VIRAL.out.blast_subset
        blast_reads_ch = BLAST_VIRAL.out.subset_reads
    } else {
        blast_subset_ch = Channel.empty()
        blast_reads_ch = Channel.empty()
    }

    // Subset reads to target number, and trim adapters
    SUBSET_TRIM(samplesheet_ch, params.n_reads_profile,
        params.adapters, single_end_ch,
        params.ont, params.random_seed)

    // Run QC on subset reads before and after adapter trimming
    RUN_QC(SUBSET_TRIM.out.subset_reads, SUBSET_TRIM.out.trimmed_subset_reads, single_end_ch)

    // Profile ribosomal and non-ribosomal reads of the subset adapter-trimmed reads
    if ( params.ont ) {
        bracken_ch = Channel.empty()
        kraken_ch = Channel.empty()
    } else {
        PROFILE(SUBSET_TRIM.out.trimmed_subset_reads, kraken_db_path, params.ref_dir, "0.4", "27", "ribo",
            params.bracken_threshold, single_end_ch)
        bracken_ch = PROFILE.out.bracken
        kraken_ch = PROFILE.out.kraken
    }

    // Get index files for publishing
    index_params_path = "${params.ref_dir}/input/index-params.json"
    index_params_path_new = "${params.base_dir}/work/params-index.json"
    index_version_path_new = "${params.base_dir}/work/pipeline-version-index.txt"
    index_min_pipeline_version_path_new = "${params.base_dir}/work/index-min-pipeline-version.txt"
    index_params_ch = Channel.fromPath(index_params_path)
        | map { file -> file.copyTo(index_params_path_new) }
    index_pipeline_version_ch = Channel.fromPath(index_version_path)
        | map { file -> file.copyTo(index_version_path_new) }
    index_compatibility_ch = Channel.fromPath(index_min_pipeline_version_path)
        | map { file -> file.copyTo(index_min_pipeline_version_path_new) }

    // Prepare other publishing variables
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
    version_ch = Channel.fromPath(pipeline_version_path)
    pipeline_compatibility_ch = Channel.fromPath(pipeline_min_index_version_path)

    // Publish results
    publish:
        // Saved inputs
        index_params_ch >> "input"
        index_pipeline_version_ch >> "logging"
        index_compatibility_ch >> "logging"
        Channel.fromPath(params.sample_sheet) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        pipeline_compatibility_ch >> "logging"
        // QC
        COUNT_TOTAL_READS.out.read_counts >> "results"
        RUN_QC.out.qc_basic >> "results"
        RUN_QC.out.qc_adapt >> "results"
        RUN_QC.out.qc_qbase >> "results"
        RUN_QC.out.qc_qseqs >> "results"
        RUN_QC.out.qc_lengths >> "results"
        // Final results
        hv_tsv_ch >> "results"
        hv_fastqs >> "intermediates"
        bracken_ch >> "results"
        kraken_ch >> "results"
        blast_subset_ch >> "results"
        blast_reads_ch >> "results"

}
