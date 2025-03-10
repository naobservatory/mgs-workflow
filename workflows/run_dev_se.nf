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


nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN_DEV_SE {
    // Setting reference paths
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"

    // Load samplesheet
    LOAD_SAMPLESHEET(params.sample_sheet, params.single_end)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch, params.single_end)

    // Extract viral reads
    if ( params.ont ) {
        EXTRACT_VIRAL_READS_ONT(samplesheet_ch, params.ref_dir, params.host_taxon)
        hv_tsv_ch = EXTRACT_VIRAL_READS_ONT.out.hv_tsv
    } else {
        EXTRACT_VIRAL_READS_SHORT(samplesheet_ch, params.ref_dir, kraken_db_path, params.bt2_score_threshold, params.adapters, params.host_taxon, "1", "24", "viral", params.bracken_threshold)
        hv_tsv_ch = EXTRACT_VIRAL_READS_SHORT.out.hv_tsv
    }

    // Subset reads to target number, and trim adapters
    SUBSET_TRIM(samplesheet_ch, params.n_reads_profile,
        params.adapters, params.single_end,
        params.ont, params.random_seed)

    // Run QC on subset reads before and after adapter trimming
    RUN_QC(SUBSET_TRIM.out.subset_reads, SUBSET_TRIM.out.trimmed_subset_reads, params.single_end)

    // Profile ribosomal and non-ribosomal reads of the subset adapter-trimmed reads
    if ( params.ont ) {
        bracken_ch = Channel.empty()
        kraken_ch = Channel.empty()
    } else {
        PROFILE(SUBSET_TRIM.out.trimmed_subset_reads, kraken_db_path, params.ref_dir, "0.4", "27", "ribo",
            params.bracken_threshold, params.single_end)
        bracken_ch = PROFILE.out.bracken
        kraken_ch = PROFILE.out.kraken
    }

    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    index_params_ch = Channel.fromPath("${params.ref_dir}/input/index-params.json")
        | map { file -> file.copyTo("${params.base_dir}/work/params-index.json") }
    index_pipeline_version_ch = Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt")
        | map { file -> file.copyTo("${params.base_dir}/work/pipeline-version-index.txt") }
    publish:
        // Saved inputs
        index_params_ch >> "input"
        index_pipeline_version_ch >> "logging"
        Channel.fromPath(params.sample_sheet) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // QC
        COUNT_TOTAL_READS.out.read_counts >> "results"
        RUN_QC.out.qc_basic >> "results"
        RUN_QC.out.qc_adapt >> "results"
        RUN_QC.out.qc_qbase >> "results"
        RUN_QC.out.qc_qseqs >> "results"
        RUN_QC.out.qc_lengths >> "results"
        // Final results
        hv_tsv_ch >> "results"
        bracken_ch >> "results"
        kraken_ch >> "results"
        // Viral BLAST results
        // EXTRACT_VIRAL_READS.out.blast_subset >> "results"
        // EXTRACT_VIRAL_READS.out.blast_reads >> "results"
}
