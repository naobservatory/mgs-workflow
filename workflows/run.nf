/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { COUNT_TOTAL_READS } from "../subworkflows/local/countTotalReads"
include { EXTRACT_VIRAL_READS } from "../subworkflows/local/extractViralReads"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
include { SUBSET_AND_TRIM_READS } from "../subworkflows/local/subsetAndTrimReads"
include { RUN_QC } from "../subworkflows/local/runQc"
include { PROFILE } from "../subworkflows/local/profile"
include { EXTRACT_RAW_READS_FROM_PROCESSED } from "../modules/local/extractRawReadsFromProcessed"
include { LOAD_SAMPLESHEET } from "../subworkflows/local/loadSampleSheet"
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN {
    // Setting reference paths
    kraken_db_path = "${params.ref_dir}/results/kraken_db"
    blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"

    // Load samplesheet
    LOAD_SAMPLESHEET(params.sample_sheet, params.grouping, params.single_end)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    group_ch = LOAD_SAMPLESHEET.out.group
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch)

    // Extract and count human-viral reads
    EXTRACT_VIRAL_READS(samplesheet_ch, group_ch, params.ref_dir, kraken_db_path, params.bt2_score_threshold, params.adapters, params.host_taxon, "1", "24", "viral", "${params.quality_encoding}", "${params.fuzzy_match_alignment_duplicates}", params.grouping, params.single_end)

    // Process intermediate output for chimera detection
    raw_processed_ch = EXTRACT_VIRAL_READS.out.bbduk_match.join(samplesheet_ch, by: 0)
    EXTRACT_RAW_READS_FROM_PROCESSED(raw_processed_ch, "raw_viral_subset")

    // BLAST validation on host-viral reads (optional)
    if ( params.blast_viral_fraction > 0 ) {
        BLAST_VIRAL(EXTRACT_VIRAL_READS.out.fasta, blast_db_path, params.blast_db_prefix, params.blast_viral_fraction)
        blast_subset_ch = BLAST_VIRAL.out.blast_subset
        blast_paired_ch = BLAST_VIRAL.out.blast_paired
    } else {
        blast_subset_ch = Channel.empty()
        blast_paired_ch = Channel.empty()
    }

    // Subset reads to target number, and trim adapters
    SUBSET_AND_TRIM_READS(samplesheet_ch, group_ch, params.n_reads_profile, params.grouping, params.adapters, params.single_end)

    // Run QC on subset reads before and after adapter trimming
    RUN_QC(SUBSET_AND_TRIM_READS.out.subset_reads, SUBSET_AND_TRIM_READS.out.trimmed_subset_reads, "2", "4 GB", params.single_end)

    // Profile ribosomal and non-ribosomal reads of the subset adapter-trimmed reads
    PROFILE(SUBSET_AND_TRIM_READS.out.trimmed_subset_reads, kraken_db_path, params.ref_dir, "0.4", "27", "ribo", params.single_end)

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
        EXTRACT_RAW_READS_FROM_PROCESSED.out.reads >> "reads_raw_viral"
        // QC
        COUNT_TOTAL_READS.out.read_counts >> "results"
        RUN_QC.out.qc_basic >> "results"
        RUN_QC.out.qc_adapt >> "results"
        RUN_QC.out.qc_qbase >> "results"
        RUN_QC.out.qc_qseqs >> "results"
        // Final results
        EXTRACT_VIRAL_READS.out.tsv >> "results"
        EXTRACT_VIRAL_READS.out.counts >> "results"
        PROFILE.out.bracken >> "results"
        PROFILE.out.kraken >> "results"
        // Validation output (if any)
        blast_subset_ch >> "results"
        blast_paired_ch >> "results"
}
