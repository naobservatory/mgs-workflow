/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONMIC PROFILING AND HUMAN VIRUS ANALYSIS ON SHORT-READ MGS DATA (EITHER SINGLE-END OR PAIRED-END) |
***********************************************************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_SAMPLESHEET } from "../subworkflows/local/loadSampleSheet"
include { COUNT_TOTAL_READS } from "../subworkflows/local/countTotalReads"
include { EXTRACT_VIRAL_READS_SHORT } from "../subworkflows/local/extractViralReadsShort"
include { EXTRACT_VIRAL_READS_ONT } from "../subworkflows/local/extractViralReadsONT"
include { SUBSET_TRIM } from "../subworkflows/local/subsetTrim"
include { RUN_QC } from "../subworkflows/local/runQc"
include { PROFILE} from "../subworkflows/local/profile"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
include { CHECK_VERSION_COMPATIBILITY } from "../subworkflows/local/checkVersionCompatibility"
include { COPY_FILE_BARE as COPY_INDEX_PARAMS } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_INDEX_PIPELINE_VERSION } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_INDEX_COMPAT } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_VERSION } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_PIPELINE_COMPAT } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_SAMPLESHEET } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_ADAPTERS } from "../modules/local/copyFile"

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow RUN {
    main:
    // Check index/pipeline version compatibility
    pipeline_version_path = file("${projectDir}/pipeline-version.txt")
    index_version_path = file("${params.ref_dir}/logging/pipeline-version.txt")
    pipeline_min_index_version_path = file("${projectDir}/pipeline-min-index-version.txt")
    index_min_pipeline_version_path = file("${params.ref_dir}/logging/index-min-pipeline-version.txt")
    CHECK_VERSION_COMPATIBILITY(pipeline_version_path, index_version_path,
        pipeline_min_index_version_path, index_min_pipeline_version_path)

    // Setting reference paths
    kraken_db_path = "${params.ref_dir}/results/kraken_db"

    // Load samplesheet and check platform
    LOAD_SAMPLESHEET(params.sample_sheet, params.platform, false)
    samplesheet_ch = LOAD_SAMPLESHEET.out.samplesheet
    start_time_str = LOAD_SAMPLESHEET.out.start_time_str
    single_end_ch = LOAD_SAMPLESHEET.out.single_end

    // Count reads in files
    COUNT_TOTAL_READS(samplesheet_ch, single_end_ch)

    // Extract and count human-viral reads
    if ( params.platform == "ont" ) {
        ont_params = [
            ref_dir: params.ref_dir,
            taxid_artificial: params.taxid_artificial
        ]
        EXTRACT_VIRAL_READS_ONT(samplesheet_ch, ont_params)
        hits_fastq = EXTRACT_VIRAL_READS_ONT.out.hits_fastq
        hits_final = EXTRACT_VIRAL_READS_ONT.out.hits_final
        inter_lca = EXTRACT_VIRAL_READS_ONT.out.inter_lca
        inter_aligner = EXTRACT_VIRAL_READS_ONT.out.inter_minimap2
        bbduk_match = Channel.empty()
        bbduk_trimmed = Channel.empty()
     } else {
        short_params = [
            ref_dir: params.ref_dir,
            aln_score_threshold: params.bt2_score_threshold,
            adapter_path: params.adapters,
            cutadapt_error_rate: params.cutadapt_error_rate,
            min_kmer_hits: "1",
            k: "24",
            bbduk_suffix: "viral",
            taxid_artificial: params.taxid_artificial
        ]
        EXTRACT_VIRAL_READS_SHORT(samplesheet_ch, short_params)
        hits_fastq = EXTRACT_VIRAL_READS_SHORT.out.hits_fastq
        hits_final = EXTRACT_VIRAL_READS_SHORT.out.hits_final
        inter_lca = EXTRACT_VIRAL_READS_SHORT.out.inter_lca
        inter_aligner = EXTRACT_VIRAL_READS_SHORT.out.inter_bowtie
        bbduk_match = EXTRACT_VIRAL_READS_SHORT.out.bbduk_match
        bbduk_trimmed = EXTRACT_VIRAL_READS_SHORT.out.bbduk_trimmed
    }
    // BLAST validation on host-viral reads (optional)
    if ( params.blast_viral_fraction > 0 ) {
        blast_viral_params = [
            ref_dir: params.ref_dir,
            blast_db_prefix: params.blast_db_prefix,
            read_fraction: params.blast_viral_fraction,
            blast_max_rank: params.blast_max_rank,
            blast_min_frac: params.blast_min_frac,
            random_seed: params.random_seed,
            perc_id: params.blast_perc_id,
            qcov_hsp_perc: params.blast_qcov_hsp_perc,
            taxid_artificial: params.taxid_artificial
        ]
        BLAST_VIRAL(hits_fastq, blast_viral_params)
        blast_subset_ch = BLAST_VIRAL.out.blast_subset
        blast_reads_ch = BLAST_VIRAL.out.subset_reads
    } else {
        blast_subset_ch = Channel.empty()
        blast_reads_ch = Channel.empty()
    }

    // Subset reads to target number, and trim adapters
    subset_trim_params = [
        n_reads: params.n_reads_profile,
        adapter_path: params.adapters,
        platform: params.platform,
        random_seed: params.random_seed
    ]
    SUBSET_TRIM(samplesheet_ch, subset_trim_params, single_end_ch)

    // Run QC on subset reads before and after adapter trimming
    RUN_QC(SUBSET_TRIM.out.subset_reads, SUBSET_TRIM.out.trimmed_subset_reads, single_end_ch)

    // Profile ribosomal and non-ribosomal reads of the subset adapter-trimmed reads
    profile_params = [
        ref_dir: params.ref_dir,
        min_kmer_fraction: "0.4",
        k: "27",
        ribo_suffix: "ribo",
        bracken_threshold: params.bracken_threshold,
        platform: params.platform
    ]
    PROFILE(SUBSET_TRIM.out.trimmed_subset_reads, kraken_db_path, profile_params, single_end_ch)

    // Get index files for publishing
    index_params_path = file("${params.ref_dir}/input/index-params.json")
    index_params_ch = COPY_INDEX_PARAMS(Channel.fromPath(index_params_path), "params-index.json")
    index_pipeline_version_ch = COPY_INDEX_PIPELINE_VERSION(Channel.fromPath(index_version_path), "pipeline-version-index.txt")
    index_min_pipeline_version_newpath = index_min_pipeline_version_path.getFileName().toString()
    index_compatibility_ch = COPY_INDEX_COMPAT(Channel.fromPath(index_min_pipeline_version_path), index_min_pipeline_version_newpath)

    // Prepare other publishing variables
    params_str = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "params-run.json")
    time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
    pipeline_version_newpath = pipeline_version_path.getFileName().toString()
    // Note: we send these input/logging files through a COPY_FILE_BASE process
    // because nextflow 25.04 now only publishes files that have passed through the working directory.
    // We first tried collectFile() as an alternative; however it intermittantly gives serialization errors.
    version_ch = COPY_VERSION(Channel.fromPath(pipeline_version_path), pipeline_version_newpath)
    pipeline_min_index_version_newpath = pipeline_min_index_version_path.getFileName().toString()
    pipeline_compatibility_ch = COPY_PIPELINE_COMPAT(Channel.fromPath(pipeline_min_index_version_path), pipeline_min_index_version_newpath)
    samplesheet_ch = COPY_SAMPLESHEET(Channel.fromPath(params.sample_sheet), "samplesheet.csv")
    adapters_ch = COPY_ADAPTERS(Channel.fromPath(params.adapters), "adapters.fasta")

    emit:
        input_run = index_params_ch.mix(samplesheet_ch, adapters_ch, params_ch)
        logging_run = index_pipeline_version_ch.mix(index_compatibility_ch, time_ch, version_ch, pipeline_compatibility_ch)
        intermediates_run = hits_fastq.mix(inter_lca, inter_aligner)
        reads_raw_viral = bbduk_match
        reads_trimmed_viral = bbduk_trimmed
        // Lots of results; split across 2 channels (QC, and other)
        qc_results_run = COUNT_TOTAL_READS.out.read_counts.mix(RUN_QC.out.qc_basic, RUN_QC.out.qc_adapt, RUN_QC.out.qc_qbase, RUN_QC.out.qc_qseqs, RUN_QC.out.qc_lengths)
        other_results_run = hits_final.mix(PROFILE.out.bracken, PROFILE.out.kraken, blast_subset_ch, blast_reads_ch)
}
