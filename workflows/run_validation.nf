/********************************************************
| WORKFLOW: POST-HOC VALIDATION OF PUTATIVE VIRAL READS |
********************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF_LABELED as EXTRACT_FASTQ } from "../modules/local/extractViralHitsToFastqNoref"
include { BLAST_VIRAL } from "../subworkflows/local/blastViral"
nextflow.preview.output = true

/****************
| MAIN WORKFLOW |
****************/

// Complete primary workflow
workflow RUN_VALIDATION {
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Get input FASTQ
        if ( params.viral_tsv == "" ) {
        // Option 1: Directly specify FASTQ path in config (interleaved/single-end)
            fastq_ch = Channel.fromPath(params.viral_fastq)
        } else {
        // Option 2: Extract read sequences from output DB from RUN workflow (default)
            // Define input
            tsv_ch = Channel.value(["viral_hits", file(params.viral_tsv)])
            fastq_out = EXTRACT_FASTQ(tsv_ch, params.drop_unpaired)
            fastq_ch = fastq_out.output.map { label, fastq -> fastq }
        }

        // BLAST validation on host-viral reads
        BLAST_VIRAL(fastq_ch, params.ref_dir, params.blast_db_prefix,
            params.blast_viral_fraction, params.blast_max_rank, params.blast_min_frac, params.random_seed,
            params.blast_perc_id, params.blast_qcov_hsp_perc, params.taxid_artificial)

        // Prepare results for publishing
        params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
        params_ch = Channel.of(params_str).collectFile(name: "params-run.json")
        time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
        version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
        index_params_ch = Channel.fromPath("${params.ref_dir}/input/index-params.json")
            | map { file -> file.copyTo("${params.base_dir}/work/params-index.json") }
        index_pipeline_version_ch = Channel.fromPath("${params.ref_dir}/logging/pipeline-version.txt")
            | map { file -> file.copyTo("${params.base_dir}/work/pipeline-version-index.txt") }

    emit:
        input_validation = [index_params_ch, params_ch]
        logging_validation = [index_pipeline_version_ch, time_ch, version_ch]
        output_validation = [BLAST_VIRAL.out.blast_subset, BLAST_VIRAL.out.subset_reads]
}
