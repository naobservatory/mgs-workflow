/***********************************************************************************************
| WORKFLOW: PREPROCESSING, TAXONOMIC PROFILING AND HUMAN-VIRUS ANALYSIS ON SHORT-READ MGS DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BLAST_HV } from "../subworkflows/local/blastHV" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")
nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow BLAST_FLOW {
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")
    // Prepare samplesheet
    samplesheet = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map{row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}

    concat_ch = samplesheet.map { sample, read1, read2 ->
            tuple(sample, [read1, read2])
        }

    // BLAST validation on human-viral reads (optional)
    if ( params.blast_hv_fraction > 0 ) {
        blast_nt_path = "${params.ref_dir}/results/nt"
        BLAST_HV(concat_ch, blast_nt_path, params.blast_hv_fraction)
    }

    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "run-params.json")
    time_ch = Channel.of(start_time_str + "\n").collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        Channel.fromPath("${params.ref_dir}/input/index-params.json") >> "input"
        Channel.fromPath("${params.ref_dir}/input/pipeline-version.txt").collectFile(name: "pipeline-version-index.txt") >> "logging"
        Channel.fromPath(params.sample_sheet) >> "input"
        Channel.fromPath(params.adapters) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
}
