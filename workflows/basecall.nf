/***********************************************************************************************
| WORKFLOW: BASECALLING NANOPORE SQUIGGLE DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BATCH_POD_5 } from "../modules/local/batchPod5"
include { BASECALL_POD_5 } from "../modules/local/dorado"
include { DEMUX_POD_5 } from "../modules/local/dorado"
include { BAM_TO_FASTQ } from "../modules/local/samtools"

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow BASECALL {
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Batching
        // batch_pod5_ch = BATCH_POD_5(params.pod_5_dir, params.batch_size)
        //    .flatten()

        // Basecalling
        bam_ch = BASECALL_POD_5(params.pod_5_dir, params.kit)

        // Demultiplexing
        demux_ch = DEMUX_POD_5(bam_ch.bam, params.kit)

        // Convert to FASTQ
        fastq_ch = BAM_TO_FASTQ(demux_ch.demux_bam, params.nanopore_run)



    publish:
        fastq_ch >> "raw"
        bam_ch.summary >> "summary"
}