/***********************************************************************************************
| WORKFLOW: BASECALLING NANOPORE SQUIGGLE DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BASECALL_POD_5} from "../modules/local/dorado"
include { BAM_TO_FASTQ} from "../modules/local/samtools"
/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow BASECALL {
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Basecalling
        // TODO: Could parallelise batching in future
        pod_5_dir = params.pod5_dir
        // calls_bam = params.calls_bam

        bam_ch = BASECALL_POD_5(pod_5_dir)

        fastq_ch = BAM_TO_FASTQ(bam_ch, params.nanopore_run)

    publish:
        fastq_ch >> "raw/${params.nanopore_run}.fastq.gz"
}