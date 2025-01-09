/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_BASIC } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_ADAPT } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QBASE } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QSEQS } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow PROCESS_OUTPUT {
    take:
        multiqc_ch
    main:
        // Collate MultiQC outputs
        multiqc_basic_ch = multiqc_ch.map{ it[0] }.collect().ifEmpty([])
        multiqc_adapt_ch = multiqc_ch.map{ it[1] }.collect().ifEmpty([])
        multiqc_qbase_ch = multiqc_ch.map{ it[2] }.collect().ifEmpty([])
        multiqc_qseqs_ch = multiqc_ch.map{ it[3] }.collect().ifEmpty([])
        // Merge MultiQC outputs
        basic_out_ch = CONCATENATE_MULTIQC_BASIC(multiqc_basic_ch, "qc_basic_stats")
        adapt_out_ch = CONCATENATE_MULTIQC_ADAPT(multiqc_adapt_ch, "qc_adapter_stats")
        qbase_out_ch = CONCATENATE_MULTIQC_QBASE(multiqc_qbase_ch, "qc_quality_base_stats")
        qseqs_out_ch = CONCATENATE_MULTIQC_QSEQS(multiqc_qseqs_ch, "qc_quality_sequence_stats")
    emit:
        basic = basic_out_ch
        adapt = adapt_out_ch
        qbase = qbase_out_ch
        qseqs = qseqs_out_ch
}
