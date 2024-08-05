/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MERGE_TSVS as MERGE_MULTIQC_BASIC } from "../../../modules/local/mergeTsvs" addParams(name: "qc_basic_stats")
include { MERGE_TSVS as MERGE_MULTIQC_ADAPT } from "../../../modules/local/mergeTsvs" addParams(name: "qc_adapter_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QBASE } from "../../../modules/local/mergeTsvs" addParams(name: "qc_quality_base_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QSEQS } from "../../../modules/local/mergeTsvs" addParams(name: "qc_quality_sequence_stats")

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
        basic_out_ch = MERGE_MULTIQC_BASIC(multiqc_basic_ch)
        adapt_out_ch = MERGE_MULTIQC_ADAPT(multiqc_adapt_ch)
        qbase_out_ch = MERGE_MULTIQC_QBASE(multiqc_qbase_ch)
        qseqs_out_ch = MERGE_MULTIQC_QSEQS(multiqc_qseqs_ch)
    emit:
        basic = basic_out_ch
        adapt = adapt_out_ch
        qbase = qbase_out_ch
        qseqs = qseqs_out_ch
}
