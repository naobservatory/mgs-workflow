/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { FASTQC_LABELED } from "../../../modules/local/fastqc" addParams(cpus: "${params.fastqc_cpus}", mem: "${params.fastqc_mem}")
include { MULTIQC_LABELED } from "../../../modules/local/multiqc"
include { SUMMARIZE_MULTIQC_PAIR } from "../../../modules/local/summarizeMultiqcPair"
include { MERGE_TSVS as MERGE_MULTIQC_BASIC } from "../../../modules/local/mergeTsvs" addParams(name: "${params.stage_label}_qc_basic_stats")
include { MERGE_TSVS as MERGE_MULTIQC_ADAPT } from "../../../modules/local/mergeTsvs" addParams(name: "${params.stage_label}_qc_adapter_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QBASE } from "../../../modules/local/mergeTsvs" addParams(name: "${params.stage_label}_qc_quality_base_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QSEQS } from "../../../modules/local/mergeTsvs" addParams(name: "${params.stage_label}_qc_quality_sequence_stats")

/***********
| WORKFLOW |
***********/

workflow QC {
    take:
        reads
    main:
        // 1. Run FASTQC on each pair of read files
        fastqc_ch = FASTQC_LABELED(reads)
        // 2. Extract data with MultiQC for each pair of read files
        multiqc_ch = MULTIQC_LABELED(params.stage_label, fastqc_ch.zip)
        // 3. Summarize MultiQC information for each pair of read files
        process_ch = SUMMARIZE_MULTIQC_PAIR(multiqc_ch.data)
        // 4. Collate MultiQC outputs
        multiqc_basic_ch = process_ch.map{ it[0] }.collect().ifEmpty([])
        multiqc_adapt_ch = process_ch.map{ it[1] }.collect().ifEmpty([])
        multiqc_qbase_ch = process_ch.map{ it[2] }.collect().ifEmpty([])
        multiqc_qseqs_ch = process_ch.map{ it[3] }.collect().ifEmpty([])
        // 5. Merge MultiQC outputs
        basic_out_ch = MERGE_MULTIQC_BASIC(multiqc_basic_ch)
        adapt_out_ch = MERGE_MULTIQC_ADAPT(multiqc_adapt_ch)
        qbase_out_ch = MERGE_MULTIQC_QBASE(multiqc_qbase_ch)
        qseqs_out_ch = MERGE_MULTIQC_QSEQS(multiqc_qseqs_ch)
        // 6. Combine outputs into a single output channel
        out_ch = basic_out_ch.combine(adapt_out_ch)
            .combine(qbase_out_ch).combine(qseqs_out_ch)
            .map({file1, file2, file3, file4 -> tuple(file1, file2, file3, file4)})
    emit:
        qc = out_ch
}
