/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { FASTQC_LABELED } from "../../../modules/local/fastqc"
include { MULTIQC_LABELED } from "../../../modules/local/multiqc"
include { SUMMARIZE_MULTIQC } from "../../../modules/local/summarizeMultiqc"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_BASIC } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_ADAPT } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QBASE } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QSEQS } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_LENGTHS } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow QC {
    take:
        reads
        fastqc_cpus
        fastqc_mem
        stage_label
        single_end
    main:
        // 1. Run FASTQC on each read file / pair of read files
        fastqc_ch = FASTQC_LABELED(reads, fastqc_cpus, fastqc_mem)
        // 2. Extract data with MultiQC for each read file / pair of read files
        multiqc_ch = MULTIQC_LABELED(stage_label, fastqc_ch.zip)
        // 3. Summarize MultiQC information for each read file / pair of read files
        process_ch = SUMMARIZE_MULTIQC(multiqc_ch.data, single_end)
        // 4. Collate MultiQC outputs
        multiqc_basic_ch = process_ch.map{ it[0] }.collect().ifEmpty([])
        multiqc_adapt_ch = process_ch.map{ it[1] }.collect().ifEmpty([])
        multiqc_qbase_ch = process_ch.map{ it[2] }.collect().ifEmpty([])
        multiqc_qseqs_ch = process_ch.map{ it[3] }.collect().ifEmpty([])
        multiqc_lengths_ch = process_ch.map{ it[4] }.collect().ifEmpty([])
        // 5. Merge MultiQC outputs
        basic_out_ch = CONCATENATE_MULTIQC_BASIC(multiqc_basic_ch, "${stage_label}_subset_qc_basic_stats")
        adapt_out_ch = CONCATENATE_MULTIQC_ADAPT(multiqc_adapt_ch, "${stage_label}_subset_qc_adapter_stats")
        qbase_out_ch = CONCATENATE_MULTIQC_QBASE(multiqc_qbase_ch, "${stage_label}_subset_qc_quality_base_stats")
        qseqs_out_ch = CONCATENATE_MULTIQC_QSEQS(multiqc_qseqs_ch, "${stage_label}_subset_qc_quality_sequence_stats")
        lengths_out_ch = CONCATENATE_MULTIQC_LENGTHS(multiqc_lengths_ch, "${stage_label}_subset_qc_length_stats")
        // 6. Combine outputs into a single output channel
        out_ch = basic_out_ch.combine(adapt_out_ch)
            .combine(qbase_out_ch).combine(qseqs_out_ch)
            .combine(lengths_out_ch)
            .map({file1, file2, file3, file4, file5 -> tuple(file1, file2, file3, file4, file5)})
    emit:
        qc = out_ch
}
