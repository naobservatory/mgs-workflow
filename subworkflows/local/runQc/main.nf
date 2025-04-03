/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC as PRE_ADAPTER_TRIM_QC } from "../../../subworkflows/local/qc"
include { QC as POST_ADAPTER_TRIM_QC } from "../../../subworkflows/local/qc"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_BASIC } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_ADAPT } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QBASE } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_QSEQS } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_MULTIQC_LENGTHS } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow RUN_QC {
    take:
      subset_reads
      trimmed_subset_reads
      single_end
    main:
      // 1. Run FASTQC before and after adapter trimming
      pre_qc_ch = PRE_ADAPTER_TRIM_QC(subset_reads, "raw", single_end)
      post_qc_ch = POST_ADAPTER_TRIM_QC(trimmed_subset_reads, "cleaned", single_end)
      // 2. Combine outputs
      qc_ch = pre_qc_ch.qc.concat(post_qc_ch.qc)
      // 3. Collate MultiQC outputs
      multiqc_basic_ch = qc_ch.map{ it[0] }.collect().ifEmpty([])
      multiqc_adapt_ch = qc_ch.map{ it[1] }.collect().ifEmpty([])
      multiqc_qbase_ch = qc_ch.map{ it[2] }.collect().ifEmpty([])
      multiqc_qseqs_ch = qc_ch.map{ it[3] }.collect().ifEmpty([])
      multiqc_lengths_ch = qc_ch.map{ it[4] }.collect().ifEmpty([])
      // 4. Merge MultiQC outputs
      basic_out_ch = CONCATENATE_MULTIQC_BASIC(multiqc_basic_ch, "subset_qc_basic_stats")
      adapt_out_ch = CONCATENATE_MULTIQC_ADAPT(multiqc_adapt_ch, "subset_qc_adapter_stats")
      qbase_out_ch = CONCATENATE_MULTIQC_QBASE(multiqc_qbase_ch, "subset_qc_quality_base_stats")
      qseqs_out_ch = CONCATENATE_MULTIQC_QSEQS(multiqc_qseqs_ch, "subset_qc_quality_sequence_stats")
      lengths_out_ch = CONCATENATE_MULTIQC_LENGTHS(multiqc_lengths_ch, "subset_qc_length_stats")
    emit:
      qc_basic = basic_out_ch.output
      qc_adapt = adapt_out_ch.output
      qc_qbase = qbase_out_ch.output
      qc_qseqs = qseqs_out_ch.output
      qc_lengths = lengths_out_ch.output
}
