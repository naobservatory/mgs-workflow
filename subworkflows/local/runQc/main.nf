/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC as PRE_ADAPTER_TRIM_QC } from "../qc"
include { QC as POST_ADAPTER_TRIM_QC } from "../qc"
include { MERGE_TSVS as MERGE_MULTIQC_BASIC } from "../../../modules/local/mergeTsvs"
include { MERGE_TSVS as MERGE_MULTIQC_ADAPT } from "../../../modules/local/mergeTsvs"
include { MERGE_TSVS as MERGE_MULTIQC_QBASE } from "../../../modules/local/mergeTsvs"
include { MERGE_TSVS as MERGE_MULTIQC_QSEQS } from "../../../modules/local/mergeTsvs"

/***********
| WORKFLOW |
***********/

workflow RUN_QC {
    take:
      subset_reads
      trimmed_subset_reads
      fastqc_cpus
      fastqc_mem
      single_end
    main:
      // 1. Run FASTQC before and after adapter trimming
      pre_qc_ch = PRE_ADAPTER_TRIM_QC(subset_reads, fastqc_cpus, fastqc_mem, "pre_trim", single_end)
      post_qc_ch = POST_ADAPTER_TRIM_QC(trimmed_subset_reads, fastqc_cpus, fastqc_mem, "post_trim", single_end)
      // 2. Combine outputs
      qc_ch = pre_qc_ch.concat(post_qc_ch)
      // 3. Collate MultiQC outputs
      multiqc_basic_ch = qc_ch.map{ it[0] }.collect().ifEmpty([])
      multiqc_adapt_ch = qc_ch.map{ it[1] }.collect().ifEmpty([])
      multiqc_qbase_ch = qc_ch.map{ it[2] }.collect().ifEmpty([])
      multiqc_qseqs_ch = multiqc_ch.map{ it[3] }.collect().ifEmpty([])
      // 4. Merge MultiQC outputs
      basic_out_ch = MERGE_MULTIQC_BASIC(multiqc_basic_ch, "qc_basic_stats")
      adapt_out_ch = MERGE_MULTIQC_ADAPT(multiqc_adapt_ch, "qc_adapter_stats")
      qbase_out_ch = MERGE_MULTIQC_QBASE(multiqc_qbase_ch, "qc_quality_base_stats")
      qseqs_out_ch = MERGE_MULTIQC_QSEQS(multiqc_qseqs_ch, "qc_quality_sequence_stats")
    emit:
      qc_basic = basic_out_ch
      qc_adapt = adapt_out_ch
      qc_qbase = qbase_out_ch
      qc_qseqs = qseqs_out_ch
}
