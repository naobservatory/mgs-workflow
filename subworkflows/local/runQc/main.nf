/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC as PRE_ADAPTER_TRIM_QC } from "../qc"
include { QC as POST_ADAPTER_TRIM_QC } from "../qc"
include { PROCESS_OUTPUT } from "../processOutput"

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
      PROCESS_OUTPUT(qc_ch)
    emit:
      qc_basic = PROCESS_OUTPUT.out.basic
      qc_adapt = PROCESS_OUTPUT.out.adapt
      qc_qbase = PROCESS_OUTPUT.out.qbase
      qc_qseqs = PROCESS_OUTPUT.out.qseqs
}
