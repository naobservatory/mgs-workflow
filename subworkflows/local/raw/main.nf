/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc"
include { TRUNCATE_CONCAT } from "../../../modules/local/truncateConcat"

/***********
| WORKFLOW |
***********/

workflow RAW {
    take:
        samplesheet_ch
        n_reads_trunc
        fastqc_cpus
        fastqc_mem
        stage_label
        single_end
    main:
        if ( n_reads_trunc > 0 ) {
            out_ch = TRUNCATE_CONCAT(samplesheet_ch, n_reads_trunc)
        } else {
            out_ch = samplesheet_ch
        }
        qc_ch = QC(out_ch, fastqc_cpus, fastqc_mem, stage_label, single_end)
    emit:
        reads = out_ch
        qc = qc_ch.qc
}
