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
        samplesheet
        n_reads_trunc
        fastqc_cpus
        fastqc_mem
        stage_label
    main:
        concat_ch = samplesheet.map { sample, read1, read2 ->
            tuple(sample, [read1, read2])
        }
        if ( n_reads_trunc > 0 ) {
            out_ch = TRUNCATE_CONCAT(concat_ch, n_reads_trunc)
        } else {
            out_ch = concat_ch
        }
        qc_ch = QC(out_ch, fastqc_cpus, fastqc_mem, stage_label)
    emit:
        reads = out_ch
        qc = qc_ch.qc
}
