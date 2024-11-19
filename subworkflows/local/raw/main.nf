/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc"
if (params.read_type == "single_end") {
    include { TRUNCATE_CONCAT_SINGLE as TRUNCATE_CONCAT } from "../../../modules/local/truncateConcat"
} else if (params.read_type == "paired_end") {
    include { TRUNCATE_CONCAT_PAIRED as TRUNCATE_CONCAT } from "../../../modules/local/truncateConcat"
}


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
    main:
        if ( n_reads_trunc > 0 ) {
            out_ch = TRUNCATE_CONCAT(samplesheet_ch, n_reads_trunc)
        } else {
            out_ch = samplesheet_ch
        }
        qc_ch = QC(out_ch, fastqc_cpus, fastqc_mem, stage_label)
    emit:
        reads = out_ch
        qc = qc_ch.qc
}
