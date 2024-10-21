/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)

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
        samplesheet
        n_reads_trunc
    main:
        if (params.read_type == "single_end") {
            concat_ch = samplesheet.map { sample, read ->
                tuple(sample, [read])
            }
        } else if (params.read_type == "paired_end") {
            concat_ch = samplesheet.map { sample, read1, read2 ->
                tuple(sample, [read1, read2])
            }
        }
        if ( n_reads_trunc > 0 ) {
                out_ch = TRUNCATE_CONCAT(concat_ch, n_reads_trunc)
            } else {
                out_ch = concat_ch
            }
            qc_ch = QC(out_ch, params.stage_label)
    emit:
        reads = out_ch
        qc = qc_ch.qc
}
