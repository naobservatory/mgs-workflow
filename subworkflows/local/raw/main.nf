/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
include { TRUNCATE_CONCAT } from "../../../modules/local/truncateConcat"

/***********
| WORKFLOW |
***********/

workflow RAW {
    take:
        libraries_ch
        raw_dir_path
        n_reads_trunc
    main:
        concat_ch = libraries_ch.map { sample, read1, read2 ->
            tuple(sample, [read1, read2])
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
