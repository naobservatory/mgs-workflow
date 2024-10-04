/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC} from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)

if ( params.read_type == "single_end" ) {
    include { CONCAT_GZIPPED_SINGLE as CONCAT_GZIPPED } from "../../../modules/local/concatGzipped" addParams(r1_suffix: params.r1_suffix, r2_suffix: params.r2_suffix)
} else if ( params.read_type == "paired_end" ) {
    include { TRUNCATE_CONCAT } from "../../../modules/local/truncateConcat"
    include { CONCAT_GZIPPED_PAIRED as CONCAT_GZIPPED } from "../../../modules/local/concatGzipped" addParams(r1_suffix: params.r1_suffix, r2_suffix: params.r2_suffix)
}

/***********
| WORKFLOW |
***********/

workflow RAW {
    take:
        libraries_ch
        raw_dir_path
        n_reads_trunc
    main:
        concat_ch = CONCAT_GZIPPED(raw_dir_path, libraries_ch)
        if ( n_reads_trunc > 0 ) {
            out_ch = TRUNCATE_CONCAT(concat_ch, n_reads_trunc)
            qc_ch = QC(trunc_ch)
        } else {
            out_ch = concat_ch
        }
        qc_ch = QC(out_ch.reads, params.stage_label)
    emit:
        reads = out_ch.reads
        qc = qc_ch.qc
}
