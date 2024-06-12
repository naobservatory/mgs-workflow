/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus)
include { CONCAT_GZIPPED } from "../modules/local/concatGzipped"
include { TRUNCATE_CONCAT } from "../modules/local/truncateConcat"

/***********
| WORKFLOW |
***********/

workflow RAW {
    take:
        libraries_ch
        raw_dir_path
        truncate_reads
        n_reads_trunc
    main:
        concat_ch = CONCAT_GZIPPED(raw_dir_path, libraries_ch)
        if ( truncate_reads ) {
            out_ch = TRUNCATE_CONCAT(concat_ch, n_reads_trunc)
            qc_ch = QC(trunc_ch)
        } else {
            out_ch = concat_ch
        }
        qc_ch = QC(out_ch)
    emit:
        data = out_ch
        qc = qc_ch
}
