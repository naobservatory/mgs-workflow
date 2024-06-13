/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
include { BBDUK } from "../../../modules/local/bbduk"

/***********
| WORKFLOW |
***********/

workflow RIBODEPLETION {
    take:
        reads_ch
        ribo_path
    main:
        bbduk_ch = BBDUK(reads_ch, ribo_path, params.min_kmer_fraction, params.k)
        qc_ch = QC(bbduk_ch.reads, params.stage_label)
    emit:
        reads = bbduk_ch.reads
        ribo = bbduk_ch.fail
        qc = qc_ch.out.qc
}
