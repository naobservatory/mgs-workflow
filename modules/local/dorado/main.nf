// Basecall Nanopore pod5 files with Dorado
process BASECALL_POD_5 {
    label "dorado"
    label "basecall"
    accelerator 1
    memory '8 GB'

    input:
        path(pod5_dir)

    output:
        path("calls.bam"), emit: bam

    shell:
        '''
        # Dorado basecalling
        dorado basecaller sup !{pod5_dir} > calls.bam
        '''
}