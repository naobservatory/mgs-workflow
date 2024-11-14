// Basecall Nanopore pod5 files with Dorado
process BASECALL_POD_5 {
    label "dorado"
    label "basecall"
    accelerator 1
    memory '8 GB'

    input:
        path(pod_5_dir)
        val kit

    output:
        path("calls.bam"), emit: bam
        path("sequencing_summary.txt"), emit: summary

    shell:
        '''
        # Extract batch number using bash
        # batch_num=$(basename !{pod_5_dir} | grep -o '[0-9]\\+') # Disabled in the absence of batching

        # Dorado basecalling
        dorado basecaller sup !{pod_5_dir} --kit-name !{kit} > calls.bam

        dorado summary calls.bam > sequencing_summary.txt
        '''
}

// Demultiplex basecalled BAM files
process DEMUX_POD_5 {
    label "dorado"
    label "demux"
    accelerator 1
    memory '8 GB'

    input:
        path(calls_bam)
        val kit
    output:
        path('demultiplexed/*'), emit: demux_bam

    shell:
        '''
        dorado demux --output-dir demultiplexed/ --kit-name !{kit} !{calls_bam}
        '''
}