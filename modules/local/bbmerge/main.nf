// Merge read pairs into a single sequence
process BBMERGE {
    label "BBTools"
    label "small"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_joined.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_bbmerge_{stats,log}.txt"), emit: log
    shell:
        '''
        # Prepare input/output for bbmerge
        in1=!{reads[0]}
        in2=!{reads[1]}
        om=!{sample}_joined.fastq.gz
        stats=!{sample}_bbmerge_stats.txt
        log=!{sample}_bbmerge_log.txt
        io="in=${in1} in2=${in2} out=${om} ihist=${stats}"
        # Execute bbmerge
        bbmerge-auto.sh ${io} mix &> ${log}
        # Check for empty output files due to errors
        if [[ ! -s ${om} ]]; then
            >&2 echo "Error: Empty output file."
            exit 1
        fi
        '''
}
