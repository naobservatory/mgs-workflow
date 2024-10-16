// Merge read pairs into a single sequence
process BBMERGE {
    label "BBTools"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_bbmerge_{merged,unmerged_1,unmerged_2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_bbmerge_{stats,log}.txt"), emit: log
    shell:
        '''
        # Prepare input/output for bbmerge
        in1=!{reads[0]}
        in2=!{reads[1]}
        ou1=!{sample}_bbmerge_unmerged_1.fastq.gz
        ou2=!{sample}_bbmerge_unmerged_2.fastq.gz
        om=!{sample}_bbmerge_merged.fastq.gz
        stats=!{sample}_bbmerge_stats.txt
        log=!{sample}_bbmerge_log.txt
        io="in=${in1} in2=${in2} out=${om} outu=${ou1} outu2=${ou2} ihist=${stats}"
        # Execute bbmerge
        bbmerge.sh ${io} &> ${log}
        # Check for empty output files due to errors
        if [[ ! -s ${ou1} ]] && [[ ! -s ${om} ]]; then
            >&2 echo "Error: Empty output files."
            exit 1
        fi
        '''
}
