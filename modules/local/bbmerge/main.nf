// Merge read pairs into a single sequence
process BBMERGE {
    label "BBTools"
    label "small"
    input:
        tuple val(sample), path(reads_interleaved)
    output:
        tuple val(sample), path("${sample}_bbmerge_{merged,unmerged}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_bbmerge_{stats,log}.txt"), emit: log
        tuple val(sample), path("${sample}_bbmerge_in.fastq.gz"), emit: input
    shell:
        '''
        set -euo pipefail
        # Prepare inputs and outputs
        ou=!{sample}_bbmerge_unmerged.fastq.gz
        om=!{sample}_bbmerge_merged.fastq.gz
        stats=!{sample}_bbmerge_stats.txt
        log=!{sample}_bbmerge_log.txt
        io="in=stdin.fastq out=${om} outu=${ou} ihist=${stats}"
        par="join interleaved t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        zcat !{reads_interleaved} \\
            | bbmerge.sh ${io} ${par} &> ${log}
        # Check for empty output files due to errors
        if [[ ! -s ${ou} ]] && [[ ! -s ${om} ]]; then
            >&2 echo "Error: Empty output files."
            exit 1
        fi
        # Link input reads for testing
        in2=!{sample}_bbmerge_in.fastq.gz
        ln -s !{reads_interleaved} ${in2}
        '''
}
