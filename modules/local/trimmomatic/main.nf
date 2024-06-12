process TRIMMOMATIC {
    label "trimmomatic"
    label "large"
    input:
        tuple val(sample), path(reads)
        path(adapters)
    output:
        tuple val(sample), path("${sample}_trimmomatic_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_trimmomatic_unpaired_{1,2}.fastq.gz"), emit: unpaired
        tuple val(sample), path("${sample}_trimmomatic_{trimlog,summary}.txt"), emit: log
    shell:
        '''
        op1=!{sample}_trimmomatic_1.fastq.gz
        op2=!{sample}_trimmomatic_2.fastq.gz
        of1=!{sample}_trimmomatic_unpaired_1.fastq.gz
        of2=!{sample}_trimmomatic_unpaired_2.fastq.gz
        log=!{sample}_trimmomatic_trimlog.txt
        sum=!{sample}_trimmomatic_summary.txt
        io="-trimlog ${log} -summary ${sum} !{reads[0]} !{reads[1]} ${op1} ${of1} ${op2} ${of2}"
        par="ILLUMINACLIP:!{adapters}:2:20:8:5 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:20"
        cmd="trimmomatic PE -threads !{task.cpus} ${io} ${par}"
        echo ${cmd}
        ${cmd}
        '''
}
