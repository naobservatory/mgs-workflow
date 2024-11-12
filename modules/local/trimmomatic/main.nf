process TRIMMOMATIC {
    label "trimmomatic"
    label "large"
    input:
        // reads is a list of two files: forward/reverse reads
        tuple val(sample), path(reads)
        path(adapters)
        val(encoding}
    output:
        tuple val(sample), path("${sample}_trimmomatic_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_trimmomatic_unpaired_{1,2}.fastq.gz"), emit: unpaired
        tuple val(sample), path("${sample}_trimmomatic_{trimlog,summary}.txt"), emit: log
    shell:
        /* Explanation of trimmomatic filters:
        ILLUMINACLIP: cut Illumina adapters
            * max 2 mismatches in a match
            * min match quality of 20 between reads for paired end read alignment
            * min match quality of 8 between adapter and read
            * adapters must be >= 5 bases
        LEADING/TRAILING:10 remove leading/trailing bases with quality < 10
        SLIDINGWINDOW:4:15 cut the read after the average quality of a 4-base
            window falls below 15
        MINLEN:20 drop reads below 20 bases (does not drop paired reads)
        */
        '''
        op1=!{sample}_trimmomatic_1.fastq.gz
        op2=!{sample}_trimmomatic_2.fastq.gz
        of1=!{sample}_trimmomatic_unpaired_1.fastq.gz
        of2=!{sample}_trimmomatic_unpaired_2.fastq.gz
        log=!{sample}_trimmomatic_trimlog.txt
        sum=!{sample}_trimmomatic_summary.txt
        io="-trimlog ${log} -summary ${sum} !{reads[0]} !{reads[1]} ${op1} ${of1} ${op2} ${of2}"
        par="ILLUMINACLIP:!{adapters}:2:20:8:5 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:20"
        cmd="trimmomatic PE -!{encoding} -threads !{task.cpus} ${io} ${par}"
        echo ${cmd}
        ${cmd}
        '''
}
