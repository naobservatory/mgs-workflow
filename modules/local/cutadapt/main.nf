process CUTADAPT {
    label "cutadapt"
    label "large"
    input:
        tuple val(sample), path(reads)
        path(adapters)
    output:
        tuple val(sample), path("${sample}_cutadapt_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_cutadapt_log.txt"), emit: log
    shell:
        '''
        par="-b file:!{adapters} -B file:!{adapters} -j !{task.cpus} -m 20 -e 0.33 --action=trim"
        out="-o !{sample}_cutadapt_1.fastq.gz -p !{sample}_cutadapt_2.fastq.gz"
        log="!{sample}_cutadapt_log.txt"
        cutadapt ${par} ${out} !{reads[0]} !{reads[1]} > ${log}
        '''
}

process CUTADAPT_MASK {
    label "cutadapt"
    label "large"
    input:
        path(seq_db)
        path(adapters)
    output:
        path("${params.label}_cutadapt_masked.fasta.gz"), emit: masked
        path("${params.label}_cutadapt_log.txt"), emit: log
    shell:
        ''' 
        out=!{params.label}_cutadapt_masked.fasta.gz
        log=!{params.label}_cutadapt_log.txt
        cutadapt --action=mask -b file:!{adapters} -j !{task.cpus} -e 0.2 -o ${out} !{seq_db} > ${log}
        '''
}
