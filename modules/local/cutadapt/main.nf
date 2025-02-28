process CUTADAPT {
    label "cutadapt"
    label "small"
    input:
        tuple val(sample), path(reads_interleaved)
        path(adapters)
    output:
        tuple val(sample), path("${sample}_cutadapt.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_cutadapt_log.txt"), emit: log
        tuple val(sample), path("${sample}_cutadapt_in.fastq.gz"), emit: input
    shell:
        /* Explanation of cutadapt parameters:
        -b (-B for R2 read) to trim any adapter in the adapters file from either end of either read
        -j number of cpu cores
        -m to drop a read pair where either element of the pair is <20 bp after trimming
        -e 0.125 To allow up to a 12.5% error rate in the matching region between an adapter
            and the read
        --action=trim to trim adapters and up/downstream sequence
        */
        '''
        output="!{sample}_cutadapt.fastq.gz"
        log="!{sample}_cutadapt_log.txt"
        par="-b file:!{adapters} -B file:!{adapters} -j !{task.cpus} -m 20 -e 0.125 --action=trim --interleaved"
        zcat !{reads_interleaved} | cutadapt ${par} - 2> ${log} | gzip -c > ${output}
        ln -s !{reads_interleaved} !{sample}_cutadapt_in.fastq.gz
        '''
}

process CUTADAPT_MASK {
    label "cutadapt"
    label "large"
    input:
        path(seq_db)
        path(adapters)
        path(label)
    output:
        path("${label}_cutadapt_masked.fasta.gz"), emit: masked
        path("${label}_cutadapt_log.txt"), emit: log
    shell:
        ''' 
        out=!{label}_cutadapt_masked.fasta.gz
        log=!{label}_cutadapt_log.txt
        cutadapt --action=mask -b file:!{adapters} -j !{task.cpus} -e 0.2 -o ${out} !{seq_db} > ${log}
        '''
}
