// Deduplicate paired reads with Clumpify
// NB: Will NOT handle reverse-complement duplicates
process CLUMPIFY_PAIRED {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_dedup_{1,2}.fastq.gz"), emit: reads
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        op1=!{sample}_dedup_1.fastq.gz
        op2=!{sample}_dedup_2.fastq.gz
        io="in=${in1} in2=${in2} out=${op1} out2=${op2}"
        # Define parameters
        par="reorder dedupe containment addcount=t t=!{task.cpus} -Xmx30g"
        # Execute
        clumpify.sh ${io} ${par}
        '''
}

// Deduplicate single/merged reads with Clumpify
// NB: Should handle RC duplicates
process CLUMPIFY_SINGLE {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_dedup.fastq.gz"), emit: reads
    shell:
        '''
        # Define input/output
        in=!{reads}
        out=!{sample}_dedup.fastq.gz
        io="in=${in} out=${out}"
        # Define parameters
        par="reorder dedupe containment rcomp passes=6 addcount=t t=!{task.cpus} -Xmx30g"
        # Execute
        clumpify.sh ${io} ${par}
        '''
}
