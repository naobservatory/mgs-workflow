process FASTQC_LABELED {
    label "FASTQC"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("*.html"), emit: html
        tuple val(sample), path("*.zip"), emit: zip
    shell:
        '''
        # Note that values over 10000 for the memory flag are invalid and cause fastqc to fail; hence the min() 
        fastqc -t !{task.cpus} --memory !{Math.min(10000, task.memory.toMega())} !{reads}
        '''
}
