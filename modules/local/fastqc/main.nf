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
        fastqc -t !{task.cpus} --memory !{task.memory.toMega()} !{reads}
        '''
}
