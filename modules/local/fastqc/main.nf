process FASTQC {
    label "FASTQC"
    cpus "${cpus}"
    memory "${mem}"
    input:
        tuple val(sample), path(reads)
        val(cpus)
        val(mem)
    output:
        path("*.html"), emit: html
        path("*.zip"), emit: zip
    shell:
        '''
        fastqc -t !{cpus} --memory !{task.memory.toMega()} !{reads}
        '''
}

process FASTQC_LABELED {
    label "FASTQC"
    cpus "${cpus}"
    memory "${mem}"
    input:
        tuple val(sample), path(reads)
        val(cpus)
        val(mem)
    output:
        tuple val(sample), path("*.html"), emit: html
        tuple val(sample), path("*.zip"), emit: zip
    shell:
        '''
        fastqc -t !{task.cpus} --memory !{task.memory.toMega()} !{reads}
        '''
}
