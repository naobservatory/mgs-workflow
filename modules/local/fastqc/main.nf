process FASTQC {
    label "FASTQC"
    cpus "${params.cpus}"
    memory "${params.mem}"
    input:
       tuple val(sample), path(reads)
    output:
        path("*.html"), emit: html
        path("*.zip"), emit: zip
    shell:
        '''
        fastqc -t !{params.cpus} !{reads}
        '''
}

process FASTQC_LABELED {
    label "FASTQC"
    cpus "${params.cpus}"
    memory "${params.mem}"
    input:
       tuple val(sample), path(reads)
    output:
        tuple val(sample), path("*.html"), emit: html
        tuple val(sample), path("*.zip"), emit: zip
    shell:
        '''
        fastqc -t !{params.cpus} !{reads}
        '''
}
