process FASTQC {
    label "FASTQC"
    cpus "${params.cpus}"
    input:
        path(reads)
    output:
        path("*.html"), emit: html
        path("*.zip"), emit: zip
    shell:
        '''
        fastqc -t !{params.cpus} !{reads}
        '''
}
