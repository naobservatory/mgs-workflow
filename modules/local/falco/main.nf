// https://docs.seqera.io/multiqc/modules/fastqc
process FALCO {
    label "FALCO"
    cpus "${params.cpus}"
    memory "${params.mem}"
    input:
       tuple val(sample), path(reads)
    output:
        path("*.zip"), emit: data
    shell:
        '''
        in1=!{reads[0]}
        in2=!{reads[1]}
        for i in ${in1} ${in2}; do
            falco -t !{params.cpus} ${i} -o ${i}_fastqc
            zip -r ${i}_fastqc.zip ${i}_fastqc
        done
        '''
}
