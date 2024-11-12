// Masking low-complexity and repetitive sequences in a reference database
process BBMASK{
    label "large"
    label "BBTools"
    input:
        path(seq_db)
        val(label)
    output:
        path("${label}_bbmask_masked.fasta.gz"), emit: masked
        path("${label}_bbmask.log.txt"), emit: log
    shell:
        '''
        # Define input/output
        in=!{seq_db}
        out=!{label}_bbmask_masked.fasta.gz
        log=!{label}_bbmask.log.txt
        par="maskrepeats=t minkr=5 maxkr=10 masklowentropy=t t=!{task.cpus} -Xmx30g"
        # Execute
        bbmask.sh in=${in} out=${out} ${par} > ${log}
        '''
}
