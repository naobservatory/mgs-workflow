// Download & concatenate ribosomal references
process JOIN_RIBO_REF {
    label "BBTools"
    label "single"
    input:
        val(ssu_url)
        val(lsu_url)
    output:
        path("ribo-ref-concat.fasta.gz"), emit: ribo_ref
    shell:
        '''
        # Download references
        wget !{ssu_url} -O ssu_ref.fasta.gz
        wget !{lsu_url} -O lsu_ref.fasta.gz
        in="ssu_ref.fasta.gz lsu_ref.fasta.gz"
        cat ${in} > ribo-ref-concat.fasta.gz
        '''
}
