// Download & concatenate ribosomal references
process JOIN_RIBO_REF {
    label "BBTools"
    label "single"
    input:
        val(ssu_url)
        val(lsu_url)
    output:
        path("ribo-ref-concat.fasta.gz"), emit: ribo_ref
        tuple path("ssu_ref_input.fasta.gz"), path("lsu_ref_input.fasta.gz"), emit: input
    shell:
        '''
        # Download references
        wget !{ssu_url} -O ssu_ref.fasta.gz
        wget !{lsu_url} -O lsu_ref.fasta.gz
        # Add suffixes to reference headers
        for ref in ssu lsu; do
            zcat ${ref}_ref.fasta.gz | awk -v suffix=$ref '
            /^>/ {
                pos = index($0, " ")
                print (pos > 0) ? substr($0,1,pos-1) "::" toupper(suffix) substr($0,pos) : $0 "::" toupper(suffix)
                next
            }
            { print }' | gzip > ${ref}_ref_suffix.fasta.gz
        done

        # Update input files for concatenation
        in="ssu_ref_suffix.fasta.gz lsu_ref_suffix.fasta.gz"
        cat ${in} > ribo-ref-concat.fasta.gz

        # Return input files for testing
        ln -s ssu_ref.fasta.gz ssu_ref_input.fasta.gz
        ln -s lsu_ref.fasta.gz lsu_ref_input.fasta.gz
        '''
}
