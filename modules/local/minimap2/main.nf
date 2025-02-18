// Create a minimap2 index
process MINIMAP2_INDEX {
    label "large"
    label "minimap2"
    input:
        path(reference_fasta)
        val(outdir)
    output:
        path("${outdir}"), emit: output
        path("input_${reference_fasta}"), emit: input

    shell:
        '''
        odir="!{outdir}"
        mkdir ${odir}
        preset="lr:hq"
        minimap2 -x ${preset} -d ${odir}/mm2_index.mmi !{reference_fasta}

        # Link input to output for testing
        ln -s !{reference_fasta} input_!{reference_fasta}
        '''
}
