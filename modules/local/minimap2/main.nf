// Create a minimap2 index
process MINIMAP2_INDEX {
    label "large"
    label "minimap2"
    input:
        path(reference_fasta)
        val(outdir)
    output:
        path("${outdir}")
    shell:
        '''
        odir="!{outdir}"
        mkdir ${odir}
        preset="lr:hq"
        minimap2 -x ${preset} -d ${odir}/mm2_index.mmi !{reference_fasta}
        '''
}
