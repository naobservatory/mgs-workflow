// Run BBMap and return mapped and unmapped reads
process BBMAP {
    label "max"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbmap_mapped_{1,2}.fastq.gz"), emit: reads_mapped
        tuple val(sample), path("${sample}_${suffix}_bbmap_unmapped_{1,2}.fastq.gz"), emit: reads_unmapped
        tuple val(sample), path("${sample}_${suffix}_bbmap.stats.txt"), emit: stats
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        ou1=!{sample}_!{suffix}_bbmap_unmapped_1.fastq.gz
        ou2=!{sample}_!{suffix}_bbmap_unmapped_2.fastq.gz
        om1=!{sample}_!{suffix}_bbmap_mapped_1.fastq.gz
        om2=!{sample}_!{suffix}_bbmap_mapped_2.fastq.gz
        stats=!{sample}_!{suffix}_bbmap.stats.txt
        io="in=${in1} in2=${in2} outu=${ou1} outu2=${ou2} outm=${om1} outm2=${om2} statsfile=${stats} path=!{index_dir}"
        # Define parameters
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx60g"
        # Execute
        bbmap.sh ${io} ${par}
        '''
}

// Generate a BBMap index from an input file
process BBMAP_INDEX {
    label "BBTools"
    label "max"
    input:
        path(reference_fasta)
        val(outdir)
    output:
        path("${outdir}")
    shell:
        '''
        odir="!{outdir}"
        mkdir ${odir}
        cp !{reference_fasta} ${odir}/reference.fasta.gz
        cd ${odir}
        bbmap.sh ref=reference.fasta.gz t=!{task.cpus} -Xmx10g
        #tar -czf human-ref-index.tar.gz human_ref_index
        '''
}
