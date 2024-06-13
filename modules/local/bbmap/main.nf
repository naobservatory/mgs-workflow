// Run BBMap and return mapped and unmapped reads
process BBMAP {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_mapped_{1,2}.fastq.gz"), emit: reads_mapped
        tuple val(sample), path("${sample}_bbmap_unmapped_{1,2}.fastq.gz"), emit: reads_unmapped
        tuple val(sample), path("${sample}_bbmap.stats.txt"), emit: stats
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        ou1=!{sample}_bbmap_unmapped_1.fastq.gz
        ou2=!{sample}_bbmap_unmapped_2.fastq.gz
        om1=!{sample}_bbmap_mapped_1.fastq.gz
        om2=!{sample}_bbmap_mapped_2.fastq.gz
        stats=!{sample}_bbmap.stats.txt
        io="in=${in1} in2=${in2} outu=${ou1} outu2=${ou2} outm=${om1} outm2=${om2} statsfile=${stats} path=!{index_dir}"
        # Define parameters
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx30g"
        # Execute
        bbmap.sh ${io} ${par}
        '''
}
