// Run Bowtie2 and return mapped and unmapped reads
process BOWTIE2 {
    label "Bowtie2"
    label "large"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
        val(index_prefix)
        val(par_string)
    output:
        tuple val(sample), path("${sample}_${params.suffix}_bowtie2_mapped.sam"), emit: sam
        tuple val(sample), path("${sample}_${params.suffix}_bowtie2_conc_{1,2}.fastq.gz"), emit: reads_conc
        tuple val(sample), path("${sample}_${params.suffix}_bowtie2_unconc_{1,2}.fastq.gz"), emit: reads_unconc
    shell:
        '''
        in1=!{reads[0]}
        in2=!{reads[1]}
        idx="!{index_dir}/!{index_prefix}"
        sam="!{sample}_!{params.suffix}_bowtie2_mapped.sam"
        alc="!{sample}_!{params.suffix}_bowtie2_conc_%.fastq.gz"
        unc="!{sample}_!{params.suffix}_bowtie2_unconc_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc} --un-conc-gz ${unc}"
        par="--threads !{task.cpus} --local --very-sensitive-local !{par_string}"
        bowtie2 ${par} ${io}
        '''
}
