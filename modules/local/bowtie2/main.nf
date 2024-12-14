// Run Bowtie2 and return mapped and unmapped reads
process BOWTIE2_PAIRED {
    label "Bowtie2"
    label "large"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
        val(par_string)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bowtie2_mapped.sam"), emit: sam
        tuple val(sample), path("${sample}_${suffix}_bowtie2_conc_{1,2}.fastq.gz"), emit: reads_conc
        tuple val(sample), path("${sample}_${suffix}_bowtie2_unconc_{1,2}.fastq.gz"), emit: reads_unconc
    shell:
        '''
        in1=!{reads[0]}
        in2=!{reads[1]}
        idx="!{index_dir}/bt2_index"
        sam="!{sample}_!{suffix}_bowtie2_mapped.sam"
        alc="!{sample}_!{suffix}_bowtie2_conc_%.fastq.gz"
        unc="!{sample}_!{suffix}_bowtie2_unconc_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc} --un-conc-gz ${unc}"
        par="--threads !{task.cpus} --local --very-sensitive-local !{par_string}"
        bowtie2 ${par} ${io}
        '''
}

process BOWTIE2_SINGLE {
    label "Bowtie2"
    label "large"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
        val(par_string)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bowtie2_mapped.sam"), emit: sam
        tuple val(sample), path("${sample}_${suffix}_bowtie2_aligned.fastq.gz"), emit: reads_conc
        tuple val(sample), path("${sample}_${suffix}_bowtie2_unaligned.fastq.gz"), emit: reads_unconc
    shell:
        '''
        in=!{reads}
        idx="!{index_dir}/bt2_index"
        sam="!{sample}_!{suffix}_bowtie2_mapped.sam"
        alc="!{sample}_!{suffix}_bowtie2_aligned.fastq.gz"
        unc="!{sample}_!{suffix}_bowtie2_unaligned.fastq.gz"
        io="-U ${in} -x ${idx} -S ${sam} --al-gz ${alc} --un-gz ${unc}"
        par="--threads !{task.cpus} --local --very-sensitive-local !{par_string}"
        bowtie2 ${par} ${io}
        '''
}

// Generate a Bowtie2 index from an input file
process BOWTIE2_INDEX {
    label "Bowtie2"
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
        bowtie2-build -f --threads !{task.cpus} !{reference_fasta} ${odir}/bt2_index
        #tar -czf bt2-human-index.tar.gz bt2_human_index
        '''
}
