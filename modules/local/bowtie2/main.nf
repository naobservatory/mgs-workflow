// Run Bowtie2 and return mapped and unmapped reads
process BOWTIE2 {
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

// Run Bowtie2 on streamed interleaved input and return mapped and unmapped reads
// NB: This handles non-concordant alignments correctly for this use case (including them with the aligned rather than unaligned reads), so we can skip some downstream processing steps
process BOWTIE2_STREAMED {
    label "bowtie2_samtools"
    label "small"
    input:
        tuple val(sample), path(reads_interleaved)
        path(index_dir)
        val(par_string)
        val(suffix)
        val(remove_sq)
        val(debug)
    output:
        tuple val(sample), path("${sample}_${suffix}_bowtie2_mapped.sam.gz"), emit: sam
        tuple val(sample), path("${sample}_${suffix}_bowtie2_mapped.fastq.gz"), emit: reads_mapped
        tuple val(sample), path("${sample}_${suffix}_bowtie2_unmapped.fastq.gz"), emit: reads_unmapped
        tuple val(sample), path("${sample}_${suffix}_bowtie2_in.fastq.gz"), emit: input
    shell:
        '''
        set -euo pipefail
        # Prepare inputs
        idx="!{index_dir}/bt2_index"
        sam="!{sample}_!{suffix}_bowtie2_mapped.sam.gz"
        al="!{sample}_!{suffix}_bowtie2_mapped.fastq.gz"
        un="!{sample}_!{suffix}_bowtie2_unmapped.fastq.gz"
        io="-x ${idx} --interleaved -"
        par="--threads !{task.cpus} --local --very-sensitive-local !{par_string}"
        # Run pipeline
        zcat !{reads_interleaved} \\
            | bowtie2 ${par} ${io} \\
            | tee \\
                !{ debug ? ">(gzip -c > test_all.sam.gz)" : "" } \\
                >(samtools view -u -f 12 - \\
                    !{ debug ? "| tee >(samtools view -h - | gzip -c > test_unmapped.sam.gz)" : "" } \\
                    | samtools fastq -1 /dev/stdout -2 /dev/stdout \\
                        -0 /dev/stdout -s /dev/stdout - \\
                    | gzip -c > ${un}) \\
                >(samtools view -u -G 12 - \\
                    !{ debug ? "| tee >(samtools view -h - | gzip -c > test_mapped.sam.gz)" : "" } \\
                    | samtools fastq -1 /dev/stdout -2 /dev/stdout \\
                        -0 /dev/stdout -s /dev/stdout - \\
                    | gzip -c > ${al}) \\
            | samtools view -h -G 12 - \\
            !{ remove_sq ? "| grep -v '^@SQ'" : "" } | gzip -c > ${sam}
        # Move input files for testing
        in2="!{sample}_!{suffix}_bowtie2_in.fastq.gz"
        mv !{reads_interleaved} ${in2}
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
