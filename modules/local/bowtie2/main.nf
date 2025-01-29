// Run Bowtie2 on streamed interleaved input and return mapped and unmapped reads
// NB: This handles non-concordant alignments correctly for this use case (including them with the aligned rather than unaligned reads), so we can skip some downstream processing steps
process BOWTIE2 {
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
        # Outputs a SAM file for all reads, which is then partitioned based on alignment status
        #   - First branch (samtools view -u -f 12 -) filters SAM to read pairs for which neither mate mapped,
        #       then extracts and saves FASTQ
        #   - Second branch (samtools view -u -G 12) filters SAM to read pairs for which either mate mapped,
        #       then extracts and saves FASTQ
        #   - Third branch (samtools view -h -G 12) also filters SAM to read pairs for which either mate mapped,
        #       optionally removes SQ header lines, then saves SAM
        # Debug statements allow saving of additional SAM files at different steps in the pipeline.
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
        ln -s !{reads_interleaved} ${in2}
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
