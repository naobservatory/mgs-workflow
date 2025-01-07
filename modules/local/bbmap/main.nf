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
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbmap.sh ${io} ${par}
        '''
}

// Run BBMap on streamed interleaved input and return mapped and unmapped reads
// NB: This is configured to work the same way as BOWTIE2_STREAMED, returning a SAM file that is then partitioned to find mapped and unmapped reads
process BBMAP_STREAMED {
    label "bbtools_samtools"
    label "small"
    input:
        tuple val(sample), path(reads_interleaved)
        path(index_dir)
        val(suffix)
        val(remove_sq)
        val(debug)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbmap_mapped.sam.gz"), emit: sam
        tuple val(sample), path("${sample}_${suffix}_bbmap_mapped.fastq.gz"), emit: reads_mapped
        tuple val(sample), path("${sample}_${suffix}_bbmap_unmapped.fastq.gz"), emit: reads_unmapped
        tuple val(sample), path("${sample}_${suffix}_bbmap_in.fastq.gz"), emit: input
        tuple val(sample), path("${sample}_${suffix}_bbmap.stats.txt"), emit: stats
    shell:
        '''
        set -euo pipefail
        # Prepare inputs
        idx="!{index_dir}"
        sam="!{sample}_!{suffix}_bbmap_mapped.sam.gz"
        al="!{sample}_!{suffix}_bbmap_mapped.fastq.gz"
        un="!{sample}_!{suffix}_bbmap_unmapped.fastq.gz"
        stats=!{sample}_!{suffix}_bbmap.stats.txt
        io="in=stdin.fastq path=${idx} statsfile=${stats} out=stdout.sam"
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch interleaved minhits=2 t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Run pipeline
        zcat !{reads_interleaved} \\
            | bbmap.sh ${io} ${par} \\
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
        # Link input files for testing
        in2="!{sample}_!{suffix}_bbmap_in.fastq.gz"
        ln -s !{reads_interleaved} ${in2}
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
        bbmap.sh ref=reference.fasta.gz t=!{task.cpus} -Xmx!{task.memory.toGiga()}g
        #tar -czf human-ref-index.tar.gz human_ref_index
        '''
}
