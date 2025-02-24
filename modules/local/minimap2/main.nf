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

// Run minimap2 on a single input FASTQ file and partition reads based on alignment status
process MINIMAP2 {
    label "large"
    label "minimap2_samtools"
    input:
        tuple val(sample), path(reads)
        path(index_dir)
        val(suffix)
        val(remove_sq)
    output:
        tuple val(sample), path("${sample}_${suffix}_minimap2_mapped.sam.gz"), emit: sam
        tuple val(sample), path("${sample}_${suffix}_minimap2_mapped.fastq.gz"), emit: reads_mapped
        tuple val(sample), path("${sample}_${suffix}_minimap2_unmapped.fastq.gz"), emit: reads_unmapped
        tuple val(sample), path("${sample}_${suffix}_minimap2_in.fastq.gz"), emit: input
    shell:
        '''
        set -euo pipefail
        # Prepare inputs
        reads="!{reads}"
        idx="!{index_dir}/mm2_index.mmi"
        sam="!{sample}_!{suffix}_minimap2_mapped.sam.gz"
        al="!{sample}_!{suffix}_minimap2_mapped.fastq.gz"
        un="!{sample}_!{suffix}_minimap2_unmapped.fastq.gz"

        # Run pipeline
        # Outputs a SAM file for all reads, which is then partitioned based on alignment status
        #   - First branch (samtools view -u -f 4 -) filters SAM to unaligned reads and saves FASTQ
        #   - Second branch (samtools view -u -F 4 -) filters SAM to aligned reads and saves FASTQ
        #   - Third branch (samtools view -h -F 4 -) also filters SAM to aligned reads and saves SAM
        zcat ${reads} \
            | minimap2 -a ${idx} /dev/fd/0 \
            | tee \
                >(samtools view -u -f 4 - \
                    | samtools fastq - | gzip -c > ${un}) \
                >(samtools view -u -F 4 - \
                    | samtools fastq - | gzip -c > ${al}) \
            | samtools view -h -F 4 - \
            !{ remove_sq ? "| grep -v '^@SQ'" : "" } | gzip -c > ${sam}
        # Link input to output for testing
        ln -s ${reads} !{sample}_!{suffix}_minimap2_in.fastq.gz
        '''
}