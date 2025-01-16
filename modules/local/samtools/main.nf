// Return reads that did not align to reference as FASTQ
process SAMTOOLS_FILTER {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        out=!{sample}_!{suffix}.fastq.gz
        var="fastq -n -f 4"
        samtools ${var} ${in} > ${out}
        '''
}
// Return reads that aligned to reference as FASTQ
process SAMTOOLS_KEEP {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        out=!{sample}_!{suffix}.fastq.gz
        var="fastq -n -F 4"
        samtools ${var} ${in} > ${out}
        '''
}

// Return reads that aligned to reference as SAM
process SAMTOOLS_KEEP_AS_SAM {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.sam"), emit: sam
    shell:
        '''
        in=!{sam}
        out=!{sample}_!{suffix}.sam
        var="view -h -F 4"
        samtools ${var} ${in} > ${out}
        '''
}

// Return matching and non-matching reads separately
process SAMTOOLS_SEPARATE {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_samtools_fail.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_samtools_pass.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        of=!{sample}_!{suffix}_samtools_fail.fastq.gz
        op=!{sample}_!{suffix}_samtools_pass.fastq.gz
        of_var="fastq -n -F 4"
        op_var="fastq -n -f 4"
        samtools ${of_var} ${in} > ${of}
        samtools ${op_var} ${in} > ${op}
        '''
}


// Merge SAM files and return as SAM | FIX!!!
process SAMTOOLS_MERGE {
    label "samtools"
    input:
        tuple val(sample), path(sam_files)
        val(prefix)
    output:
        path("${prefix}.sam")
    shell:
        '''
        # Get list of input files
        in_paths=(${sam_files})

        # Merge SAM files
        samtools merge -h "${name}.sam" ${in_paths}
        '''
}