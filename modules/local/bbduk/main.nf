// Detection and removal of contaminant reads
process BBDUK_PAIRED {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(min_kmer_fraction)
        val(k)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbduk_pass_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${suffix}_bbduk_fail_{1,2}.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_bbduk.stats.txt"), emit: log
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        op1=!{sample}_!{suffix}_bbduk_pass_1.fastq.gz
        op2=!{sample}_!{suffix}_bbduk_pass_2.fastq.gz
        of1=!{sample}_!{suffix}_bbduk_fail_1.fastq.gz
        of2=!{sample}_!{suffix}_bbduk_fail_2.fastq.gz
        stats=!{sample}_!{suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        io="in=${in1} in2=${in2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats}"
        # Define parameters
        par="minkmerfraction=!{min_kmer_fraction} k=!{k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

process BBDUK_SINGLE {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(min_kmer_fraction)
        val(k)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbduk_pass.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${suffix}_bbduk_fail.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_bbduk.stats.txt"), emit: log
    shell:
        '''
        # Define input/output
        in=!{reads}
        op=!{sample}_!{suffix}_bbduk_pass.fastq.gz
        of=!{sample}_!{suffix}_bbduk_fail.fastq.gz
        stats=!{sample}_!{suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        io="in=${in} ref=${ref} out=${op} outm=${of} stats=${stats}"
        # Define parameters
        par="minkmerfraction=!{min_kmer_fraction} k=!{k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

// Detection and removal of contaminant reads (use minkmerhits instead of minkmerfraction)
process BBDUK_HITS {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(min_kmer_hits)
        val(k)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbduk_pass_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${suffix}_bbduk_fail_{1,2}.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_bbduk.stats.txt"), emit: log
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        op1=!{sample}_!{suffix}_bbduk_pass_1.fastq.gz
        op2=!{sample}_!{suffix}_bbduk_pass_2.fastq.gz
        of1=!{sample}_!{suffix}_bbduk_fail_1.fastq.gz
        of2=!{sample}_!{suffix}_bbduk_fail_2.fastq.gz
        stats=!{sample}_!{suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        io="in=${in1} in2=${in2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats}"
        # Define parameters
        par="minkmerhits=!{min_kmer_hits} k=!{k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

// Streamed version of BBDUK_HITS that returns an interleaved file
process BBDUK_HITS_STREAMED {
    label "small"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(min_kmer_hits)
        val(k)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_bbduk_pass.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${suffix}_bbduk_fail.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_bbduk.stats.txt"), emit: log
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        op=!{sample}_!{suffix}_bbduk_pass.fastq.gz
        of=!{sample}_!{suffix}_bbduk_fail.fastq.gz
        stats=!{sample}_!{suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        io="in=<(cat ${in1}) in2=<(cat ${in2}) ref=${ref} out=${op1} outm=${of1} stats=${stats}"
        # Define parameters
        par="minkmerhits=!{min_kmer_hits} k=!{k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

// Masking contaminant kmers in a sequence database
process BBDUK_MASK {
    label "large"
    label "BBTools"
    input:
        path(seq_db)
        path(contaminant_ref)
        val(k)
        val(label)
    output:
        path("${label}_masked.fasta.gz"), emit: masked
        path("${label}_mask.stats.txt"), emit: log
    shell:
        '''
        # Define input/output
        in=!{seq_db}
        out=!{label}_masked.fasta.gz
        stats=!{label}_mask.stats.txt
        ref=!{contaminant_ref}
        par="k=!{k} hdist=1 mink=8 mm=f rcomp=t maskmiddle=t mask=N t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        bbduk.sh in=${in} out=${out} ref=${ref} stats=${stats} ${par}
        '''
}
