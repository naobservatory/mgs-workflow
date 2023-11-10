
/*****************
| CORE PROCESSES |
*****************/

process CONCAT_GZIPPED {
    cpus 1
    publishDir "${projectDir}/raw_concat", mode: "symlink"
    errorStrategy "finish"
    input:
        val raw_files_directory
        tuple val(sample), val(libraries)
    output:
        tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz")
    shell:
        '''
        # Preamble
        read_dir=!{projectDir}/!{raw_files_directory}
        echo Raw files directory: $read_dir
        echo Sample: !{sample}
        echo Libraries: !{libraries.join(" ")}
        # Get file paths from library IDs
        r1=""
        r2=""
        for l in !{libraries.join(" ")}; do
            L1=$(ls ${read_dir}/*${l}*_1.fastq.gz);
            L2=$(ls ${read_dir}/*${l}*_2.fastq.gz);
            r1="${r1} ${L1}"
            r2="${r2} ${L2}"
            done
        # Concatenate
        echo Read 1 files to concatenate: ${r1}
        cat ${r1} > !{sample}_1.fastq.gz
        echo Read 2 files to concatenate: ${r2}
        cat ${r2} > !{sample}_2.fastq.gz
        '''
}

process PREPROCESS_FASTP {
    cpus 16
    publishDir "${projectDir}/preproc", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(read1), path(read2)
        val adapters
    output:
        tuple val(sample), path("${sample}_fastp_unmerged_1.fastq.gz"), path("${sample}_fastp_unmerged_2.fastq.gz"), path("${sample}_fastp_merged.fastq.gz"), path("${sample}_fastp_failed.fastq.gz"), path("${sample}_fastp.json"), path("${sample}_fastp.html")
    shell:
        '''
        # Define paths and subcommands
        o1=!{sample}_fastp_unmerged_1.fastq.gz
        o2=!{sample}_fastp_unmerged_2.fastq.gz
        om=!{sample}_fastp_merged.fastq.gz
        of=!{sample}_fastp_failed.fastq.gz
        oj=!{sample}_fastp.json
        oh=!{sample}_fastp.html
        ad=!{projectDir}/!{adapters}
        io="--in1 !{read1} --in2 !{read2} --out1 ${o1} --out2 ${o2} --merged_out ${om} --failed_out ${of} --html ${oh} --json ${oj} --adapter_fasta ${ad}"
        par="--cut_tail --correction --detect_adapter_for_pe --merge --trim_poly_x --verbose --dont_eval_duplication --thread !{task.cpus}"
        # Execute
        fastp ${io} ${par}
        '''
}

process RIBO_BBDUK {
    cpus 16
    publishDir "${projectDir}/ribo", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged1), path(unmerged2), path(merged), path(failed), path(json), path(html)
        val ribo_ref
    output:
        tuple val(sample), path("${sample}_bbduk_unmerged_noribo_1.fastq.gz"), path("${sample}_bbduk_unmerged_noribo_2.fastq.gz"), path("${sample}_bbduk_merged_noribo.fastq.gz"), path("${sample}_bbduk_unmerged_ribo_1.fastq.gz"), path("${sample}_bbduk_unmerged_ribo_2.fastq.gz"), path("${sample}_bbduk_merged_ribo.fastq.gz"), path("${sample}_bbduk_unmerged_stats.txt"), path("${sample}_bbduk_merged_stats.txt")
    shell:
        '''
        # Define input/output
        op1=!{sample}_bbduk_unmerged_noribo_1.fastq.gz
        op2=!{sample}_bbduk_unmerged_noribo_2.fastq.gz
        opm=!{sample}_bbduk_merged_noribo.fastq.gz
        of1=!{sample}_bbduk_unmerged_ribo_1.fastq.gz
        of2=!{sample}_bbduk_unmerged_ribo_2.fastq.gz
        ofm=!{sample}_bbduk_merged_ribo.fastq.gz
        stats_unmerged=!{sample}_bbduk_unmerged_stats.txt
        stats_merged=!{sample}_bbduk_merged_stats.txt
        ref=!{projectDir}/!{ribo_ref}
        io_unmerged="in=!{unmerged1} in2=!{unmerged2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats_unmerged}"
        io_merged="in=!{merged} ref=${ref} out=${opm} outm=${ofm} stats=${stats_merged}"
        # Define parameters
        par="k=33 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io_unmerged} ${par}
        bbduk.sh ${io_merged} ${par}
        '''
}

process DEDUP_CLUMPIFY_PRE {
    cpus 16
    publishDir "${projectDir}/dedup_pre", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(merged_noribo), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_ribo), path(unmerged_stats), path(merged_stats)
    output:
        tuple val(sample), path("${sample}_clumpify_unmerged_noribo_1.fastq.gz"), path("${sample}_clumpify_unmerged_noribo_2.fastq.gz"), path("${sample}_clumpify_unmerged_ribo_1.fastq.gz"), path("${sample}_clumpify_unmerged_ribo_2.fastq.gz"), path("${sample}_clumpify_merged_noribo.fastq.gz"), path("${sample}_clumpify_merged_ribo.fastq.gz")
    shell:
        '''
        # Define input/output
        op1=!{sample}_clumpify_unmerged_noribo_1.fastq.gz
        op2=!{sample}_clumpify_unmerged_noribo_2.fastq.gz
        opm=!{sample}_clumpify_merged_noribo.fastq.gz
        of1=!{sample}_clumpify_unmerged_ribo_1.fastq.gz
        of2=!{sample}_clumpify_unmerged_ribo_2.fastq.gz
        ofm=!{sample}_clumpify_merged_ribo.fastq.gz
        io_unmerged_noribo="in=!{unmerged_noribo_1} in2=!{unmerged_noribo_2} out=${op1} out2=${op2}"
        io_unmerged_ribo="in=!{unmerged_ribo_1} in2=!{unmerged_ribo_2} out=${of1} out2=${of2}"
        io_merged_noribo="in=!{merged_noribo} out=${opm}"
        io_merged_ribo="in=!{merged_ribo}, out=${ofm}"
        # Define parameters
        par="reorder dedup containment t=!{task.cpus}"
        # Execute
        clumpify.sh ${io_unmerged_noribo} ${par}
        clumpify.sh ${io_unmerged_ribo} ${par}
        clumpify.sh ${io_merged_noribo} ${par}
        clumpify.sh ${io_merged_ribo} ${par}
        '''
}

process KRAKEN {
    cpus 16
    publishDir "${projectDir}/dedup_pre", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_noribo), path(merged_ribo)
        val db_path
    output:
        tuple val(sample), path("${sample}_ribo.output.gz"), path("${sample}_noribo.output.gz")
    shell:
        '''
        # Define input/output
        db=!{projectDir}/!{db_path}
        out_unmerged_noribo=!{sample}_unmerged_noribo.output
        out_unmerged_ribo=!{sample}_unmerged_ribo.output
        out_merged_noribo=!{sample}_merged_noribo.output
        out_merged_ribo=!{sample}_merged_ribo.output
        io_unmerged_noribo="--paired --output ${out_unmerged_noribo} !{unmerged_noribo_1} !{unmerged_noribo_2}"
        io_unmerged_ribo="--paired --output ${out_unmerged_ribo} !{unmerged_ribo_1} !{unmerged_ribo_2}"
        io_merged_noribo="--paired --output ${out_merged_noribo} !{merged_noribo}
        io_merged_ribo="--paired --output ${out_merged_ribo} !{merged_ribo}
        # Define parameters
        par="--db ${db} --use-names --threads !{task.cps}"
        # Run Kraken
        kraken2 ${par} ${io_unmerged_noribo}
        kraken2 ${par} ${io_unmerged_ribo}
        kraken2 ${par} ${io_merged_noribo}
        kraken2 ${par} ${io_merged_ribo}
        # Concatenate
        out_noribo=!{sample}_noribo.output.gz
        out_ribo=!{sample}_ribo.output.gz
        cat ${out_unmerged_noribo} ${out_merged_noribo} | gzip > ${out_noribo}
        cat ${out_unmerged_ribo} ${out_merged_ribo} | gzip > ${out_ribo}
        # Cleanup
        rm ${out_unmerged_noribo} ${out_unmerged_ribo} ${out_merged_noribo} ${out_merged_ribo}
        '''
}

/***************
| QC PROCESSES |
***************/

process FASTQC_CONCAT {
    cpus 2
    publishDir "${projectDir}/qc", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(read1), path(read2)
    output:
        tuple val(sample), path("${sample}_1_fastqc.zip"), path("${sample}_1_fastqc.html"), path("${sample}_2_fastqc.zip"), path("${sample}_2_fastqc.html")
    shell:
        '''
        fastqc -t !{task.cpus} !{read1} !{read2}
        '''
}

process FASTQC_PREPROC {
    cpus 3
    publishDir "${projectDir}/qc", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged1), path(unmerged2), path(merged), path(failed), path(json), path(html)
    output:
        tuple val(sample), path("${sample}_fastp_unmerged_1_fastqc.zip"), path("${sample}_fastp_unmerged_1_fastqc.html"), path("${sample}_fastp_unmerged_2_fastqc.zip"), path("${sample}_fastp_unmerged_2_fastqc.html"), path("${sample}_fastp_merged_fastqc.zip"), path("${sample}_fastp_merged_fastqc.html")
    shell:
        '''
        fastqc -t !{task.cpus} !{unmerged1} !{unmerged2} !{merged}
        '''
}

process FASTQC_RIBO {
    cpus 6
    publishDir "${projectDir}/qc", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(merged_noribo), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_ribo), path(unmerged_stats), path(merged_stats)
    output:
        tuple val(sample), path("${sample}_bbduk_unmerged_noribo_1_fastqc.zip"), path("${sample}_bbduk_unmerged_noribo_1_fastqc.html"), path("${sample}_bbduk_unmerged_noribo_2_fastqc.zip"), path("${sample}_bbduk_unmerged_noribo_2_fastqc.html"), path("${sample}_bbduk_merged_noribo_fastqc.zip"), path("${sample}_bbduk_merged_noribo_fastqc.html"), path("${sample}_bbduk_unmerged_ribo_1_fastqc.zip"), path("${sample}_bbduk_unmerged_ribo_1_fastqc.html"), path("${sample}_bbduk_unmerged_ribo_2_fastqc.zip"), path("${sample}_bbduk_unmerged_ribo_2_fastqc.html"), path("${sample}_bbduk_merged_ribo_fastqc.zip"), path("${sample}_bbduk_merged_ribo_fastqc.html")
    shell:
        '''
        fastqc -t !{task.cpus} !{unmerged_noribo_1} !{unmerged_noribo_2} !{merged_noribo} !{unmerged_ribo_1} !{unmerged_ribo_2} !{merged_ribo}
        '''
}

process FASTQC_DEDUP_PRE {
    cpus 6
    publishDir "${projectDir}/qc", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_noribo), path(merged_ribo)
    output:
        tuple val(sample), path("${sample}_clumpify_unmerged_noribo_1_fastqc.zip"), path("${sample}_clumpify_unmerged_noribo_1_fastqc.html"), path("${sample}_clumpify_unmerged_noribo_2_fastqc.zip"), path("${sample}_clumpify_unmerged_noribo_2_fastqc.html"), path("${sample}_clumpify_merged_noribo_fastqc.zip"), path("${sample}_clumpify_merged_noribo_fastqc.html"), path("${sample}_clumpify_unmerged_ribo_1_fastqc.zip"), path("${sample}_clumpify_unmerged_ribo_1_fastqc.html"), path("${sample}_clumpify_unmerged_ribo_2_fastqc.zip"), path("${sample}_clumpify_unmerged_ribo_2_fastqc.html"), path("${sample}_clumpify_merged_ribo_fastqc.zip"), path("${sample}_clumpify_merged_ribo_fastqc.html")
    shell:
        '''
        fastqc -t !{task.cpus} !{unmerged_noribo_1} !{unmerged_noribo_2} !{merged_noribo} !{unmerged_ribo_1} !{unmerged_ribo_2} !{merged_ribo}
        '''
}

/***********
| WORKFLOW |
***********/

workflow {
    // 1. Concatenation
    libraries_ch = Channel
        .fromPath(params.library_tab)
        .splitCsv(header: true)
        .map{row -> [row.sample, row.library]}
        .groupTuple()
    concat_ch = CONCAT_GZIPPED(params.raw_dir, libraries_ch)
    // 2. Preprocessing, ribodetection & deduplication
    preproc_ch = PREPROCESS_FASTP(concat_ch, params.adapters)
    ribo_ch = RIBO_BBDUK(preproc_ch, params.ribo_ref)
    dedup_pre_ch = DEDUP_CLUMPIFY_PRE(ribo_ch)
    // 3. Kraken taxonomic assignment + processing
    kraken_ch = KRAKEN(dedup_pre_ch, params.kraken_db)
    // X. QC
    qc_concat_ch = FASTQC_CONCAT(concat_ch)
    qc_preproc_ch = FASTQC_PREPROC(preproc_ch)
    qc_ribo_ch = FASTQC_RIBO(ribo_ch)
    qc_dedup_pre_ch = FASTQC_DEDUP_PRE(dedup_pre_ch)
}
