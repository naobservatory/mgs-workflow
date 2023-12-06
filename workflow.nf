/**************
| 0. PREAMBLE |
**************/

// TODO

/****************************
| 1. RAW READ HANDLING & QC |
****************************/

// 1.1. Concatenate split files from same sample together
process CONCAT_GZIPPED {
    cpus 1
    publishDir "${projectDir}/output/raw_concat", mode: "symlink"
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

// 1.2. FASTQC
process FASTQC_CONCAT {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/raw_concat", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(read1), path(read2)
    output:
        path("${sample}_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{read1} !{read2}
        '''
}

// 1.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_CONCAT {
    publishDir "${projectDir}/output/qc/multiqc/raw_concat", mode: "symlink"
    errorStrategy "finish"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

libraries_ch = Channel
    .fromPath(params.library_tab)
    .splitCsv(header: true)
    .map{row -> [row.sample, row.library]}
    .groupTuple()

workflow HANDLE_RAW_READS {
    take:
        libraries_ch
    main:
        concat_ch = CONCAT_GZIPPED(params.raw_dir, libraries_ch)
        fastqc_concat_ch = FASTQC_CONCAT(concat_ch)
        multiqc_concat_ch = MULTIQC_CONCAT(fastqc_concat_ch.collect().ifEmpty([]))
    emit:
        data = concat_ch
        fastqc = fastqc_concat_ch
        multiqc_report = multiqc_concat_ch[0]
        multiqc_data = multiqc_concat_ch[1]
}

/**************************
| 2. TRIMMING & FILTERING |
**************************/

// 2.1. Trim & filter by length & quality
// NB: Does not include deduplication
// TODO: Investigate alternative tools & parameter settings
process PREPROCESS_FASTP {
    cpus 16
    publishDir "${projectDir}/output/cleaned", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(read1), path(read2)
        val adapters
    output:
        tuple val(sample), path("${sample}_fastp_{merged,unmerged_1,unmerged_2}.fastq.gz"), path("${sample}_fastp_failed.fastq.gz"), path("${sample}_fastp.{json,html}")
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
        par="--cut_front --cut_tail --correction --detect_adapter_for_pe --merge --trim_poly_x --average_qual 20 --verbose --dont_eval_duplication --thread !{task.cpus}"
        # Execute
        fastp ${io} ${par}
        '''
}

// 2.2. FASTQC
process FASTQC_CLEANED {
    cpus 3
    publishDir "${projectDir}/output/qc/fastqc/cleaned", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(failed), path(reports)
    output:
        path("${sample}_fastp_{merged,unmerged_1,unmerged_2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads}
        '''
}

// 2.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_CLEANED {
    publishDir "${projectDir}/output/qc/multiqc/cleaned", mode: "symlink"
    errorStrategy "finish"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow CLEAN_READS {
    take:
        concat_ch
    main:
        clean_ch = PREPROCESS_FASTP(concat_ch, params.adapters)
        fastqc_cleaned_ch = FASTQC_CLEANED(clean_ch)
        multiqc_cleaned_ch = MULTIQC_CLEANED(fastqc_cleaned_ch.collect().ifEmpty([]))
    emit:
        data = clean_ch
        fastqc = fastqc_cleaned_ch
        multiqc_report = multiqc_cleaned_ch[0]
        multiqc_data = multiqc_cleaned_ch[1]
}

/*******************
| 3. DEDUPLICATION |
*******************/

// 3.1. Deduplicate with Clumpify
// NB: Will NOT handle RC duplicates
// TODO: Investigate alternative tools & approaches
process DEDUP_CLUMPIFY {
    cpus 16
    publishDir "${projectDir}/output/dedup", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(failed), path(reports)
    output:
        tuple val(sample), path("${sample}_dedup_{unmerged_1,unmerged_2,merged}.fastq.gz")
    shell:
        '''
        # Define input/output
        merged=!{reads[0]}
        unmerged1=!{reads[1]}
        unmerged2=!{reads[2]}
        op1=!{sample}_dedup_unmerged_1.fastq.gz
        op2=!{sample}_dedup_unmerged_2.fastq.gz
        opm=!{sample}_dedup_merged.fastq.gz
        io_unmerged="in=${unmerged1} in2=${unmerged2} out=${op1} out2=${op2}"
        io_merged="in=${merged} out=${opm}"
        # Define parameters
        par="reorder dedupe containment t=!{task.cpus}"
        # Execute
        clumpify.sh ${io_unmerged} ${par}
        clumpify.sh ${io_merged} ${par}
        '''
}

// 3.2. FASTQC
process FASTQC_DEDUP {
    cpus 3
    publishDir "${projectDir}/output/qc/fastqc/dedup", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads)
    output:
        path("${sample}_dedup_{unmerged_1,unmerged_2,merged}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads}
        '''
}

// 3.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_DEDUP {
    publishDir "${projectDir}/output/qc/multiqc/dedup", mode: "symlink"
    errorStrategy "finish"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow DEDUP_READS {
    take:
        clean_ch
    main:
        dedup_ch = DEDUP_CLUMPIFY(clean_ch)
        fastqc_dedup_ch = FASTQC_DEDUP(dedup_ch)
        multiqc_dedup_ch = MULTIQC_DEDUP(fastqc_dedup_ch.collect().ifEmpty([]))
    emit:
        data = dedup_ch
        fastqc = fastqc_dedup_ch
        multiqc_report = multiqc_dedup_ch[0]
        multiqc_data = multiqc_dedup_ch[1]
}

/***************************
| 4. INITIAL RIBODEPLETION |
***************************/

// 4.1. Initial detection and removal of ribosomal reads
// NB: Using quite stringent parameters here to reduce false positives
process BBDUK_RIBO_INITIAL {
    cpus 16
    publishDir "${projectDir}/ribo_initial", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads)
        val ribo_ref
    output:
        tuple val(sample), path("${sample}_bbduk_noribo_{unmerged_1,unmerged_2,merged}.fastq.gz"), path("${sample}_bbduk_ribo_{unmerged_1,unmerged_2,merged}.fastq.gz"), path("${sample}_bbduk_{merged,unmerged}.stats.txt")
    shell:
        '''
        # Define input/output
        merged=!{reads[0]}
        unmerged1=!{reads[1]}
        unmerged2=!{reads[2]}
        op1=!{sample}_bbduk_noribo_unmerged_1.fastq.gz
        op2=!{sample}_bbduk_noribo_unmerged_2.fastq.gz
        opm=!{sample}_bbduk_noribo_merged.fastq.gz
        of1=!{sample}_bbduk_ribo_unmerged_1.fastq.gz
        of2=!{sample}_bbduk_ribo_unmerged_2.fastq.gz
        ofm=!{sample}_bbduk_ribo_merged.fastq.gz
        stats_unmerged=!{sample}_bbduk_unmerged.stats.txt
        stats_merged=!{sample}_bbduk_merged.stats.txt
        ref=!{projectDir}/!{ribo_ref}
        io_unmerged="in=${unmerged1} in2=${unmerged2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats_unmerged}"
        io_merged="in=${merged} ref=${ref} out=${opm} outm=${ofm} stats=${stats_merged}"
        # Define parameters
        par="minkmerfraction=0.5 k=43 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io_unmerged} ${par}
        bbduk.sh ${io_merged} ${par}
        '''
}

// 4.2. FASTQC
process FASTQC_RIBO_INITIAL {
    cpus 3
    publishDir "${projectDir}/output/qc/fastqc/ribo_initial", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
    output:
        path("${sample}_bbduk_noribo_{unmerged_1,unmerged_2,merged}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_noribo}
        '''
}

// 4.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_RIBO_INITIAL {
    publishDir "${projectDir}/output/qc/multiqc/ribo_initial", mode: "symlink"
    errorStrategy "finish"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_RIBO_INITIAL {
    take:
        dedup_ch
    main:
        ribo_ch = BBDUK_RIBO_INITIAL(dedup_ch, params.ribo_ref)
        fastqc_ribo_ch = FASTQC_RIBO_INITIAL(ribo_ch)
        multiqc_ribo_ch = MULTIQC_RIBO_INITIAL(fastqc_ribo_ch.collect().ifEmpty([]))
    emit:
        data = ribo_ch
        fastqc = fastqc_ribo_ch
        multiqc_report = multiqc_ribo_ch[0]
        multiqc_data = multiqc_ribo_ch[1]
}

/********************
| 5. HOST DEPLETION |
********************/

// 5.1. Build human index from reference genome
// NB: Currently uses a masked reference from 2014 from JGI
// TODO: Replace with updated masked human genome
process BBMAP_INDEX_HOST {
    cpus 16
    publishDir "${projectDir}/output/host/index", mode: "symlink"
    errorStrategy "finish"
    input:
        path(reference)
    output:
        tuple path("host_ref.fasta.gz"), path("genome"), path("index")
    shell:
        '''
        cp !{reference} host_ref.fasta.gz
        bbmap.sh ref=host_ref.fasta.gz -Xmx23g
        '''
}

// 5.2. Segregate human reads with bbmap
process BBMAP_HOST_DEPLETION {
    cpus 16
    publishDir "${projectDir}/output/host", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        tuple path(index_fasta), path(index_genome), path(index_index)
    output:
        tuple val(sample), path("${sample}_bbmap_nohost_{unmerged_1,unmerged_2,merged}.fastq.gz"), path("${sample}_bbmap_host_{unmerged_1,unmerged_2,merged}.fastq.gz"), path("${sample}_bbmap_{merged,unmerged}.stats.txt")
    shell:
        '''
        # Define input/output
        merged=!{reads_noribo[0]}
        unmerged1=!{reads_noribo[1]}
        unmerged2=!{reads_noribo[2]}
        op1=!{sample}_bbmap_nohost_unmerged_1.fastq.gz
        op2=!{sample}_bbmap_nohost_unmerged_2.fastq.gz
        opm=!{sample}_bbmap_nohost_merged.fastq.gz
        of1=!{sample}_bbmap_host_unmerged_1.fastq.gz
        of2=!{sample}_bbmap_host_unmerged_2.fastq.gz
        ofm=!{sample}_bbmap_host_merged.fastq.gz
        stats_unmerged=!{sample}_bbmap_unmerged_stats.txt
        stats_merged=!{sample}_bbmap_merged_stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} refstats=${stats_unmerged} path=!{index_fasta}"
        io_merged="in=${merged} outu=${opm} outm=${ofm} refstats=${stats_merged} path=!{index_fasta}"
        # Define parameters (copied from Brian Bushnell)
        par="minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim -Xmx23g t=!{task.cpus}"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        bbmap.sh ${io_merged} ${par}
        '''
}

// 5.3. FASTQC
process FASTQC_HOST {
    cpus 3
    publishDir "${projectDir}/output/qc/fastqc/host", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_nohost), path(reads_host), path(stats)
    output:
        path("${sample}_bbmap_nohost_{unmerged_1,unmerged_2,merged}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_nohost}
        '''
}

// 5.4. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_HOST {
    publishDir "${projectDir}/output/qc/multiqc/host", mode: "symlink"
    errorStrategy "finish"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_HOST {
    take:
        ribo_ch
    main:
        index_ch = BBMAP_INDEX_HOST(params.host_ref)
        host_ch = BBMAP_HOST_DEPLETION(ribo_ch, index_ch)
        fastqc_host_ch = FASTQC_HOST(host_ch)
        multiqc_host_ch = MULTIQC_HOST(fastqc_host_ch.collect().ifEmpty([]))
    emit:
        data = host_ch
        fastqc = fastqc_host_ch
        multiqc_report = multiqc_host_ch[0]
        multiqc_data = multiqc_host_ch[1]
}

workflow {
    // Core pipeline
    HANDLE_RAW_READS(libraries_ch)
    CLEAN_READS(HANDLE_RAW_READS.out.data)
    DEDUP_READS(CLEAN_READS.out.data)
    REMOVE_RIBO_INITIAL(DEDUP_READS.out.data)
    //REMOVE_HOST(REMOVE_RIBO_INITIAL.out.data)
    // TODO: Collate QC
}

/***************************************
| EXTRA PROCESSES, NOT YET IN PIPELINE |
***************************************/

//workflow {
//    // 3. Kraken taxonomic assignment + processing
//    kraken_ch = KRAKEN(dedup_pre_ch, params.kraken_db)
//    filter_ch = EXTRACT_VIRAL_HITS(kraken_ch, params.script_dir, params.virus_db)
//    // X. QC
//}

process KRAKEN {
    cpus 16
    publishDir "${projectDir}/kraken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_noribo), path(merged_ribo)
        val db_path
    output:
        tuple val(sample), path("${sample}_ribo.output"), path("${sample}_noribo.output")
    shell:
        '''
        # Define input/output
        db=!{projectDir}/!{db_path}
        out_unmerged_noribo=!{sample}_unmerged_noribo.output
        out_unmerged_ribo=!{sample}_unmerged_ribo.output
        out_merged_noribo=!{sample}_merged_noribo.output
        out_merged_ribo=!{sample}_merged_ribo.output
        io_unmerged_noribo="--paired !{unmerged_noribo_1} !{unmerged_noribo_2}"
        io_unmerged_ribo="--paired !{unmerged_ribo_1} !{unmerged_ribo_2}"
        io_merged_noribo="!{merged_noribo}"
        io_merged_ribo="!{merged_ribo}"
        # Define parameters
        par="--db ${db} --use-names --threads !{task.cpus}"
        # Run Kraken
        kraken2 ${par} ${io_unmerged_noribo} > ${out_unmerged_noribo}
        kraken2 ${par} ${io_unmerged_ribo} > ${out_unmerged_ribo}
        kraken2 ${par} ${io_merged_noribo} > ${out_merged_noribo}
        kraken2 ${par} ${io_merged_ribo} > ${out_merged_ribo}
        # Concatenate
        out_noribo=!{sample}_noribo.output
        out_ribo=!{sample}_ribo.output
        cat ${out_unmerged_noribo} ${out_merged_noribo} > ${out_noribo}
        cat ${out_unmerged_ribo} ${out_merged_ribo} > ${out_ribo}
        # Cleanup
        rm ${out_unmerged_noribo} ${out_unmerged_ribo} ${out_merged_noribo} ${out_merged_ribo}
        '''
}
//process KRAKEN {
//    cpus 16
//    publishDir "${projectDir}/kraken", mode: "symlink"
//    errorStrategy "finish"
//    input:
//        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_noribo), path(merged_ribo)
//        val db_path
//    output:
//        tuple val(sample), path("${sample}_ribo.output.gz"), path("${sample}_noribo.output.gz")
//    shell:
//        '''
//        # Define input/output
//        db=!{projectDir}/!{db_path}
//        out_unmerged_noribo=!{sample}_unmerged_noribo.output.gz
//        out_unmerged_ribo=!{sample}_unmerged_ribo.output.gz
//        out_merged_noribo=!{sample}_merged_noribo.output.gz
//        out_merged_ribo=!{sample}_merged_ribo.output.gz
//        io_unmerged_noribo="--paired !{unmerged_noribo_1} !{unmerged_noribo_2}"
//        io_unmerged_ribo="--paired !{unmerged_ribo_1} !{unmerged_ribo_2}"
//        io_merged_noribo="!{merged_noribo}"
//        io_merged_ribo="!{merged_ribo}"
//        # Define parameters
//        par="--db ${db} --use-names --threads !{task.cpus}"
//        # Run Kraken
//        kraken2 ${par} ${io_unmerged_noribo} | gzip > ${out_unmerged_noribo}
//        kraken2 ${par} ${io_unmerged_ribo} | gzip > ${out_unmerged_ribo}
//        kraken2 ${par} ${io_merged_noribo} | gzip > ${out_merged_noribo}
//        kraken2 ${par} ${io_merged_ribo} | gzip > ${out_merged_ribo}
//        # Concatenate
//        out_noribo=!{sample}_noribo.output.gz
//        out_ribo=!{sample}_ribo.output.gz
//        cat ${out_unmerged_noribo} ${out_merged_noribo} > ${out_noribo}
//        cat ${out_unmerged_ribo} ${out_merged_ribo} > ${out_ribo}
//        # Cleanup
//        rm ${out_unmerged_noribo} ${out_unmerged_ribo} ${out_merged_noribo} ${out_merged_ribo}
//        '''
//}

/* TOO SLOW - NEED TO OPTIMIZE */
process EXTRACT_VIRAL_HITS {
    cpus 16
    publishDir "${projectDir}/hv_hits", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(kraken_ribo), path(kraken_noribo)
        val scriptDir
        val virusPath
    output:
        tuple val(sample), path("${sample}_hv_hits.output")
    shell:
        '''
        # Define paths
        script_path=!{projectDir}/!{scriptDir}/filter_kraken.py
        virus_path=!{projectDir}/!{virusPath}
        out_ribo="!{sample}_hv_hits_ribo.output"
        out_noribo="!{sample}_hv_hits_noribo.output"
        ${script_path} -t !{task.cpus} !{kraken_ribo} ${virus_path} ${out_ribo}
        ${script_path} -t !{task.cpus} !{kraken_noribo} ${virus_path} ${out_noribo}
        # Concatenate & cleanup
        out="!{sample}_hv_hits.output"
        cat ${out_noribo} ${out_ribo} > ${out}
        rm ${out_ribo} ${out_noribo}
        '''
}

process EXTRACT_VIRAL_READS {
    cpus 16
    publishDir "${projectDir}/hv_reads_putative", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(unmerged_noribo_1), path(unmerged_noribo_2), path(unmerged_ribo_1), path(unmerged_ribo_2), path(merged_noribo), path(merged_ribo)
        tuple val(sample), path(hv_hits)
    output:
        tuple val(sample), path("${sample}_hv_unmerged_noribo_1.fastq.gz"), path("${sample}_hv_unmerged_noribo_2.fastq.gz"), path("${sample}_hv_unmerged_ribo_1.fastq.gz"), path("${sample}_hv_unmerged_ribo_2.fastq.gz"), path("${sample}_hv_merged_noribo.fastq.gz"), path("${sample}_hv_merged_ribo.fastq.gz")
    shell:
        '''
        # Generate names file
        names_path=test
        # Define input/output
        op1=!{sample}_hv_unmerged_noribo_1.fastq.gz
        op2=!{sample}_hv_unmerged_noribo_2.fastq.gz
        opm=!{sample}_hv_merged_noribo.fastq.gz
        of1=!{sample}_hv_unmerged_ribo_1.fastq.gz
        of2=!{sample}_hv_unmerged_ribo_2.fastq.gz
        ofm=!{sample}_hv_merged_ribo.fastq.gz
        io_unmerged_noribo="in=!{unmerged_noribo_1} in2=!{unmerged_noribo_2} out=${op1} out2=${op2}"
        io_unmerged_ribo="in=!{unmerged_ribo_1} in2=!{unmerged_ribo_2} out=${of1} out2=${of2}"
        io_merged_noribo="in=!{merged_noribo} out=${opm}"
        io_merged_ribo="in=!{merged_ribo} out=${ofm}"
        # Define parameters
        par="names=${names_path} include=t t=!{task.cpus}"
        # Execute
        filterbyname.sh ${io_unmerged_noribo} ${par}
        filterbyname.sh ${io_unmerged_ribo} ${par}
        filterbyname.sh ${io_merged_noribo} ${par}
        filterbyname.sh ${io_merged_ribo} ${par}
        '''
}
