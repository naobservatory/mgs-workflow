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
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/raw_concat", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/raw_concat", mode: "copy", overwrite: "true"
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
    publishDir "${projectDir}/output/preprocess/cleaned", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(read1), path(read2)
        val adapters
    output:
        tuple val(sample), path("${sample}_fastp_{1,2}.fastq.gz"), path("${sample}_fastp_failed.fastq.gz"), path("${sample}_fastp.{json,html}")
    shell:
        '''
        # Define paths and subcommands
        o1=!{sample}_fastp_1.fastq.gz
        o2=!{sample}_fastp_2.fastq.gz
        of=!{sample}_fastp_failed.fastq.gz
        oj=!{sample}_fastp.json
        oh=!{sample}_fastp.html
        ad=!{projectDir}/!{adapters}
        io="--in1 !{read1} --in2 !{read2} --out1 ${o1} --out2 ${o2} --failed_out ${of} --html ${oh} --json ${oj} --adapter_fasta ${ad}"
        par="--cut_front --cut_tail --correction --detect_adapter_for_pe --trim_poly_x --cut_mean_quality 25 --average_qual 25 --qualified_quality_phred 20 --verbose --dont_eval_duplication --thread !{task.cpus} --low_complexity_filter"
        # Execute
        fastp ${io} ${par}
        '''
}

// 2.2. FASTQC
process FASTQC_CLEANED {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/cleaned", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(failed), path(reports)
    output:
        path("${sample}_fastp_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads}
        '''
}

// 2.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_CLEANED {
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/cleaned", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/cleaned", mode: "copy", overwrite: "true"
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
    publishDir "${projectDir}/output/preprocess/dedup", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(failed), path(reports)
    output:
        tuple val(sample), path("${sample}_dedup_{1,2}.fastq.gz")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads[0]}
        unmerged2=!{reads[1]}
        op1=!{sample}_dedup_1.fastq.gz
        op2=!{sample}_dedup_2.fastq.gz
        io_unmerged="in=${unmerged1} in2=${unmerged2} out=${op1} out2=${op2}"
        # Define parameters
        par="reorder dedupe containment t=!{task.cpus}"
        # Execute
        clumpify.sh ${io_unmerged} ${par}
        '''
}

// 3.2. FASTQC
process FASTQC_DEDUP {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/dedup", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads)
    output:
        path("${sample}_dedup_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads}
        '''
}

// 3.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_DEDUP {
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/dedup", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/dedup", mode: "copy", overwrite: "true"
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
    publishDir "${projectDir}/output/preprocess/ribo_initial", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads)
        val ribo_ref
    output:
        tuple val(sample), path("${sample}_bbduk_noribo_{1,2}.fastq.gz"), path("${sample}_bbduk_ribo_{1,2}.fastq.gz"), path("${sample}_bbduk.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads[0]}
        unmerged2=!{reads[1]}
        op1=!{sample}_bbduk_noribo_1.fastq.gz
        op2=!{sample}_bbduk_noribo_2.fastq.gz
        of1=!{sample}_bbduk_ribo_1.fastq.gz
        of2=!{sample}_bbduk_ribo_2.fastq.gz
        stats_unmerged=!{sample}_bbduk.stats.txt
        ref=!{projectDir}/!{ribo_ref}
        io_unmerged="in=${unmerged1} in2=${unmerged2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats_unmerged}"
        # Define parameters
        par="minkmerfraction=0.6 k=43 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io_unmerged} ${par}
        '''
}

// 4.2. FASTQC
process FASTQC_RIBO_INITIAL {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/ribo_initial", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
    output:
        path("${sample}_bbduk_noribo_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_noribo}
        '''
}

// 4.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_RIBO_INITIAL {
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/ribo_initial", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/ribo_initial", mode: "copy", overwrite: "true"
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
    publishDir "${projectDir}/output/preprocess/host/index", mode: "symlink"
    errorStrategy "finish"
    input:
        val(reference)
    output:
        path("ref_index")
    shell:
        '''
        mkdir ref_index
        ref=!{projectDir}/!{reference}
        cp ${ref} ref_index/host_ref.fasta.gz
        cd ref_index
        bbmap.sh ref=host_ref.fasta.gz -Xmx23g t=!{task.cpus}
        '''
}

// 5.2. Segregate human reads with bbmap
// TODO: Consider loosening parameters
process BBMAP_HOST_DEPLETION {
    cpus 16
    publishDir "${projectDir}/output/preprocess/host", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        path(index_ref_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_nohost_{1,2}.fastq.gz"), path("${sample}_bbmap_host_{1,2}.fastq.gz"), path("${sample}_bbmap.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads_noribo[0]}
        unmerged2=!{reads_noribo[1]}
        op1=!{sample}_bbmap_nohost_1.fastq.gz
        op2=!{sample}_bbmap_nohost_2.fastq.gz
        of1=!{sample}_bbmap_host_1.fastq.gz
        of2=!{sample}_bbmap_host_2.fastq.gz
        stats=!{sample}_bbmap.stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} statsfile=${stats} path=!{index_ref_dir}"
        # Define parameters (copied from Brian Bushnell)
        par="minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim -Xmx23g t=!{task.cpus}"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        '''
}

// 5.3. FASTQC
process FASTQC_HOST {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/host", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_nohost), path(reads_host), path(stats)
    output:
        path("${sample}_bbmap_nohost_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_nohost}
        '''
}

// 5.4. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_HOST {
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/host", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/host", mode: "copy", overwrite: "true"
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

/****************************************
| 6. HUMAN VIRUS DETECTION WITH BOWTIE2 |
****************************************/

// TODO: Write processes for Bowtie2 index construction, alignment, data parsing
// TODO: Write workflow

// 6.1. Mask collated virus genomes
process MASK_HV_GENOMES {
    conda "${projectDir}/${params.env_dir}/bowtie.yaml"
    cpus 1
    publishDir "${projectDir}/output/hviral/index", mode: "symlink"
    errorStrategy "finish"
    input:
        path(collected_genomes)
    output:
        path("masked_genomes.fasta.gz")
    shell:
        '''
        zcat -f !{collected_genomes} | dustmasker -out "masked_genomes.fasta" -outfmt fasta
        sed -i '/^>/!s/[a-z]/x/g' masked_genomes.fasta
        gzip masked_genomes.fasta
        '''
}

// 6.2. Build Bowtie2 index from masked genomes
process BUILD_BOWTIE2_DB {
    conda "${projectDir}/${params.env_dir}/bowtie.yaml"
    cpus 16
    publishDir "${projectDir}/output/hviral/index", mode: "symlink"
    errorStrategy "finish"
    input:
        path(masked_genomes)
    output:
        path("bt2_hv_index") // Output directory for index files
    shell:
        '''
        mkdir bt2_hv_index
        bowtie2-build -f --threads !{task.cpus} !{masked_genomes} bt2_hv_index/hv_index
        '''
}

// 6.3. Run Bowtie2 and return mapped HV reads
// TODO: Return unmapped reads for downstream processing? Currently assuming these are rare enough not to matter.
// TODO: If return unmapped reads, add FASTQC/MultiQC processes
process RUN_BOWTIE2 {
    conda "${projectDir}/${params.env_dir}/bowtie.yaml"
    cpus 16
    publishDir "${projectDir}/output/hviral/bowtie", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_nohost), path(reads_host), path(stats)
        path(index_dir)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped.sam"), path("${sample}_bowtie2_mapped_{1,2}.fastq.gz")
    shell:
        '''
        in1=!{reads_nohost[0]}
        in2=!{reads_nohost[1]}
        idx="!{index_dir}/hv_index"
        sam="!{sample}_bowtie2_mapped.sam"
        alc="!{sample}_bowtie2_mapped_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc}"
        par="--threads !{task.cpus} --no-unal --no-sq --local --very-sensitive-local --score-min G,1,0 --mp 4,1"
        bowtie2 ${par} ${io}
        '''
}

// 6.4. Merge-join deduplicated Bowtie2 output for Kraken processing
process MERGE_JOIN_BOWTIE {
    cpus 1
    publishDir "${projectDir}/output/hviral/merged", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(sam_out), path(reads_out)
        val scriptDir
    output:
        tuple val(sample), path("${sample}_bowtie2_mjc.fastq.gz"), path("${sample}_bowtie2_bbmerge_stats.txt")
    shell:
        '''
        # Prepare input/output for bbmerge
        in1=!{reads_out[0]}
        in2=!{reads_out[1]}
        ou1=!{sample}_bowtie2_bbmerge_unmerged_1.fastq.gz
        ou2=!{sample}_bowtie2_bbmerge_unmerged_2.fastq.gz
        om=!{sample}_bowtie2_bbmerge_merged.fastq.gz
        stats=!{sample}_bowtie2_bbmerge_stats.txt
        io="in=${in1} in2=${in2} out=${om} outu=${ou1} outu2=${ou2} ihist=${stats}"
        # Execute bbmerge
        bbmerge.sh ${io}
        # Prepare to join unmerged read pairs
        script_path=!{projectDir}/!{scriptDir}/join_fastq.py
        oj=!{sample}_bowtie2_bbmerge_unmerged_joined.fastq.gz
        # Join unmerged read pairs
        ${script_path} ${ou1} ${ou2} ${oj}
        # Concatenate single output file
        oo=!{sample}_bowtie2_mjc.fastq.gz
        cat ${om} ${oj} > ${oo}
        '''
}

// 6.5. Perform taxonomic assignment with Kraken2
process KRAKEN_BOWTIE {
    cpus 16
    publishDir "${projectDir}/output/hviral/kraken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(stats) // Single input file, merged, joined & concatenated
        val db_path
    output:
        tuple val(sample), path("${sample}_bowtie2_kraken.output"), path("${sample}_bowtie2_kraken.report")
    shell:
        '''
        # Define input/output
        db=!{projectDir}/!{db_path}
        in=!{reads}
        out=!{sample}_bowtie2_kraken.output
        report=!{sample}_bowtie2_kraken.report
        io="--output ${out} --report ${report} ${in}"
        # Define parameters
        par="--db ${db} --use-names --threads !{task.cpus}"
        # Run Kraken
        kraken2 ${par} ${io}
        '''
}

// 6.6. Process Kraken2 output and identify HV- and non-HV-assigned reads
process PROCESS_KRAKEN_BOWTIE {
    cpus 1
    publishDir "${projectDir}/output/hviral/kraken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(output), path(report)
        val scriptDir
        val nodes_path
        val hv_path
    output:
        tuple val(sample), path("${sample}_bowtie2_kraken_processed.tsv")
    shell:
        '''
        in=!{output}
        out=!{sample}_bowtie2_kraken_processed.tsv
        nodes=!{projectDir}/!{nodes_path}
        hv=!{projectDir}/!{hv_path}
        script_path=!{projectDir}/!{scriptDir}/process_kraken_hv.py
        ${script_path} ${in} ${hv} ${nodes} ${out}
        '''
}

// 6.7. Process Bowtie2 SAM output
// NB: Currently paired, need to update if switch to merged
process PROCESS_BOWTIE_SAM {
    cpus 1
    publishDir "${projectDir}/output/hviral/bowtie", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(sam), path(reads)
        val scriptDir
        val genomeid_taxid_map_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv")
    shell:
        '''
        in=!{sam}
        out=!{sample}_bowtie2_sam_processed.tsv
        map=!{projectDir}/!{genomeid_taxid_map_path}
        script_path=!{projectDir}/!{scriptDir}/process_sam_hv.py
        ${script_path} ${in} ${map} ${out}
        '''
}

// 6.8. Merge processed SAM and Kraken TSVs and compute length-normalized alignment scores
process MERGE_SAM_KRAKEN {
    cpus 1
    publishDir "${projectDir}/output/hviral/hits", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(kraken_processed), path(sam_processed)
    output:
        path("${sample}_hv_hits_putative.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        sam <- read_tsv("!{sam_processed}", show_col_types = FALSE)
        krk <- read_tsv("!{kraken_processed}", show_col_types = FALSE)
        mrg <- sam %>% rename(seq_id = query_name) %>% inner_join(krk, by="seq_id") %>%
            mutate(adj_score_fwd = best_alignment_score_fwd/log(query_len_fwd),
                   adj_score_rev = best_alignment_score_rev/log(query_len_rev),
                   sample="!{sample}")
        cat(nrow(sam), nrow(krk), nrow(mrg))
        write_tsv(mrg, "!{sample}_hv_hits_putative.tsv")
        '''
}

// 6.9. Collapse outputs from different samples into one TSV
process MERGE_SAMPLES_HV {
    cpus 1
    publishDir "${projectDir}/output/hviral/hits", mode: "symlink"
    publishDir "${projectDir}/output/results", mode: "copy", overwrite: "true"
    errorStrategy "finish"
    input:
        path(tsvs)
    output:
        path("hv_hits_putative.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        in_paths <- str_split("!{tsvs}", " ")[[1]]
        print(in_paths)
        tabs <- lapply(in_paths, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        for (t in tabs) print(dim(t))
        tab_out <- bind_rows(tabs)
        print(dim(tab_out))
        sapply(tabs, nrow) %>% sum %>% print
        write_tsv(tab_out, "hv_hits_putative.tsv")
        '''
}

workflow MAP_HUMAN_VIRUSES {
    take:
        host_ch
    main:
        mask_ch = MASK_HV_GENOMES("${projectDir}/${params.hv_genomes}")
        index_ch = BUILD_BOWTIE2_DB(mask_ch)
        bowtie2_ch = RUN_BOWTIE2(host_ch, index_ch)
        merge_ch = MERGE_JOIN_BOWTIE(bowtie2_ch, params.script_dir)
        kraken_ch = KRAKEN_BOWTIE(merge_ch, params.kraken_db)
        kraken_processed_ch = PROCESS_KRAKEN_BOWTIE(kraken_ch, params.script_dir, params.nodes, params.virus_db)
        bowtie2_processed_ch = PROCESS_BOWTIE_SAM(bowtie2_ch, params.script_dir, params.genomeid_map)
        merged_input_ch = kraken_processed_ch.combine(bowtie2_processed_ch, by: 0)
        merged_ch = MERGE_SAM_KRAKEN(merged_input_ch)
        merged_ch_2 = MERGE_SAMPLES_HV(merged_ch.collect().ifEmpty([]))
    emit:
        data = merged_ch_2
}

/*****************************
| 7. SECONDARY RIBODEPLETION |
*****************************/

// 7.1. Secondary detection and removal of ribosomal reads
// NB: Using more liberal parameters here since very high specificity less important
// TODO: Adjust input structure after writing section 6
process BBDUK_RIBO_SECONDARY {
    cpus 16
    publishDir "${projectDir}/output/preprocess/ribo_secondary", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_nohost), path(reads_host), path(stats)
        val ribo_ref
    output:
        tuple val(sample), path("${sample}_bbduk_noribo_{1,2}.fastq.gz"), path("${sample}_bbduk_ribo_{1,2}.fastq.gz"), path("${sample}_bbduk.stats.txt")
    shell:
        '''
        # Define input/output
        in1=!{reads_nohost[0]}
        in2=!{reads_nohost[1]}
        op1=!{sample}_bbduk_noribo_1.fastq.gz
        op2=!{sample}_bbduk_noribo_2.fastq.gz
        of1=!{sample}_bbduk_ribo_1.fastq.gz
        of2=!{sample}_bbduk_ribo_2.fastq.gz
        stats=!{sample}_bbduk.stats.txt
        ref=!{projectDir}/!{ribo_ref}
        io="in=${in1} in2=${in2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats}"
        # Define parameters
        par="minkmerfraction=0.4 k=27 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

// 7.2. FASTQC
process FASTQC_RIBO_SECONDARY {
    cpus 2
    publishDir "${projectDir}/output/qc/fastqc/ribo_secondary", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
    output:
        path("${sample}_bbduk_noribo_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_noribo}
        '''
}

// 7.3. MultiQC
// TODO: Customize output for later aggregation
process MULTIQC_RIBO_SECONDARY {
    cpus 1
    publishDir "${projectDir}/output/qc/multiqc/ribo_secondary", mode: "symlink"
    publishDir "${projectDir}/output/results/multiqc/ribo_secondary", mode: "copy", overwrite: "true"
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

workflow REMOVE_RIBO_SECONDARY {
    take:
        novirus_ch
    main:
        ribo_ch = BBDUK_RIBO_SECONDARY(novirus_ch, params.ribo_ref)
        fastqc_ribo_ch = FASTQC_RIBO_SECONDARY(ribo_ch)
        multiqc_ribo_ch = MULTIQC_RIBO_SECONDARY(fastqc_ribo_ch.collect().ifEmpty([]))
    emit:
        data = ribo_ch
        fastqc = fastqc_ribo_ch
        multiqc_report = multiqc_ribo_ch[0]
        multiqc_data = multiqc_ribo_ch[1]
}

/**************************************
| 8. TAXONOMIC ASSIGNMENT WITH KRAKEN |
**************************************/

// 8.1. Merge-join reads for Kraken processing
process MERGE_JOIN {
    cpus 1
    publishDir "${projectDir}/output/taxonomy/merged", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        val scriptDir
    output:
        tuple val(sample), path("${sample}_mjc.fastq.gz"), path("${sample}_bbmerge_stats.txt")
    shell:
        '''
        # Prepare input/output for bbmerge
        in1=!{reads_noribo[0]}
        in2=!{reads_noribo[1]}
        ou1=!{sample}_bbmerge_unmerged_1.fastq.gz
        ou2=!{sample}_bbmerge_unmerged_2.fastq.gz
        om=!{sample}_bbmerge_merged.fastq.gz
        stats=!{sample}_bbmerge_stats.txt
        io="in=${in1} in2=${in2} out=${om} outu=${ou1} outu2=${ou2} ihist=${stats}"
        # Execute bbmerge
        bbmerge.sh ${io}
        # Prepare to join unmerged read pairs
        script_path=!{projectDir}/!{scriptDir}/join_fastq.py
        oj=!{sample}_bbmerge_unmerged_joined.fastq.gz
        # Join unmerged read pairs
        ${script_path} ${ou1} ${ou2} ${oj}
        # Concatenate single output file
        oo=!{sample}_mjc.fastq.gz
        cat ${om} ${oj} > ${oo}
        '''
}

// 8.2. Perform taxonomic assignment with Kraken2
// TODO: Check & update unclassified_out file configuration
process KRAKEN {
    cpus 16
    publishDir "${projectDir}/output/taxonomy/kraken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(reads), path(stats)
        val db_path
    output:
        tuple val(sample), path("${sample}.output"), path("${sample}.report"), path("${sample}_unclassified.fastq.gz")
    shell:
        '''
        # Define input/output
        db=!{projectDir}/!{db_path}
        in=!{reads}
        out=!{sample}.output
        report=!{sample}.report
        unc=!{sample}_unclassified.fastq
        io="--output ${out} --report ${report} --unclassified-out ${unc} ${in}"
        # Define parameters
        par="--db ${db} --use-names --threads !{task.cpus}"
        # Run Kraken
        kraken2 ${par} ${io}
        # Gzip output
        gzip ${unc}
        '''
}

// 8.3. Summarize Kraken output with Bracken
process BRACKEN_DOMAINS {
    cpus 1
    conda "${projectDir}/${params.env_dir}/bracken.yaml"
    publishDir "${projectDir}/output/taxonomy/bracken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(output), path(report), path(unc_reads)
        val db_path
    output:
        tuple val(sample), path("${sample}.bracken")
    shell:
        '''
        # Define input/output
        db=!{projectDir}/!{db_path}
        in=!{report}
        out=!{sample}.bracken
        io="-d ${db} -i ${in} -o ${out}"
        # Define parameters
        par="-l D"
        # Run Bracken
        bracken ${io} ${par}
        '''
}

// 8.4. Label Bracken files with sample IDs
process LABEL_BRACKEN {
    cpus 1
    publishDir "${projectDir}/output/taxonomy/bracken", mode: "symlink"
    errorStrategy "finish"
    input:
        tuple val(sample), path(bracken_output)
    output:
        path("${sample}_labeled.bracken")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        tab <- read_tsv("!{bracken_output}", show_col_types = FALSE) %>%
            mutate(sample="!{sample}")
        write_tsv(tab, "!{sample}_labeled.bracken")
        '''
}

// 8.5. Combine Bracken files into a single output file
process MERGE_BRACKEN {
    cpus 1
    publishDir "${projectDir}/output/taxonomy/results", mode: "symlink"
    publishDir "${projectDir}/output/results", mode: "copy", overwrite: "true"
    errorStrategy "finish"
    input:
        path(tsvs)
    output:
        path("bracken_counts.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        in_paths <- str_split("!{tsvs}", " ")[[1]]
        print(in_paths)
        tabs <- lapply(in_paths, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        for (t in tabs) print(dim(t))
        tab_out <- bind_rows(tabs)
        print(dim(tab_out))
        sapply(tabs, nrow) %>% sum %>% print
        write_tsv(tab_out, "bracken_counts.tsv")
        '''
}

workflow CLASSIFY_READS {
    take:
        data_ch
    main:
        merged_ch = MERGE_JOIN(data_ch, params.script_dir)
        kraken_ch = KRAKEN(merged_ch, params.kraken_db)
        bracken_ch = BRACKEN_DOMAINS(kraken_ch, params.kraken_db)
        label_ch = LABEL_BRACKEN(bracken_ch)
        merged_ch_2 = MERGE_BRACKEN(label_ch.collect().ifEmpty([]))
    emit:
        kraken = kraken_ch
        bracken = merged_ch_2
}

workflow {
    // Preprocessing
    HANDLE_RAW_READS(libraries_ch)
    CLEAN_READS(HANDLE_RAW_READS.out.data)
    DEDUP_READS(CLEAN_READS.out.data)
    REMOVE_RIBO_INITIAL(DEDUP_READS.out.data)
    REMOVE_HOST(REMOVE_RIBO_INITIAL.out.data)
    REMOVE_RIBO_SECONDARY(REMOVE_HOST.out.data)
    // TODO: Collate QC metrics from preprocessing
    // Human viral reads
    MAP_HUMAN_VIRUSES(REMOVE_HOST.out.data)
    // Broad taxonomic profiling
    CLASSIFY_READS(REMOVE_RIBO_SECONDARY.out.data)
}
