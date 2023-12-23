/**************
| 0. PREAMBLE |
**************/

// Set up directories
rawDir = "${projectDir}/${params.raw_dir}"
envDir = "${projectDir}/${params.env_dir}"
scriptDir = "${projectDir}/${params.script_dir}"
pubDir = "${projectDir}/${params.pub_dir}"

// Generate paths from param files
adapter_ref_path = "${projectDir}/${params.adapters}"
ribo_ref_path = "${projectDir}/${params.ribo_ref}"
hv_genomes = "${projectDir}/${params.hv_genomes}"
kraken_db_path = "${projectDir}/${params.kraken_db}"
virus_db_path = "${projectDir}/${params.virus_db}"
nodes_path = "${projectDir}/${params.nodes}"
genomeid_map_path = "${projectDir}/${params.genomeid_map}"

/****************************
| 1. RAW READ HANDLING & QC |
****************************/

// 1.1. Concatenate split files from same sample together
process CONCAT_GZIPPED {
    label "single"
    publishDir "${pubDir}/raw_concat", mode: "symlink"
    input:
        val raw_files_directory
        tuple val(sample), val(libraries)
    output:
        tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz")
    shell:
        '''
        # Preamble
        read_dir=!{rawDir}
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
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/raw_concat", mode: "symlink"
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
process MULTIQC_CONCAT {
    conda "${envDir}/qc.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/raw_concat", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("raw_concat"), path("multiqc_data")
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
        concat_ch = CONCAT_GZIPPED(rawDir, libraries_ch)
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
process PREPROCESS_FASTP {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/cleaned", mode: "symlink"
    input:
        tuple val(sample), path(read1), path(read2)
        path adapters
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
        ad=!{adapters}
        io="--in1 !{read1} --in2 !{read2} --out1 ${o1} --out2 ${o2} --failed_out ${of} --html ${oh} --json ${oj} --adapter_fasta ${ad}"
        par="--cut_front --cut_tail --correction --detect_adapter_for_pe --trim_poly_x --cut_mean_quality 25 --average_qual 25 --qualified_quality_phred 20 --verbose --dont_eval_duplication --thread !{task.cpus} --low_complexity_filter"
        # Execute
        fastp ${io} ${par}
        '''
}

// 2.2. FASTQC
process FASTQC_CLEANED {
    cpus 2
    conda "${envDir}/qc.yaml"
    publishDir "${pubDir}/qc/fastqc/cleaned", mode: "symlink"
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
process MULTIQC_CLEANED {
    label "single"
    conda "${envDir}/qc.yaml"
    publishDir "${pubDir}/qc/multiqc/cleaned", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("cleaned"), path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow CLEAN_READS {
    take:
        concat_ch
    main:
        clean_ch = PREPROCESS_FASTP(concat_ch, adapter_ref_path)
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
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/dedup", mode: "symlink"
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
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/dedup", mode: "symlink"
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
process MULTIQC_DEDUP {
    conda "${envDir}/qc.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/dedup", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("dedup"), path("multiqc_data")
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
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/ribo_initial", mode: "symlink"
    input:
        tuple val(sample), path(reads)
        path ribo_ref
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
        ref=!{ribo_ref}
        io_unmerged="in=${unmerged1} in2=${unmerged2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats_unmerged}"
        # Define parameters
        par="minkmerfraction=0.6 k=43 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io_unmerged} ${par}
        '''
}

// 4.2. FASTQC
process FASTQC_RIBO_INITIAL {
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/ribo_initial", mode: "symlink"
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
process MULTIQC_RIBO_INITIAL {
    conda "${envDir}/qc.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/ribo_initial", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("ribo_initial"), path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_RIBO_INITIAL {
    take:
        dedup_ch
    main:
        ribo_ch = BBDUK_RIBO_INITIAL(dedup_ch, ribo_ref_path)
        fastqc_ribo_ch = FASTQC_RIBO_INITIAL(ribo_ch)
        multiqc_ribo_ch = MULTIQC_RIBO_INITIAL(fastqc_ribo_ch.collect().ifEmpty([]))
    emit:
        data = ribo_ch
        fastqc = fastqc_ribo_ch
        multiqc_report = multiqc_ribo_ch[0]
        multiqc_data = multiqc_ribo_ch[1]
}

/*********************
| 5. HUMAN DEPLETION |
*********************/

// 5.1. Mask human genome in preparation for indexing
process BBMASK_HUMAN {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_human/mask", mode: "symlink"
    input:
        val(reference_path)
    output:
        path("human_ref.fasta.gz")
        path("human_ref_masked.fasta.gz")
    shell:
        '''
        # Get human genome reference
        in=human_ref.fasta.gz
        out=human_ref_masked.fasta.gz
        wget !{reference_path} -O ${in}
        io="in=${in} out=${out}"
        par="threads=!{task.cpus} maskrepeats masklowentropy"
        bbmask.sh ${par} ${io}
        '''
}

// 5.2. Build human index from reference genome
process BBMAP_INDEX_HUMAN {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_human/index", mode: "symlink"
    input:
        path(masked_reference)
    output:
        path("ref_index")
    shell:
        '''
        mkdir ref_index
        cp !{masked_reference} ref_index/host_ref.fasta.gz
        cd ref_index
        bbmap.sh ref=host_ref.fasta.gz t=!{task.cpus}
        '''
}

// 5.3. Segregate human reads with bbmap
process BBMAP_HUMAN_DEPLETION {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_human", mode: "symlink"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        path(index_ref_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_nohuman_{1,2}.fastq.gz"), path("${sample}_bbmap_human_{1,2}.fastq.gz"), path("${sample}_bbmap_human.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads_noribo[0]}
        unmerged2=!{reads_noribo[1]}
        op1=!{sample}_bbmap_nohuman_1.fastq.gz
        op2=!{sample}_bbmap_nohuman_2.fastq.gz
        of1=!{sample}_bbmap_human_1.fastq.gz
        of2=!{sample}_bbmap_human_2.fastq.gz
        stats=!{sample}_bbmap_human.stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} statsfile=${stats} path=!{index_ref_dir}"
        # Define parameters
        par="minid=0.9 maxindel=3 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx30g"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        '''
}

// 5.4. FASTQC
process FASTQC_HUMAN {
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/remove_human", mode: "symlink"
    input:
        tuple val(sample), path(reads_nohuman), path(reads_human), path(stats)
    output:
        path("${sample}_bbmap_nohuman_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_nohuman}
        '''
}

// 5.5. MultiQC
process MULTIQC_HUMAN {
    conda "${envDir}/qc.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/remove_human", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("remove_human"), path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_HUMAN {
    take:
        ribo_ch
    main:
        mask_ch = BBMASK_HUMAN(params.human_url)
        index_ch = BBMAP_INDEX_HUMAN(mask_ch[1])
        deplete_ch = BBMAP_HUMAN_DEPLETION(ribo_ch, index_ch)
        fastqc_deplete_ch = FASTQC_HUMAN(deplete_ch)
        multiqc_deplete_ch = MULTIQC_HUMAN(fastqc_deplete_ch.collect().ifEmpty([]))
    emit:
        data = deplete_ch
        fastqc = fastqc_deplete_ch
        multiqc_report = multiqc_deplete_ch[0]
        multiqc_data = multiqc_deplete_ch[1]
}

/*********************
| 6. OTHER DEPLETION |
*********************/

// 6.1. Join reference sequences together
// TODO: Replace hardcoded references with an extensible list
process JOIN_OTHER_REF {
    label "single"
    publishDir "${pubDir}/preprocess/remove_other/ref", mode: "symlink"
    input:
        val(cow_url)
        val(pig_url)
    output:
        path("other_ref_concat.fasta.gz")
    shell:
        '''
        # Download references
        wget !{cow_url} -O cow_ref.fasta.gz
        wget !{pig_url} -O pig_ref.fasta.gz
        in="cow_ref.fasta.gz pig_ref.fasta.gz"
        cat ${in} > other_ref_concat.fasta.gz # TODO: Make robust to differences in gzip status
        '''
}

// 6.2. Mask reference genomes in preparation for indexing
process BBMASK_REFERENCES {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_other/mask", mode: "symlink"
    input:
        path(concat_references)
    output:
        path("other_ref_masked.fasta.gz")
    shell:
        '''
        in=!{concat_references}
        out=other_ref_masked.fasta.gz
        io="in=${in} out=${out}"
        par="threads=!{task.cpus} maskrepeats masklowentropy"
        bbmask.sh ${io} ${par}
        '''

}

// 6.3. Build index from reference genomes
process BBMAP_INDEX_REFERENCES {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_other/index", mode: "symlink"
    input:
        path(masked_reference)
    output:
        path("ref_index")
    shell:
        '''
        mkdir ref_index
        cp !{masked_reference} ref_index/other_ref.fasta.gz
        cd ref_index
        bbmap.sh ref=other_ref.fasta.gz t=!{task.cpus} -Xmx30g
        '''
}

// 6.4. Segregate reference reads with bbmap
process BBMAP_REFERENCE_DEPLETION {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/remove_other", mode: "symlink"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        path(index_ref_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_noref_{1,2}.fastq.gz"), path("${sample}_bbmap_ref_{1,2}.fastq.gz"), path("${sample}_bbmap_other.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads_noribo[0]}
        unmerged2=!{reads_noribo[1]}
        op1=!{sample}_bbmap_noref_1.fastq.gz
        op2=!{sample}_bbmap_noref_2.fastq.gz
        of1=!{sample}_bbmap_ref_1.fastq.gz
        of2=!{sample}_bbmap_ref_2.fastq.gz
        stats=!{sample}_bbmap_other.stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} statsfile=${stats} path=!{index_ref_dir}"
        # Define parameters
        par="minid=0.9 maxindel=3 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx30g"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        '''
}

// 6.5. FASTQC
process FASTQC_REMOVE_OTHER {
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/remove_other", mode: "symlink"
    input:
        tuple val(sample), path(reads_noref), path(reads_ref), path(stats)
    output:
        path("${sample}_bbmap_noref_{1,2}_fastqc.{zip,html}")
    shell:
        '''
        fastqc -t !{task.cpus} !{reads_noref}
        '''
}

// 6.6. MultiQC
process MULTIQC_REMOVE_OTHER {
    conda "${envDir}/qc.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/remove_other", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("remove_other"), path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_OTHER {
    take:
        ribo_ch
    main:
        join_ch = JOIN_OTHER_REF(params.cow_url, params.pig_url) // TODO: Replace hardcoded references with extensible list
        mask_ch = BBMASK_REFERENCES(join_ch)
        index_ch = BBMAP_INDEX_REFERENCES(mask_ch)
        deplete_ch = BBMAP_REFERENCE_DEPLETION(ribo_ch, index_ch)
        fastqc_deplete_ch = FASTQC_REMOVE_OTHER(deplete_ch)
        multiqc_deplete_ch = MULTIQC_REMOVE_OTHER(fastqc_deplete_ch.collect().ifEmpty([]))
    emit:
        data = deplete_ch
        fastqc = fastqc_deplete_ch
        multiqc_report = multiqc_deplete_ch[0]
        multiqc_data = multiqc_deplete_ch[1]
}

/****************************************
| 6. HUMAN VIRUS DETECTION WITH BOWTIE2 |
****************************************/

// TODO: Write processes for Bowtie2 index construction, alignment, data parsing
// TODO: Write workflow

// 6.1. Mask collated virus genomes
// TODO: Replace with bbmask?
process MASK_HV_GENOMES {
    conda "${envDir}/bowtie.yaml"
    label "single"
    publishDir "${pubDir}/hviral/index", mode: "symlink"
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
    conda "${envDir}/bowtie.yaml"
    label "large"
    publishDir "${pubDir}/hviral/index", mode: "symlink"
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
process RUN_BOWTIE2 {
    conda "${envDir}/bowtie.yaml"
    label "large"
    publishDir "${pubDir}/hviral/bowtie", mode: "symlink"
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
    label "single"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/hviral/merged", mode: "symlink"
    input:
        tuple val(sample), path(sam_out), path(reads_out)
        path scriptDir
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
        script_path=!{scriptDir}/join_fastq.py
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
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/hviral/kraken", mode: "symlink"
    input:
        tuple val(sample), path(reads), path(stats) // Single input file, merged, joined & concatenated
        path db_path
    output:
        tuple val(sample), path("${sample}_bowtie2_kraken.output"), path("${sample}_bowtie2_kraken.report")
    shell:
        '''
        # Define input/output
        db=!{db_path}
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
    label "single"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/hviral/kraken", mode: "symlink"
    input:
        tuple val(sample), path(output), path(report)
        path scriptDir
        path nodes_path
        path hv_path
    output:
        tuple val(sample), path("${sample}_bowtie2_kraken_processed.tsv")
    shell:
        '''
        in=!{output}
        out=!{sample}_bowtie2_kraken_processed.tsv
        nodes=!{nodes_path}
        hv=!{hv_path}
        script_path=!{scriptDir}/process_kraken_hv.py
        ${script_path} ${in} ${hv} ${nodes} ${out}
        '''
}

// 6.7. Process Bowtie2 SAM output
// NB: Currently paired, need to update if switch to merged
process PROCESS_BOWTIE_SAM {
    label "single"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/hviral/bowtie", mode: "symlink"
    input:
        tuple val(sample), path(sam), path(reads)
        path scriptDir
        path genomeid_taxid_map_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv")
    shell:
        '''
        in=!{sam}
        out=!{sample}_bowtie2_sam_processed.tsv
        map=!{genomeid_taxid_map_path}
        script_path=!{scriptDir}/process_sam_hv.py
        ${script_path} ${in} ${map} ${out}
        '''
}

// 6.8. Merge processed SAM and Kraken TSVs and compute length-normalized alignment scores
process MERGE_SAM_KRAKEN {
    label "single"
    conda "${envDir}/r.yaml"
    publishDir "${pubDir}/hviral/hits", mode: "symlink"
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
    label "single"
    conda "${envDir}/r.yaml"
    publishDir "${pubDir}/hviral/hits", mode: "symlink"
    publishDir "${pubDir}/results", mode: "copy", overwrite: "true"
    input:
        path(tsvs)
    output:
        path("hv_hits_putative_all.tsv")
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
        write_tsv(tab_out, "hv_hits_putative_all.tsv")
        '''
}

// 6.10. Perform initial HV read filtering
process FILTER_HV {
    label "single"
    conda "${envDir}/r.yaml"
    publishDir "${pubDir}/hviral/hits", mode: "symlink"
    publishDir "${pubDir}/results", mode: "copy", overwrite: "true"
    input:
        path(hv_hits)
    output:
        path("hv_hits_putative_filtered.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        score_threshold <- 15 # TODO: Make a parameter
        data <- read_tsv("!{hv_hits}", col_names = TRUE, show_col_types = FALSE)
        filtered <- mutate(data, hit_hv = as.logical(!is.na(str_match(encoded_hits, paste0(" ", as.character(taxid), ":"))))) %>%
            filter((!classified) | assigned_hv) %>% 
            filter(adj_score_fwd > score_threshold | adj_score_rev > score_threshold | assigned_hv | hit_hv)
        print(dim(data))
        print(dim(filtered))
        write_tsv(filtered, "hv_hits_putative_filtered.tsv")
        '''
}

workflow MAP_HUMAN_VIRUSES {
    take:
        host_ch
    main:
        mask_ch = MASK_HV_GENOMES(hv_genomes)
        index_ch = BUILD_BOWTIE2_DB(mask_ch)
        bowtie2_ch = RUN_BOWTIE2(host_ch, index_ch)
        merge_ch = MERGE_JOIN_BOWTIE(bowtie2_ch, scriptDir)
        kraken_ch = KRAKEN_BOWTIE(merge_ch, kraken_db_path)
        kraken_processed_ch = PROCESS_KRAKEN_BOWTIE(kraken_ch, scriptDir, nodes_path, virus_db_path)
        bowtie2_processed_ch = PROCESS_BOWTIE_SAM(bowtie2_ch, scriptDir, genomeid_map_path)
        merged_input_ch = kraken_processed_ch.combine(bowtie2_processed_ch, by: 0)
        merged_ch = MERGE_SAM_KRAKEN(merged_input_ch)
        merged_ch_2 = MERGE_SAMPLES_HV(merged_ch.collect().ifEmpty([]))
        filtered_ch = FILTER_HV(merged_ch_2)
    emit:
        data_all = merged_ch_2
        data_filtered = filtered_ch
}

/*****************************
| 7. SECONDARY RIBODEPLETION |
*****************************/

// 7.1. Secondary detection and removal of ribosomal reads
// NB: Using more liberal parameters here since very high specificity less important
process BBDUK_RIBO_SECONDARY {
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/preprocess/ribo_secondary", mode: "symlink"
    input:
        tuple val(sample), path(reads_nohost), path(reads_host), path(stats)
        path ribo_ref
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
        ref=!{ribo_ref}
        io="in=${in1} in2=${in2} ref=${ref} out=${op1} out2=${op2} outm=${of1} outm2=${of2} stats=${stats}"
        # Define parameters
        par="minkmerfraction=0.4 k=27 t=!{task.cpus}"
        # Execute
        bbduk.sh ${io} ${par}
        '''
}

// 7.2. FASTQC
process FASTQC_RIBO_SECONDARY {
    conda "${envDir}/qc.yaml"
    cpus 2
    publishDir "${pubDir}/qc/fastqc/ribo_secondary", mode: "symlink"
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
process MULTIQC_RIBO_SECONDARY {
    label "single"
    conda "${envDir}/qc.yaml"
    publishDir "${pubDir}/qc/multiqc/ribo_secondary", mode: "symlink"
    input:
        path("*")
    output:
        path("multiqc_report.html")
        tuple val("ribo_secondary"), path("multiqc_data")
    shell:
        '''
        multiqc .
        '''
}

workflow REMOVE_RIBO_SECONDARY {
    take:
        novirus_ch
    main:
        ribo_ch = BBDUK_RIBO_SECONDARY(novirus_ch, ribo_ref_path)
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
    label "single"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/taxonomy/merged", mode: "symlink"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        path scriptDir
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
        script_path=!{scriptDir}/join_fastq.py
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
    label "large"
    conda "${envDir}/main.yaml"
    publishDir "${pubDir}/taxonomy/kraken", mode: "symlink"
    input:
        tuple val(sample), path(reads), path(stats)
        path db_path
    output:
        tuple val(sample), path("${sample}.output"), path("${sample}.report"), path("${sample}_unclassified.fastq.gz")
    shell:
        '''
        # Define input/output
        db=!{db_path}
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
    label "single"
    conda "${envDir}/bracken.yaml"
    publishDir "${pubDir}/taxonomy/bracken", mode: "symlink"
    input:
        tuple val(sample), path(output), path(report), path(unc_reads)
        path db_path
    output:
        tuple val(sample), path("${sample}.bracken")
    shell:
        '''
        # Define input/output
        db=!{db_path}
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
    label "single"
    conda "${envDir}/r.yaml"
    publishDir "${pubDir}/taxonomy/bracken", mode: "symlink"
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
    label "single"
    conda "${envDir}/r.yaml"
    publishDir "${pubDir}/taxonomy/results", mode: "symlink"
    publishDir "${pubDir}/results", mode: "copy", overwrite: "true"
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
        merged_ch = MERGE_JOIN(data_ch, scriptDir)
        kraken_ch = KRAKEN(merged_ch, kraken_db_path)
        bracken_ch = BRACKEN_DOMAINS(kraken_ch, kraken_db_path)
        label_ch = LABEL_BRACKEN(bracken_ch)
        merged_ch_2 = MERGE_BRACKEN(label_ch.collect().ifEmpty([]))
    emit:
        kraken = kraken_ch
        bracken = merged_ch_2
}

/*********************************
| 9. COLLATE AND PROCESS RESULTS |
*********************************/

// 9.1. Extract MultiQC data into more usable forms
process SUMMARIZE_MULTIQC_SINGLE {
    conda "${envDir}/r.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/summaries", mode: "symlink"
    input:
        tuple val(stage), path(multiqc_data)
        path(scriptDir)
    output:
        path("${stage}_qc_basic_stats.tsv")
        path("${stage}_qc_adapter_stats.tsv")
        path("${stage}_qc_quality_base_stats.tsv")
        path("${stage}_qc_quality_sequence_stats.tsv")
    shell:
        '''
        script_path=!{scriptDir}/summarize-multiqc-single.R
        ${script_path} -i !{multiqc_data} -s !{stage} -o ${PWD}
        '''
}

// 9.2. Combine MultiQC summary files across workflow stages
process MERGE_MULTIQC {
    conda "${envDir}/r.yaml"
    label "single"
    publishDir "${pubDir}/qc/multiqc/summaries", mode: "symlink"
    publishDir "${pubDir}/results", mode: "copy", overwrite: "true"
    input:
        path(basic_stats_tsvs)
        path(adapter_tsvs)
        path(base_quality_tsvs)
        path(sequence_quality_tsvs)
    output:
        path("qc_basic_stats.tsv")
        path("qc_adapter_stats.tsv")
        path("qc_quality_base_stats.tsv")
        path("qc_quality_sequence_stats.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        in_paths_basic <- str_split("!{basic_stats_tsvs}", " ")[[1]]
        in_paths_adapt <- str_split("!{adapter_tsvs}", " ")[[1]]
        in_paths_qbase <- str_split("!{base_quality_tsvs}", " ")[[1]]
        in_paths_qseqs <- str_split("!{sequence_quality_tsvs}", " ")[[1]]
        tabs_basic <- lapply(in_paths_basic, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_adapt <- lapply(in_paths_adapt, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qbase <- lapply(in_paths_qbase, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qseqs <- lapply(in_paths_qseqs, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        # Bind rows
        tab_basic_out <- bind_rows(tabs_basic)
        tab_adapt_out <- bind_rows(tabs_adapt)
        tab_qbase_out <- bind_rows(tabs_qbase)
        tab_qseqs_out <- bind_rows(tabs_qseqs)
        # Write output
        write_tsv(tab_basic_out, "qc_basic_stats.tsv")
        write_tsv(tab_adapt_out, "qc_adapter_stats.tsv")
        write_tsv(tab_qbase_out, "qc_quality_base_stats.tsv")
        write_tsv(tab_qseqs_out, "qc_quality_sequence_stats.tsv")
        '''
}

// 9.3. Summarize taxonomic composition from Bracken and MultiQC output
process SUMMARIZE_COMPOSITION {
    conda "${envDir}/r.yaml"
    label "single"
    publishDir "${pubDir}/taxonomy/results", mode: "symlink"
    publishDir "${pubDir}/results", mode: "copy", overwrite: "true"
    input:
        path(bracken_merged)
        path(multiqc_basic_merged)
    output:
        path("taxonomic_composition.tsv")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        bracken <- read_tsv("!{bracken_merged}", show_col_types = FALSE)
        basic   <- read_tsv("!{multiqc_basic_merged}", show_col_types = FALSE)
        total_assigned <- bracken %>% group_by(sample) %>% summarize(
          name = "Total",
          kraken_assigned_reads = sum(kraken_assigned_reads),
          added_reads = sum(added_reads),
          new_est_reads = sum(new_est_reads),
          fraction_total_reads = sum(fraction_total_reads)
        )
        bracken_spread <- bracken %>% select(name, sample, new_est_reads) %>%
          mutate(name = tolower(name)) %>%
          pivot_wider(id_cols = "sample", names_from = "name", values_from = "new_est_reads")
        # Count reads
        read_counts_preproc <- basic %>% select(sample, stage, n_read_pairs) %>%
          pivot_wider(id_cols = c("sample"), names_from="stage", values_from="n_read_pairs")
        read_counts <- read_counts_preproc %>%
          inner_join(total_assigned %>% select(sample, new_est_reads), by = "sample") %>%
          rename(assigned = new_est_reads) %>%
          inner_join(bracken_spread, by="sample")
        # Assess composition
        read_comp <- transmute(read_counts, sample=sample,
                               n_filtered = raw_concat-cleaned,
                               n_duplicate = cleaned-dedup,
                               n_ribosomal = (dedup-ribo_initial) + (remove_human-ribo_secondary),
                               n_unassigned = ribo_secondary-assigned,
                               n_bacterial = bacteria,
                               n_archaeal = archaea,
                               n_viral = viruses,
                               n_human = (ribo_initial-remove_human) + eukaryota,
                               n_other_filtered = remove_human-remove_other)
        read_comp_long <- pivot_longer(read_comp, -(sample), names_to = "classification",
                                       names_prefix = "n_", values_to = "n_reads") %>%
          mutate(classification = fct_inorder(str_to_sentence(classification))) %>%
          group_by(sample) %>% mutate(p_reads = n_reads/sum(n_reads))
        # Write output
        write_tsv(read_comp_long, "taxonomic_composition.tsv")
        '''
}

workflow PROCESS_OUTPUT {
    take:
        multiqc_data_raw_concat
        multiqc_data_cleaned
        multiqc_data_dedup
        multiqc_data_ribo_initial
        multiqc_data_human
        multiqc_data_other
        multiqc_data_ribo_secondary
        bracken_merged
    main:
        // Summarize each MultiQC directory separately
        multiqc_single = multiqc_data_raw_concat.mix(multiqc_data_cleaned, multiqc_data_dedup, multiqc_data_ribo_initial, multiqc_data_human, multiqc_data_other, multiqc_data_ribo_secondary)
        multiqc_summ   = SUMMARIZE_MULTIQC_SINGLE(multiqc_single, scriptDir)
        multiqc_basic = multiqc_summ[0].collect().ifEmpty([])
        multiqc_adapt = multiqc_summ[1].collect().ifEmpty([])
        multiqc_qbase = multiqc_summ[2].collect().ifEmpty([])
        multiqc_qseqs = multiqc_summ[3].collect().ifEmpty([])
        // Merge MultiQC outputs
        multiqc_merged = MERGE_MULTIQC(multiqc_basic, multiqc_adapt, multiqc_qbase, multiqc_qseqs)
        // Summarize taxonomic composition
        taxo = SUMMARIZE_COMPOSITION(bracken_merged, multiqc_merged[0])
    emit:
        basic = multiqc_merged[0]
        adapt = multiqc_merged[1]
        qbase = multiqc_merged[2]
        qseqs = multiqc_merged[3]
        composition = taxo
}

workflow {
    // Preprocessing
    HANDLE_RAW_READS(libraries_ch)
    CLEAN_READS(HANDLE_RAW_READS.out.data)
    DEDUP_READS(CLEAN_READS.out.data)
    REMOVE_RIBO_INITIAL(DEDUP_READS.out.data)
    REMOVE_HUMAN(REMOVE_RIBO_INITIAL.out.data)
    REMOVE_OTHER(REMOVE_HUMAN.out.data)
    REMOVE_RIBO_SECONDARY(REMOVE_OTHER.out.data)
    // TODO: Collate QC metrics from preprocessing
    // Human viral reads
    MAP_HUMAN_VIRUSES(REMOVE_OTHER.out.data)
    // Broad taxonomic profiling
    CLASSIFY_READS(REMOVE_RIBO_SECONDARY.out.data)
    // Process output
    PROCESS_OUTPUT(HANDLE_RAW_READS.out.multiqc_data, CLEAN_READS.out.multiqc_data, DEDUP_READS.out.multiqc_data, REMOVE_RIBO_INITIAL.out.multiqc_data, REMOVE_HUMAN.out.multiqc_data, REMOVE_OTHER.out.multiqc_data, REMOVE_RIBO_SECONDARY.out.multiqc_data, CLASSIFY_READS.out.bracken)
}
