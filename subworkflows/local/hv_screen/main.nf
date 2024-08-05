/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BOWTIE2 as BOWTIE2_HV } from "../../../modules/local/bowtie2" addParams(suffix: "hv")
include { EXTRACT_UNCONC_READ_IDS } from "../../../modules/local/extractUnconcReadIDs"
include { EXTRACT_UNCONC_READS } from "../../../modules/local/extractUnconcReads"
include { COMBINE_MAPPED_BOWTIE2_READS } from "../../../modules/local/combineMappedBowtie2Reads"
include { FASTP } from "../../../modules/local/fastp"
include { BBDUK } from "../../../modules/local/bbduk" addParams(suffix: params.bbduk_suffix)
include { PROCESS_BOWTIE2_SAM_PAIRED } from "../../../modules/local/processBowtie2Sam"

/***********
| WORKFLOW |
***********/

workflow HV_SCREEN_1 {
    // Directly run Bowtie2 on raw reads
    take:
        reads_ch
        ref_dir
    main:
        // Get reference paths
        bt2_hv_index_path = "${ref_dir}/results/bt2-hv-index"
        genomeid_map_path = "${ref_dir}/results/genomeid-to-taxid.json"
        // Run Bowtie2 against an HV database and process output
        bowtie2_ch = BOWTIE2_HV(reads_ch, bt2_hv_index_path, "--no-unal --no-sq --score-min G,5,11")
        bowtie2_unconc_ids_ch = EXTRACT_UNCONC_READ_IDS(bowtie2_ch.sam)
        bowtie2_unconc_reads_ch = EXTRACT_UNCONC_READS(bowtie2_ch.reads_unconc.combine(bowtie2_unconc_ids_ch, by: 0))
        bowtie2_reads_combined_ch = COMBINE_MAPPED_BOWTIE2_READS(bowtie2_ch.reads_conc.combine(bowtie2_unconc_reads_ch, by: 0))
        bowtie2_sam_ch = PROCESS_BOWTIE2_SAM_PAIRED(bowtie2_ch.sam, genomeid_map_path)
        // Concatenate across samples
        file1_ch = bowtie2_reads_combined_ch.map{ it[1][0] }.collectFile(name: "screen_1_1.fastq.gz")
        file2_ch = bowtie2_reads_combined_ch.map{ it[1][1] }.collectFile(name: "screen_1_2.fastq.gz")
        out_ch = file1_ch.mix(file2_ch)
    emit:
        fastq = out_ch
        sam = bowtie2_sam_ch
}

workflow HV_SCREEN_2 {
    // Run FASTP on raw reads, then Bowtie2 on FASTP output
    take:
        reads_ch
        ref_dir
    main:
        // Get reference paths
        bt2_hv_index_path = "${ref_dir}/results/bt2-hv-index"
        genomeid_map_path = "${ref_dir}/results/genomeid-to-taxid.json"
        // Run FASTP on raw reads
        fastp_ch = FASTP(reads_ch, params.adapter_path)
        // Run Bowtie2 against an HV database and process output
        bowtie2_ch = BOWTIE2_HV(fastp_ch.reads, bt2_hv_index_path, "--no-unal --no-sq --score-min G,5,11")
        bowtie2_unconc_ids_ch = EXTRACT_UNCONC_READ_IDS(bowtie2_ch.sam)
        bowtie2_unconc_reads_ch = EXTRACT_UNCONC_READS(bowtie2_ch.reads_unconc.combine(bowtie2_unconc_ids_ch, by: 0))
        bowtie2_reads_combined_ch = COMBINE_MAPPED_BOWTIE2_READS(bowtie2_ch.reads_conc.combine(bowtie2_unconc_reads_ch, by: 0))
        bowtie2_sam_ch = PROCESS_BOWTIE2_SAM_PAIRED(bowtie2_ch.sam, genomeid_map_path)
        // Concatenate across samples
        file1_ch = bowtie2_reads_combined_ch.map{ it[1][0] }.collectFile(name: "screen_2_1.fastq.gz")
        file2_ch = bowtie2_reads_combined_ch.map{ it[1][1] }.collectFile(name: "screen_2_2.fastq.gz")
        out_ch = file1_ch.mix(file2_ch)
    emit:
        fastq = out_ch
        sam = bowtie2_sam_ch
}

workflow HV_SCREEN_3 {
    // Run BBDuk without Bowtie2 or FASTP
    take:
        reads_ch
        ref_dir
    main:
        // Get reference paths
        hv_ref_path = "${ref_dir}/results/human-viral-genomes-filtered.fasta.gz"
        // Run BBDUK
        bbduk_ch = BBDUK(reads_ch, hv_ref_path, params.min_kmer_fraction, params.k)
        // Concatenate across samples
        reads_ch = bbduk_ch.fail
        file1_ch = reads_ch.map{ it[1][0] }.collectFile(name: "screen_3_${params.min_kmer_fraction}_${params.k}_1.fastq.gz")
        file2_ch = reads_ch.map{ it[1][1] }.collectFile(name: "screen_3_${params.min_kmer_fraction}_${params.k}_2.fastq.gz")
        out_ch = file1_ch.mix(file2_ch)
        stats_ch = bbduk_ch.log.map{ it[1] }.collectFile(name: "bbduk_raw_stats_collected.txt")
    emit:
        fastq = out_ch
        stats = stats_ch
}

workflow HV_SCREEN_4 {
    // Run BBDuk after FASTP (but without Bowtie2)
    take:
        reads_ch
        ref_dir
    main:
        // Get reference paths
        hv_ref_path = "${ref_dir}/results/human-viral-genomes-filtered.fasta.gz"
        // Run FASTP on raw reads
        fastp_ch = FASTP(reads_ch, params.adapter_path)
        // Run BBDUK
        bbduk_ch = BBDUK(fastp_ch.reads, hv_ref_path, params.min_kmer_fraction, params.k)
        // Concatenate across samples
        reads_ch = bbduk_ch.fail
        file1_ch = reads_ch.map{ it[1][0] }.collectFile(name: "screen_4_${params.min_kmer_fraction}_${params.k}_1.fastq.gz")
        file2_ch = reads_ch.map{ it[1][1] }.collectFile(name: "screen_4_${params.min_kmer_fraction}_${params.k}_2.fastq.gz")
        out_ch = file1_ch.mix(file2_ch)
        stats_ch = bbduk_ch.log.map{ it[1] }.collectFile(name: "bbduk_cleaned_stats_collected.txt")
    emit:
        fastq = out_ch
        stats = stats_ch
}

workflow HV_SCREEN_5 {
    // Run BBDuk using a masked HV reference
    take:
        reads_ch
        masked_db
    main:
        // Run BBDUK
        bbduk_ch = BBDUK(reads_ch, masked_db, params.min_kmer_fraction, params.k)
        // Concatenate across samples
        reads_ch = bbduk_ch.fail
        file1_ch = reads_ch.map{ it[1][0] }.collectFile(name: "screen_5_${params.min_kmer_fraction}_${params.k}_1.fastq.gz")
        file2_ch = reads_ch.map{ it[1][1] }.collectFile(name: "screen_5_${params.min_kmer_fraction}_${params.k}_2.fastq.gz")
        out_ch = file1_ch.mix(file2_ch)
        stats_ch = bbduk_ch.log.map{ it[1] }.collectFile(name: "bbduk_masked_stats_collected.txt")
    emit:
        fastq = out_ch
        stats = stats_ch
}
