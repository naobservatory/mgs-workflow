/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS } from "../../../modules/local/bbduk" addParams(suffix: params.bbduk_suffix)
include { CUTADAPT } from "../../../modules/local/cutadapt"
include { TRIMMOMATIC } from "../../../modules/local/trimmomatic" addParams(encoding: params.encoding)
include { BOWTIE2 as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2" addParams(suffix: "virus")
include { BOWTIE2 as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2" addParams(suffix: "human")
include { BOWTIE2 as BOWTIE2_OTHER } from "../../../modules/local/bowtie2" addParams(suffix: "other")
include { PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2Sam"
include { COUNT_ALIGNMENT_DUPLICATES } from "../../../modules/local/countAlignmentDuplicates" addParams(fuzzy_match: params.fuzzy_match.toInteger())
include { EXTRACT_UNCONC_READ_IDS } from "../../../modules/local/extractUnconcReadIDs"
include { EXTRACT_UNCONC_READS } from "../../../modules/local/extractUnconcReads"
include { COMBINE_MAPPED_BOWTIE2_READS } from "../../../modules/local/combineMappedBowtie2Reads"
include { BBMAP as BBMAP_HUMAN } from "../../../modules/local/bbmap" addParams(suffix: "human")
include { BBMAP as BBMAP_OTHER } from "../../../modules/local/bbmap" addParams(suffix: "other")
include { TAXONOMY } from "../../../subworkflows/local/taxonomy" addParams(dedup_rc: true, classification_level: "F", read_fraction: 1)
include { PROCESS_KRAKEN_VIRAL } from "../../../modules/local/processKrakenViral"
include { MERGE_SAM_KRAKEN } from "../../../modules/local/mergeSamKraken"
include { MERGE_TSVS as MERGE_TSVS_BOWTIE2_KRAKEN } from "../../../modules/local/mergeTsvs" addParams(name: "bowtie2_kraken_merged")
include { MERGE_TSVS as MERGE_TSVS_BBMERGE } from "../../../modules/local/mergeTsvs" addParams(name: "bbmerge")
include { MERGE_TSVS as MERGE_TSVS_DEDUP } from "../../../modules/local/mergeTsvs" addParams(name: "dedup")
include { MERGE_TSVS as MERGE_TSVS_ALIGNMENT_DUPLICATES } from "../../../modules/local/mergeTsvs" addParams(name: "alignment_duplicates")
include { FILTER_VIRUS_READS } from "../../../modules/local/filterVirusReads"
include { COLLAPSE_VIRUS_READS } from "../../../modules/local/collapseVirusReads"
include { ADD_FRAG_DUP_TO_VIRUS_READS } from "../../../modules/local/addFragDupToVirusReads"
include { MAKE_VIRUS_READS_FASTA } from "../../../modules/local/makeVirusReadsFasta"
include { COUNT_VIRUS_CLADES } from "../../../modules/local/countVirusClades"
include { CONCAT_GROUP } from "../../../modules/local/concatGroup"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS {
    take:
        reads_ch
        group_ch
        ref_dir
        kraken_db_ch
        aln_score_threshold
        adapter_path
        host_taxon
    main:
        // Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-filtered.fasta.gz"
        genome_meta_path  = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        bbm_human_index_path = "${ref_dir}/results/bbm-human-index"
        bbm_other_index_path = "${ref_dir}/results/bbm-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        // Run initial screen against viral genomes with BBDuk
        bbduk_ch = BBDUK_HITS(reads_ch, viral_genome_path, params.min_kmer_hits, params.k)
        // Carry out stringent adapter removal with Cutadapt and Trimmomatic
        adapt_ch = CUTADAPT(bbduk_ch.fail, adapter_path)
        trim_ch = TRIMMOMATIC(adapt_ch.reads, adapter_path)
        // Grouping for deduplication
        if (params.grouping) {
            // Join samplesheet with trimmed_reads and update fastq files
            trim_group_ch = group_ch.join(trim_ch.reads, by: 0)
            .map { sample, group, reads -> tuple(sample, reads[0], reads[1], group) }
            .groupTuple(by: 3)
            // Split into multi-sample and single-sample groups
            multi_sample_groups = trim_group_ch.filter { it[0].size() > 1 }
            single_sample_groups = trim_group_ch.filter { it[0].size() == 1 }
                .map { samples, fwd_list, rev_list, group -> tuple(group, [fwd_list[0], rev_list[0]]) }
            grouped_ch = CONCAT_GROUP(multi_sample_groups).mix(single_sample_groups)
        } else {
            grouped_ch = trim_ch.reads
        }
        // Run Bowtie2 against a viral database and process output
        bowtie2_ch = BOWTIE2_VIRUS(grouped_ch, bt2_virus_index_path, "--no-unal --no-sq --score-min G,1,1")
        bowtie2_sam_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_ch.sam, genome_meta_path, virus_db_path)
        alignment_dup_summary = COUNT_ALIGNMENT_DUPLICATES(bowtie2_sam_ch)
        bowtie2_unconc_ids_ch = EXTRACT_UNCONC_READ_IDS(bowtie2_ch.sam)
        bowtie2_unconc_reads_ch = EXTRACT_UNCONC_READS(bowtie2_ch.reads_unconc.combine(bowtie2_unconc_ids_ch, by: 0))
        bowtie2_reads_combined_ch = COMBINE_MAPPED_BOWTIE2_READS(bowtie2_ch.reads_conc.combine(bowtie2_unconc_reads_ch, by: 0))
        // Filter contaminants
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_reads_combined_ch, bt2_human_index_path, "")
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unconc, bt2_other_index_path, "")
        human_bbm_ch = BBMAP_HUMAN(other_bt2_ch.reads_unconc, bbm_human_index_path)
        other_bbm_ch = BBMAP_OTHER(human_bbm_ch.reads_unmapped, bbm_other_index_path)
        // Run Kraken on filtered viral candidates
        tax_ch = TAXONOMY(other_bbm_ch.reads_unmapped, kraken_db_ch)
        // Process Kraken output and merge with Bowtie2 output across samples
        kraken_output_ch = PROCESS_KRAKEN_VIRAL(tax_ch.kraken_output, virus_db_path, host_taxon)
        bowtie2_kraken_merged_ch = MERGE_SAM_KRAKEN(kraken_output_ch.combine(bowtie2_sam_ch, by: 0))
        merged_ch = MERGE_TSVS_BOWTIE2_KRAKEN(bowtie2_kraken_merged_ch.collect().ifEmpty([]))
        merged_bbmerge_results = MERGE_TSVS_BBMERGE(tax_ch.bbmerge_summary.collect().ifEmpty([]))
        merged_dedup_results = MERGE_TSVS_DEDUP(tax_ch.dedup_summary.collect().ifEmpty([]))
        merged_alignment_dup_results = MERGE_TSVS_ALIGNMENT_DUPLICATES(alignment_dup_summary.collect().ifEmpty([]))
        // Filter and process putative hit TSV
        filtered_ch = FILTER_VIRUS_READS(merged_ch, aln_score_threshold)
        collapsed_ch = COLLAPSE_VIRUS_READS(filtered_ch)
        collapsed_frag_dup_ch = ADD_FRAG_DUP_TO_VIRUS_READS(collapsed_ch, merged_bbmerge_results, merged_dedup_results, merged_alignment_dup_results)
        fasta_ch = MAKE_VIRUS_READS_FASTA(collapsed_frag_dup_ch)
        // Count clades
        count_ch = COUNT_VIRUS_CLADES(collapsed_frag_dup_ch, virus_db_path)
    emit:
        tsv = collapsed_frag_dup_ch
        fasta = fasta_ch
        counts = count_ch
}
