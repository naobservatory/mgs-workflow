// Version of EXTRACT_VIRAL_READS that uses streaming and interleaved files to minimize memory requirements and loading times

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS_STREAMED } from "../../../modules/local/bbduk"
include { CUTADAPT_STREAMED } from "../../../modules/local/cutadapt"
include { BOWTIE2_STREAMED as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_OTHER } from "../../../modules/local/bowtie2"
include { BBMAP_STREAMED as BBMAP_HUMAN } from "../../../modules/local/bbmap"
include { BBMAP_STREAMED as BBMAP_OTHER } from "../../../modules/local/bbmap"
include { TAXONOMY_STREAMED as TAXONOMY } from "../../../subworkflows/local/taxonomyStreamed"
include { PROCESS_VIRAL_BOWTIE2_SAM_2 as PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2Sam" // NB: Already streamed
include { PROCESS_KRAKEN_VIRAL_2 as PROCESS_KRAKEN_VIRAL } from "../../../modules/local/processKrakenViral" // NB: Already streamed
include { SORT_TSV as SORT_KRAKEN_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { REHEAD_TSV as REHEAD_BOWTIE_VIRAL } from "../../../modules/local/reheadTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"

include { CONCATENATE_TSVS as CONCATENATE_TSVS_BOWTIE2_KRAKEN } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_TSVS_BBMERGE } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_TSVS_DEDUP } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_TSVS_ALIGNMENT_DUPLICATES } from "../../../modules/local/concatenateTsvs"
include { FILTER_VIRUS_READS } from "../../../modules/local/filterVirusReads"
include { COLLAPSE_VIRUS_READS } from "../../../modules/local/collapseVirusReads"
include { ADD_FRAG_DUP_TO_VIRUS_READS } from "../../../modules/local/addFragDupToVirusReads"
include { MAKE_VIRUS_READS_FASTA } from "../../../modules/local/makeVirusReadsFasta"
include { COUNT_VIRUS_CLADES } from "../../../modules/local/countVirusClades"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_STREAMED {
    take:
        reads_ch
        group_ch
        ref_dir
        kraken_db_ch
        aln_score_threshold
        adapter_path
        host_taxon
        min_kmer_hits
        k
        bbduk_suffix
        encoding
        fuzzy_match
        grouping
        single_end
    main:
        // 0. Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-filtered.fasta.gz"
        genome_meta_path  = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        bbm_human_index_path = "${ref_dir}/results/bbm-human-index"
        bbm_other_index_path = "${ref_dir}/results/bbm-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        // 1. Run initial screen against viral genomes with BBDuk
        bbduk_ch = BBDUK_HITS_STREAMED(reads_ch, viral_genome_path, min_kmer_hits, k, bbduk_suffix)
        // 2. Carry out stringent adapter removal with Cutadapt and Trimmomatic
        adapt_ch = CUTADAPT_STREAMED(bbduk_ch.fail, adapter_path)
        // NB: Dropping Trimmomatic here for now; note that this will produce additional false positives until this is replaced (in active development in another branch)
        // TODO: Replace Trimmomatic with a different tool (once Harmon is done investigating)
        trim_ch = adapt_ch
        // NB: No grouping, all readwise (i.e. no dedup)
        // 3. Run Bowtie2 against a viral database and process output
        bowtie2_ch = BOWTIE2_VIRUS(trim_ch.reads, bt2_virus_index_path, "--score-min G,1,1", "virus", true, false)
        // TODO: Process Bowtie2 SAM output for integration into output TSV (PROCESS_VIRAL_BOWTIE2_SAM)
        // 4. Filter contaminants
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_ch.reads_mapped, bt2_human_index_path, "", "human", false, false)
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unmapped, bt2_other_index_path, "", "other", false, false)
        human_bbm_ch = BBMAP_HUMAN(other_bt2_ch.reads_unmapped, bbm_human_index_path, "human", false, false)
        other_bbm_ch = BBMAP_OTHER(human_bbm_ch.reads_unmapped, bbm_other_index_path, "other", false, false)
        // 5. Run Kraken on filtered viral candidates (via taxonomy subworkflow)
        tax_ch = TAXONOMY(other_bbm_ch.reads_unmapped, kraken_db_ch, "F", single_end)
        // 6. Process and combine Kraken and Bowtie2 output
        bowtie2_sam_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_ch.sam, genome_meta_path, virus_db_path)
        kraken_output_ch = PROCESS_KRAKEN_VIRAL(tax_ch.kraken_output, virus_db_path, host_taxon)
        bowtie2_sam_rehead_ch = REHEAD_BOWTIE_VIRAL(bowtie2_sam_ch.output, "query_name", "seq_id", "bowtie2_viral")
        bowtie2_sam_sorted_ch = SORT_BOWTIE_VIRAL(bowtie2_sam_rehead_ch.output, "seq_id", "bowtie2_viral")
        kraken_sorted_ch = SORT_KRAKEN_VIRAL(kraken_output_ch.output, "seq_id", "kraken_viral")
        out_combined_ch = bowtie2_sam_sorted_ch.sorted.combine(kraken_sorted_ch.sorted, by: 0)
        out_joined_ch = JOIN_TSVS(out_combined_ch, "seq_id", "inner", "bowtie2_kraken_viral")
        out_labeled_ch = ADD_SAMPLE_COLUMN(out_joined_ch, "sample", "bowtie2_kraken_viral")
    emit:
        bbduk_match = bbduk_ch.fail
        reads_test  = other_bbm_ch.reads_unmapped
        kraken_test = tax_ch.kraken_output
}

//workflow EXTRACT_VIRAL_READS {
//    main:
//        trim_ch = TRIMMOMATIC(adapt_ch.reads, adapter_path, encoding)
//        // Process Kraken output and merge with Bowtie2 output across samples
//        bowtie2_kraken_merged_ch = MERGE_SAM_KRAKEN(kraken_output_ch.combine(bowtie2_sam_ch, by: 0))
//        merged_ch = CONCATENATE_TSVS_BOWTIE2_KRAKEN(bowtie2_kraken_merged_ch.collect().ifEmpty([]), "bowtie2_kraken_merged")
//        merged_bbmerge_results = CONCATENATE_TSVS_BBMERGE(tax_ch.bbmerge_summary.collect().ifEmpty([]), "bbmerge")
//        merged_dedup_results = CONCATENATE_TSVS_DEDUP(tax_ch.dedup_summary.collect().ifEmpty([]), "dedup")
//        merged_alignment_dup_results = CONCATENATE_TSVS_ALIGNMENT_DUPLICATES(alignment_dup_summary.collect().ifEmpty([]), "alignment_duplicates")
//        // Filter and process putative hit TSV
//        filtered_ch = FILTER_VIRUS_READS(merged_ch, aln_score_threshold)
//        collapsed_ch = COLLAPSE_VIRUS_READS(filtered_ch)
//        collapsed_frag_dup_ch = ADD_FRAG_DUP_TO_VIRUS_READS(collapsed_ch, merged_bbmerge_results, merged_dedup_results, merged_alignment_dup_results)
//        fasta_ch = MAKE_VIRUS_READS_FASTA(collapsed_frag_dup_ch)
//        // Count clades
//        count_ch = COUNT_VIRUS_CLADES(collapsed_frag_dup_ch, virus_db_path)
//    emit:
//        tsv = collapsed_frag_dup_ch
//        fasta = fasta_ch
//        counts = count_ch
//        bbduk_match = bbduk_ch.fail
//}
