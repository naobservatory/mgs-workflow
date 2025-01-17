// Version of EXTRACT_VIRAL_READS that uses streaming and interleaved files to minimize memory requirements and loading times

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS_STREAMED } from "../../../modules/local/bbduk"
include { CUTADAPT_STREAMED as CUTADAPT } from "../../../modules/local/cutadapt"
include { FASTP_PAIRED_STREAMED as FASTP } from "../../../modules/local/fastp"
include { BOWTIE2_STREAMED as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_OTHER } from "../../../modules/local/bowtie2"
include { TAXONOMY_STREAMED as TAXONOMY } from "../../../subworkflows/local/taxonomyStreamed"
include { PROCESS_VIRAL_BOWTIE2_SAM_2 as PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2Sam" // NB: Already streamed
include { PROCESS_KRAKEN_VIRAL_2 as PROCESS_KRAKEN_VIRAL } from "../../../modules/local/processKrakenViral" // NB: Already streamed
include { SORT_TSV as SORT_KRAKEN_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BBMERGE } from "../../../modules/local/sortTsv"
include { REHEAD_TSV as REHEAD_BOWTIE_VIRAL } from "../../../modules/local/reheadTsv"
include { JOIN_TSVS as JOIN_KRAKEN_BOWTIE } from "../../../modules/local/joinTsvs"
include { JOIN_TSVS as JOIN_BBMERGE } from "../../../modules/local/joinTsvs"
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS_STREAMED as CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"
include { FILTER_VIRUS_READS_STREAMED as FILTER_VIRUS_READS } from "../../../modules/local/filterVirusReads"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { EXTRACT_VIRAL_HITS_TO_FASTQ } from "../../../modules/local/extractViralHitsToFASTQ"

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
        // TODO: Check and remove unused inputs (e.g. grouping)
    main:
        // 0. Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-filtered.fasta.gz"
        genome_meta_path  = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        // 1. Run initial screen against viral genomes with BBDuk
        bbduk_ch = BBDUK_HITS_STREAMED(reads_ch, viral_genome_path, min_kmer_hits, k, bbduk_suffix)
        // 2. Carry out stringent adapter removal with FASTP and Cutadapt
        fastp_ch = FASTP(bbduk_ch.fail, adapter_path)
        adapt_ch = CUTADAPT(fastp_ch.reads, adapter_path)
        // NB: No grouping, all readwise (i.e. no dedup)
        // 3. Run Bowtie2 against a viral database and process output
        bowtie2_ch = BOWTIE2_VIRUS(adapt_ch.reads, bt2_virus_index_path, "--score-min G,1,1", "virus", true, false)
        // 4. Filter contaminants
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_ch.reads_mapped, bt2_human_index_path, "", "human", false, false)
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unmapped, bt2_other_index_path, "", "other", false, false)
        // 5. Run Kraken on filtered viral candidates (via taxonomy subworkflow)
        tax_ch = TAXONOMY(other_bt2_ch.reads_unmapped, kraken_db_ch, "F", single_end)
        // 6. Process and combine Kraken and Bowtie2 output
        bowtie2_sam_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_ch.sam, genome_meta_path, virus_db_path)
        kraken_output_ch = PROCESS_KRAKEN_VIRAL(tax_ch.kraken_output, virus_db_path, host_taxon)
        bowtie2_sam_rehead_ch = REHEAD_BOWTIE_VIRAL(bowtie2_sam_ch.output, "query_name", "seq_id", "bowtie2_viral")
        bowtie2_sam_sorted_ch = SORT_BOWTIE_VIRAL(bowtie2_sam_rehead_ch.output, "seq_id", "bowtie2_viral")
        kraken_sorted_ch = SORT_KRAKEN_VIRAL(kraken_output_ch.output, "seq_id", "kraken_viral")
        out_combined_ch = bowtie2_sam_sorted_ch.sorted.combine(kraken_sorted_ch.sorted, by: 0)
        out_joined_ch = JOIN_KRAKEN_BOWTIE(out_combined_ch, "seq_id", "inner", "bowtie2_kraken_viral")
        // 7. Add BBMerge summary information
        bbmerge_sorted_ch = SORT_BBMERGE(tax_ch.bbmerge_summary, "seq_id", "bbmerge_summary")
        bbmerge_combined_ch = out_joined_ch.output.combine(bbmerge_sorted_ch.sorted, by: 0)
        out_joined_ch_2 = JOIN_BBMERGE(bbmerge_combined_ch, "seq_id", "left", "viral_bbmerge")
        out_labeled_ch = ADD_SAMPLE_COLUMN(out_joined_ch_2.output, "sample", "viral_bbmerge")
        // 8. Concatenate across reads
        label_combined_ch = out_labeled_ch.output.map{ sample, file -> file }.collect().ifEmpty([])
        concat_ch = CONCATENATE_TSVS(label_combined_ch, "virus_hits_all")
        // 9. Filter by length-normalized alignment score
        filter_ch = FILTER_VIRUS_READS(concat_ch.output, aln_score_threshold)
        // 10. Extract filtered virus hits in FASTQ format
        fastq_unfiltered_collect = other_bt2_ch.reads_unmapped.map{ sample, file -> file }.collect().ifEmpty([])
        fastq_unfiltered_concat = CONCATENATE_FILES(fastq_unfiltered_collect, "reads_unfiltered", "fastq.gz")
        fastq_ch = EXTRACT_VIRAL_HITS_TO_FASTQ(filter_ch.output, fastq_unfiltered_concat.output)
    emit:
        bbduk_match = bbduk_ch.fail
        hits_all = concat_ch.output
        hits_filtered = filter_ch.output
        hits_fastq = fastq_ch.fastq
        test_reads  = other_bt2_ch.reads_unmapped
        test_kraken = kraken_output_ch.output
        test_bowtie = bowtie2_sam_ch.output
        test_joined = out_labeled_ch.output
}

// Removed functionality, to be moved to a new script or workflow
// - Grouping for deduplication
// - Deduplication with Clumpify
// - Duplicate annotation with Bowtie2
// - Addition of duplicate information to output TSV
// - Clade counting
