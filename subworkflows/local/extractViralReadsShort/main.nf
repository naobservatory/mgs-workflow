// Short-read version of EXTRACT_VIRAL_READS that uses streaming and interleaved files to minimize memory requirements and loading times

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS_INTERLEAVE as BBDUK_HITS } from "../../../modules/local/bbduk"
include { CUTADAPT } from "../../../modules/local/cutadapt"
include { FASTP } from "../../../modules/local/fastp"
include { BOWTIE2 as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2"
include { BOWTIE2 as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2"
include { BOWTIE2 as BOWTIE2_OTHER } from "../../../modules/local/bowtie2"
include { TAXONOMY } from "../../../subworkflows/local/taxonomy"
include { PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2Sam"
include { PROCESS_KRAKEN_VIRAL } from "../../../modules/local/processKrakenViral"
include { SORT_TSV as SORT_KRAKEN_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BBMERGE } from "../../../modules/local/sortTsv"
include { REHEAD_TSV as REHEAD_BOWTIE_VIRAL } from "../../../modules/local/reheadTsv"
include { JOIN_TSVS as JOIN_KRAKEN_BOWTIE } from "../../../modules/local/joinTsvs"
include { JOIN_TSVS as JOIN_BBMERGE } from "../../../modules/local/joinTsvs"
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"
include { FILTER_VIRUS_READS } from "../../../modules/local/filterVirusReads"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { EXTRACT_VIRAL_HITS_TO_FASTQ } from "../../../modules/local/extractViralHitsToFastq"


include { LCA_TSV } from "../../../modules/local/lcaTsv"
include { FILTER_VIRAL_SAM } from "../../../modules/local/filterViralSam"
include { PROCESS_VIRAL_BOWTIE2_SAM as SAM_TO_TSV } from "../../../modules/local/processViralBowtie2Sam"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_SHORT {
    take:
        reads_ch
        ref_dir
        kraken_db_ch
        aln_score_threshold
        adapter_path
        host_taxon
        cutadapt_error_rate
        min_kmer_hits
        k
        bbduk_suffix
        bracken_threshold
    main:
        // Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-masked.fasta.gz"
        genome_meta_path  = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        nodes_db = "${ref_dir}/results/taxonomy-nodes.dmp"
        names_db = "${ref_dir}/results/taxonomy-names.dmp"

        // 1. Run initial screen against viral genomes with BBDuk
        bbduk_ch = BBDUK_HITS(reads_ch, viral_genome_path, min_kmer_hits, k, bbduk_suffix)
        // 2. Carry out stringent adapter removal with FASTP and Cutadapt
        fastp_ch = FASTP(bbduk_ch.fail, adapter_path, true)
        adapt_ch = CUTADAPT(fastp_ch.reads, adapter_path, cutadapt_error_rate)
        // 3. Run Bowtie2 against a viral database and process output
        par_virus = "--local --very-sensitive-local --score-min G,0.1,19 -k 10"
        bowtie2_ch = BOWTIE2_VIRUS(adapt_ch.reads, bt2_virus_index_path,
            par_virus, "virus", true, false)
        // 4. Filter contaminants
        par_contaminants = "--local --very-sensitive-local"
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_ch.reads_mapped, bt2_human_index_path,
            par_contaminants, "human", false, false)
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unmapped, bt2_other_index_path,
            par_contaminants, "other", false, false)

        // Consolidated viral SAM filtering: removes contaminants, applies score threshold, adds missing mates
        bowtie2_ch_combined = bowtie2_ch.sam.combine(other_bt2_ch.reads_unmapped, by: 0)
        bowtie2_filtered_ch = FILTER_VIRAL_SAM(bowtie2_ch_combined, aln_score_threshold)

        // Convert SAM to TSV
        bowtie2_tsv_ch = SAM_TO_TSV(bowtie2_filtered_ch.sam, genome_meta_path, virus_db_path)

        // Run LCA
        lca_ch = LCA_TSV(bowtie2_tsv_ch.output, nodes_db, names_db,
            "query_name", "bowtie2_taxid_best", "bowtie2_length_normalized_score_max", 1,
            "bowtie2")

        // 5. Run Kraken on filtered viral candidates (via taxonomy subworkflow)
//      tax_ch = TAXONOMY(other_bt2_ch.reads_unmapped, kraken_db_ch, "F", bracken_threshold, channel.value(false))
        // 6. Process and combine Kraken and Bowtie2 output
//      bowtie2_sam_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_ch.sam, genome_meta_path, virus_db_path)
//      kraken_output_ch = PROCESS_KRAKEN_VIRAL(tax_ch.kraken_output, virus_db_path, host_taxon)
//      bowtie2_sam_rehead_ch = REHEAD_BOWTIE_VIRAL(bowtie2_sam_ch.output, "query_name", "seq_id")
//      bowtie2_sam_sorted_ch = SORT_BOWTIE_VIRAL(bowtie2_sam_rehead_ch.output, "seq_id")
//      kraken_sorted_ch = SORT_KRAKEN_VIRAL(kraken_output_ch.output, "seq_id")
//      out_combined_ch = bowtie2_sam_sorted_ch.sorted.combine(kraken_sorted_ch.sorted, by: 0)
//      out_joined_ch = JOIN_KRAKEN_BOWTIE(out_combined_ch, "seq_id", "inner", "bowtie2_kraken_viral")
//      // 7. Add BBMerge summary information
//      bbmerge_sorted_ch = SORT_BBMERGE(tax_ch.bbmerge_summary, "seq_id")
//      bbmerge_combined_ch = out_joined_ch.output.combine(bbmerge_sorted_ch.sorted, by: 0)
//      out_joined_ch_2 = JOIN_BBMERGE(bbmerge_combined_ch, "seq_id", "left", "viral_bbmerge")
//      out_labeled_ch = ADD_SAMPLE_COLUMN(out_joined_ch_2.output, "sample", "viral_bbmerge")
//      // 8. Concatenate across reads
//      label_combined_ch = out_labeled_ch.output.map{ sample, file -> file }.collect().ifEmpty([])
//      concat_ch = CONCATENATE_TSVS(label_combined_ch, "virus_hits_unfiltered")
//      // 9. Filter by length-normalized alignment score
//      filter_ch = FILTER_VIRUS_READS(concat_ch.output, aln_score_threshold, "virus_hits_final")
//      // 10. Extract filtered virus hits in FASTQ format
//      fastq_unfiltered_collect = other_bt2_ch.reads_unmapped.map{ sample, file -> file }.collect().ifEmpty([])
//      fastq_unfiltered_concat = CONCATENATE_FILES(fastq_unfiltered_collect, "reads_unfiltered", "fastq.gz")
//      fastq_ch = EXTRACT_VIRAL_HITS_TO_FASTQ(filter_ch.output, fastq_unfiltered_concat.output)
    emit:
        bbduk_match = bbduk_ch.fail
        bbduk_trimmed = adapt_ch.reads
//      hits_unfiltered = concat_ch.output
//      hits_final = filter_ch.output
//      hits_fastq = fastq_ch.fastq
//      test_reads  = other_bt2_ch.reads_unmapped
//      test_kraken = kraken_output_ch.output
//      test_bowtie = bowtie2_sam_ch.output
//      test_joined = out_labeled_ch.output
}
