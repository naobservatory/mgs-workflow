// Short-read version of EXTRACT_VIRAL_READS that uses streaming and interleaved files to minimize memory requirements and loading times
// This is a temporary workflow that will replace EXTRACT_VIRAL_READS_SHORT once the LCA implementation is complete
/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS_INTERLEAVE as BBDUK_HITS } from "../../../modules/local/bbduk"
include { CUTADAPT } from "../../../modules/local/cutadapt"
include { FASTP } from "../../../modules/local/fastp"
include { BOWTIE2 as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2"
include { BOWTIE2 as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2"
include { BOWTIE2 as BOWTIE2_OTHER } from "../../../modules/local/bowtie2"
include { PROCESS_VIRAL_BOWTIE2_SAM_LCA as PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2SamLca"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { EXTRACT_VIRAL_HITS_TO_FASTQ } from "../../../modules/local/extractViralHitsToFastq"
include { LCA_TSV } from "../../../modules/local/lcaTsv"
include { SORT_FASTQ } from "../../../modules/local/sortFastq"
include { SORT_FILE } from "../../../modules/local/sortFile"
include { FILTER_VIRAL_SAM } from "../../../modules/local/filterViralSam"
include { PROCESS_LCA_BOWTIE_COLUMNS } from "../processLcaBowtieColumns"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_SHORT_LCA {
    take:
        reads_ch
        ref_dir
        aln_score_threshold
        adapter_path
        cutadapt_error_rate
        min_kmer_hits
        k
        bbduk_suffix
        taxid_artificial
    main:
        // Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-masked.fasta.gz"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        nodes_db = "${ref_dir}/results/taxonomy-nodes.dmp"
        names_db = "${ref_dir}/results/taxonomy-names.dmp"
        // Define columns without duplication
        col_keep_no_prefix = ["seq_id", "aligner_taxid_lca", "aligner_taxid_top", 
                              "aligner_length_normalized_score_mean", "aligner_taxid_lca_natural",
                              "aligner_n_assignments_natural", "aligner_length_normalized_score_mean_natural",
                              "aligner_taxid_lca_artificial", "aligner_n_assignments_artificial", 
                              "aligner_length_normalized_score_mean_artificial"]
        col_keep_add_prefix = ["genome_id_all", "taxid_all", "fragment_length", 
                               "best_alignment_score", "best_alignment_score_rev",
                               "edit_distance", "edit_distance_rev", "ref_start", 
                               "ref_start_rev", "query_len", "query_len_rev",
                               "query_seq", "query_seq_rev", "query_rc", 
                               "query_rc_rev", "query_qual", "query_qual_rev", 
                               "pair_status"]
        // 1. Run initial screen against viral genomes with BBDuk
        bbduk_ch = BBDUK_HITS(reads_ch, viral_genome_path, min_kmer_hits, k, bbduk_suffix)
        // 2. Carry out stringent adapter removal with FASTP and Cutadapt
        fastp_ch = FASTP(bbduk_ch.fail, adapter_path, true)
        adapt_ch = CUTADAPT(fastp_ch.reads, adapter_path, cutadapt_error_rate)
        // 3. Run Bowtie2 against a viral database and process output
        par_virus = "--local --very-sensitive-local --score-min G,0.1,19 -k 10"
        bowtie2_ch = BOWTIE2_VIRUS(adapt_ch.reads, bt2_virus_index_path, par_virus, "virus", true, false, true)
        // 4. Filter contaminants
        par_contaminants = "--local --very-sensitive-local"
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_ch.reads_mapped, bt2_human_index_path,
            par_contaminants, "human", false, false, true)
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unmapped, bt2_other_index_path,
            par_contaminants, "other", false, false, true)
        // 5. Sort SAM and FASTQ files before filtering
        bowtie2_sam_sorted_ch = SORT_FILE(bowtie2_ch.sam, "-t\$\'\\t\' -k1,1", "sam")
        other_fastq_sorted_ch = SORT_FASTQ(other_bt2_ch.reads_unmapped)
        // 6. Consolidated viral SAM filtering: keep contaminant-free reads, applies score threshold, adds missing mates
        bowtie2_ch_combined = bowtie2_sam_sorted_ch.output.combine(other_fastq_sorted_ch.output, by: 0)
        bowtie2_filtered_ch = FILTER_VIRAL_SAM(bowtie2_ch_combined, aln_score_threshold)
        // 7. Convert SAM to TSV
        bowtie2_tsv_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_filtered_ch.sam, genome_meta_path, virus_db_path, true)
        // 8. Run LCA
        lca_ch = LCA_TSV(bowtie2_tsv_ch.output, nodes_db, names_db,
            "seq_id", "taxid", "length_normalized_score", taxid_artificial,
            "aligner")
        // 9. Process LCA and Bowtie2 columns
        processed_ch = PROCESS_LCA_BOWTIE_COLUMNS(
            lca_ch.output,
            bowtie2_tsv_ch.output,
            col_keep_no_prefix,
            col_keep_add_prefix,
            "prim_align_"
        )
        // 10. Add sample to column
        out_labeled_ch = ADD_SAMPLE_COLUMN(processed_ch.output, "sample", "viral_bowtie2")
        // 11. Concatenate across reads
        label_combined_ch = out_labeled_ch.output.map{ _sample, file -> file }.collect().ifEmpty([])
        concat_ch = CONCATENATE_TSVS(label_combined_ch, "virus_hits_final")
        // 12. Extract filtered virus hits in FASTQ format
        fastq_unfiltered_collect = other_bt2_ch.reads_unmapped.map{ _sample, file -> file }.collect().ifEmpty([])
        fastq_unfiltered_concat = CONCATENATE_FILES(fastq_unfiltered_collect, "reads_unfiltered", "fastq.gz")
        fastq_ch = EXTRACT_VIRAL_HITS_TO_FASTQ(concat_ch.output, fastq_unfiltered_concat.output)
    emit:
        bbduk_match = bbduk_ch.fail
        bbduk_trimmed = adapt_ch.reads
        hits_final = concat_ch.output
        hits_prelca = bowtie2_tsv_ch.output
        hits_fastq = fastq_ch.fastq
        test_reads = other_bt2_ch.reads_unmapped
        test_filt_bowtie = bowtie2_filtered_ch.sam
        test_unfilt_bowtie = bowtie2_ch.sam
}
