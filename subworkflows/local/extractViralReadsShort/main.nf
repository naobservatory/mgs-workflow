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
include { PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2Sam"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { EXTRACT_VIRAL_HITS_TO_FASTQ } from "../../../modules/local/extractViralHitsToFastq"
include { LCA_TSV } from "../../../modules/local/lcaTsv"
include { SORT_FASTQ } from "../../../modules/local/sortFastq"
include { SORT_FILE } from "../../../modules/local/sortFile"
include { FILTER_VIRAL_SAM } from "../../../modules/local/filterViralSam"
include { PROCESS_LCA_ALIGNER_OUTPUT } from "../../../subworkflows/local/processLcaAlignerOutput/"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_SHORT {
    take:
        reads_ch
        params_map
    main:
        // Extract parameters from map
        ref_dir = params_map.ref_dir
        aln_score_threshold = params_map.aln_score_threshold
        adapter_path = params_map.adapter_path
        cutadapt_error_rate = params_map.cutadapt_error_rate
        min_kmer_hits = params_map.min_kmer_hits
        k = params_map.k
        bbduk_suffix = params_map.bbduk_suffix
        taxid_artificial = params_map.taxid_artificial
        
        // Get reference paths
        viral_genome_path = "${ref_dir}/results/virus-genomes-masked.fasta.gz"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        bt2_virus_index_path = "${ref_dir}/results/bt2-virus-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        nodes_db = "${ref_dir}/results/taxonomy-nodes.dmp"
        names_db = "${ref_dir}/results/taxonomy-names.dmp"
       // Define columns to keep, separating by ones to prefix and ones to not
        col_keep_no_prefix = ["seq_id", "sample", "aligner_taxid_lca", "aligner_taxid_top", 
                              "aligner_length_normalized_score_mean", "aligner_taxid_lca_combined",
                              "aligner_n_assignments_combined", "aligner_length_normalized_score_mean_combined",
                              "aligner_taxid_lca_artificial", "aligner_n_assignments_artificial", 
                              "aligner_length_normalized_score_mean_artificial", "query_len", "query_len_rev",
                              "query_seq", "query_seq_rev", "query_qual", "query_qual_rev"]
        col_keep_add_prefix = ["genome_id_all", "taxid_all", "fragment_length", 
                               "best_alignment_score", "best_alignment_score_rev",
                               "edit_distance", "edit_distance_rev", "ref_start", 
                               "ref_start_rev", "query_rc", "query_rc_rev", "pair_status"]
         // 1. Run initial screen against viral genomes with BBDuk
        bbduk_params = [
            min_kmer_hits: min_kmer_hits,
            k: k,
            suffix: bbduk_suffix
        ]
        bbduk_ch = BBDUK_HITS(reads_ch, viral_genome_path, bbduk_params)
        // 2. Carry out stringent adapter removal with FASTP and Cutadapt
        fastp_ch = FASTP(bbduk_ch.fail, adapter_path, true)
        adapt_ch = CUTADAPT(fastp_ch.reads, adapter_path, cutadapt_error_rate)
        // 3. Run Bowtie2 against a viral database and process output
        par_virus = "--local --very-sensitive-local --score-min G,0.1,19 -k 10"
        bowtie2_virus_params = [
            par_string: par_virus,
            suffix: "virus",
            remove_sq: true,
            debug: false,
            interleaved: true
        ]
        bowtie2_ch = BOWTIE2_VIRUS(adapt_ch.reads, bt2_virus_index_path, bowtie2_virus_params)
        // 4. Filter contaminants
        par_contaminants = "--local --very-sensitive-local"
        bowtie2_human_params = [
            par_string: par_contaminants,
            suffix: "human",
            remove_sq: false,
            debug: false,
            interleaved: true
        ]
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_ch.reads_mapped, bt2_human_index_path, bowtie2_human_params)
        bowtie2_other_params = [
            par_string: par_contaminants,
            suffix: "other",
            remove_sq: false,
            debug: false,
            interleaved: true
        ]
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unmapped, bt2_other_index_path, bowtie2_other_params)
        // 5. Sort SAM and FASTQ files before filtering
        bowtie2_sam_sorted_ch = SORT_FILE(bowtie2_ch.sam, "-t\$\'\\t\' -k1,1", "sam")
        other_fastq_sorted_ch = SORT_FASTQ(other_bt2_ch.reads_unmapped)
        // 6. Consolidated viral SAM filtering: keep contaminant-free reads, applies score threshold, adds missing mates
        bowtie2_ch_combined = bowtie2_sam_sorted_ch.output.combine(other_fastq_sorted_ch.output, by: 0)
        bowtie2_filtered_ch = FILTER_VIRAL_SAM(bowtie2_ch_combined, aln_score_threshold)
        // 7. Convert SAM to TSV
        bowtie2_tsv_ch = PROCESS_VIRAL_BOWTIE2_SAM(bowtie2_filtered_ch.sam, genome_meta_path, virus_db_path, true)
        // 8. Run LCA
        lca_params = [
            group_field: "seq_id",
            taxid_field: "taxid",
            score_field: "length_normalized_score",
            taxid_artificial: taxid_artificial,
            prefix: "aligner"
        ]
        lca_ch = LCA_TSV(bowtie2_tsv_ch.output, nodes_db, names_db, lca_params)
        // 9. Process LCA and Bowtie2 columns
        lca_aligner_params = [
            col_keep_no_prefix: col_keep_no_prefix,
            col_keep_add_prefix: col_keep_add_prefix,
            column_prefix: "prim_align_"
        ]
        processed_ch = PROCESS_LCA_ALIGNER_OUTPUT(
            lca_ch.output,
            bowtie2_tsv_ch.output,
            lca_aligner_params
        )
        // 10. Extract filtered virus hits in FASTQ format
        fastq_unfiltered_collect = other_bt2_ch.reads_unmapped.map{ _sample, file -> file }.collect().ifEmpty([])
        fastq_unfiltered_concat = CONCATENATE_FILES(fastq_unfiltered_collect, "reads_unfiltered", "fastq.gz")
        fastq_ch = EXTRACT_VIRAL_HITS_TO_FASTQ(processed_ch.viral_hits_tsv, fastq_unfiltered_concat.output)
    emit:
        bbduk_match = bbduk_ch.fail
        bbduk_trimmed = adapt_ch.reads
        hits_final = processed_ch.viral_hits_tsv
        inter_lca = processed_ch.lca_tsv
        inter_bowtie = processed_ch.aligner_tsv
        hits_prelca = bowtie2_tsv_ch.output
        hits_fastq = fastq_ch.fastq
        test_reads = other_bt2_ch.reads_unmapped
        test_filt_bowtie = bowtie2_filtered_ch.sam
        test_unfilt_bowtie = bowtie2_ch.sam
}
