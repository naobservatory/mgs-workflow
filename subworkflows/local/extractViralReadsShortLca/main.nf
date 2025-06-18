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
include { PROCESS_VIRAL_BOWTIE2_SAM_LCA as PROCESS_VIRAL_BOWTIE2_SAM } from "../../../modules/local/processViralBowtie2SamLca"
include { SORT_TSV as SORT_BOWTIE_VIRAL } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_LCA } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_BOWTIE2 } from "../../../modules/local/sortTsv"
include { ADD_SAMPLE_COLUMN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { EXTRACT_VIRAL_HITS_TO_FASTQ } from "../../../modules/local/extractViralHitsToFastq"
include { LCA_TSV } from "../../../modules/local/lcaTsv"
include { SORT_FASTQ } from "../../../modules/local/sortFastq"
include { SORT_FILE } from "../../../modules/local/sortFile"
include { FILTER_VIRAL_SAM } from "../../../modules/local/filterViralSam"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { FILTER_TO_PRIMARY_ALIGNMENTS } from "../../../modules/local/filterToPrimaryAlignments"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_SHORT_LCA {
    take:
        reads_ch
        ref_dir
        aln_score_threshold
        adapter_path
        host_taxon
        cutadapt_error_rate
        min_kmer_hits
        k
        bbduk_suffix
        taxid_artificial
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
            par_virus, "virus", true, false, true)
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
            "seq_id", "aligner_taxid", "aligner_length_normalized_score", taxid_artificial,
            "bowtie2")
        // 9. Sort both TSV files by seq_id for joining
        lca_sorted_ch = SORT_LCA(lca_ch.output, "seq_id")
        bowtie2_sorted_ch = SORT_BOWTIE2(bowtie2_tsv_ch.output, "seq_id")
        // 10. Join LCA and Bowtie2 TSV on seq_id to combine taxonomic and alignment data
        joined_input_ch = lca_sorted_ch.sorted.join(bowtie2_sorted_ch.sorted, by: 0)
        joined_ch = JOIN_TSVS(joined_input_ch, "seq_id", "inner", "lca_bowtie2")
        // 11. Filter to keep only primary alignments (is_secondary=False)
        filtered_ch = FILTER_TO_PRIMARY_ALIGNMENTS(joined_ch.output)
        out_labeled_ch = ADD_SAMPLE_COLUMN(filtered_ch.output, "sample", "viral_bowtie2")
        // 12. Concatenate across reads
        label_combined_ch = out_labeled_ch.output.map{ sample, file -> file }.collect().ifEmpty([])
        concat_ch = CONCATENATE_TSVS(label_combined_ch, "virus_hits_final")
        // 13. Extract filtered virus hits in FASTQ format
        fastq_unfiltered_collect = other_bt2_ch.reads_unmapped.map{ sample, file -> file }.collect().ifEmpty([])
        fastq_unfiltered_concat = CONCATENATE_FILES(fastq_unfiltered_collect, "reads_unfiltered", "fastq.gz")
        fastq_ch = EXTRACT_VIRAL_HITS_TO_FASTQ(concat_ch.output, fastq_unfiltered_concat.output)
    emit:
        bbduk_match = bbduk_ch.fail
        bbduk_trimmed = adapt_ch.reads
        hits_final = concat_ch.output
        hits_prelca = bowtie2_tsv_ch.output
        hits_fastq = fastq_ch.fastq
        test_reads  = other_bt2_ch.reads_unmapped
        test_filt_bowtie = bowtie2_filtered_ch.sam
        test_unfilt_bowtie = bowtie2_ch.sam
}
