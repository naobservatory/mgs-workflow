// Version of EXTRACT_VIRAL_READS that uses streaming and interleaved files to minimize memory requirements and loading times

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_HITS_STREAMED } from "../../../modules/local/bbduk"
include { CUTADAPT_STREAMED } from "../../../modules/local/cutadapt"
include { BOWTIE2_STREAMED as BOWTIE2_VIRUS } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2"
include { BOWTIE2_STREAMED as BOWTIE2_OTHER } from "../../../modules/local/bowtie2"

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
    emit:
        test_out = other_bt2_ch.reads_unmapped
}
