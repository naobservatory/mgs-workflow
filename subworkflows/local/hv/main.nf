/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BOWTIE2 as BOWTIE2_HV } from "../../../modules/local/bowtie2" addParams(suffix: "hv")
include { BOWTIE2 as BOWTIE2_HUMAN } from "../../../modules/local/bowtie2" addParams(suffix: "human")
include { BOWTIE2 as BOWTIE2_OTHER } from "../../../modules/local/bowtie2" addParams(suffix: "other")
include { PROCESS_BOWTIE2_SAM_PAIRED } from "../../../modules/local/processBowtie2Sam"
include { EXTRACT_UNCONC_READ_IDS } from "../../../modules/local/extractUnconcReadIDs"
include { EXTRACT_UNCONC_READS } from "../../../modules/local/extractUnconcReads"
include { COMBINE_MAPPED_BOWTIE2_READS } from "../../../modules/local/combineMappedBowtie2Reads"
include { BBMAP as BBMAP_HUMAN } from "../../../modules/local/bbmap" addParams(suffix: "human")
include { BBMAP as BBMAP_OTHER } from "../../../modules/local/bbmap" addParams(suffix: "other")
include { TAXONOMY } from "../../../subworkflows/local/taxonomy" addParams(dedup_rc: true, classification_level: "F", read_fraction: 1)
include { PROCESS_KRAKEN_HV } from "../../../modules/local/processKrakenHV"
include { MERGE_SAM_KRAKEN } from "../../../modules/local/mergeSamKraken"
include { MERGE_TSVS } from "../../../modules/local/mergeTsvs" addParams(name: "bowtie2_kraken_merged")
include { FILTER_HV } from "../../../modules/local/filterHV"
include { COLLAPSE_HV } from "../../../modules/local/collapseHV"
include { MAKE_HV_FASTA } from "../../../modules/local/makeHvFasta"
include { COUNT_HV_CLADES } from "../../../modules/local/countHvClades"
include { BBDUK_HITS } from "../../../modules/local/bbduk" addParams(suffix: params.bbduk_suffix)
include { CUTADAPT } from "../../../modules/local/cutadapt"
include { TRIMMOMATIC } from "../../../modules/local/trimmomatic" addParams(encoding: params.encoding)

/***********
| WORKFLOW |
***********/

workflow HV {
    take:
        reads_ch
        ref_dir
        kraken_db_ch
        aln_score_threshold
        adapter_path
    main:
        // Get reference paths
        hv_ref_path = "${ref_dir}/results/human-viral-genomes-filtered.fasta.gz"
        genomeid_map_path = "${ref_dir}/results/genomeid-to-taxid.json"
        bt2_hv_index_path = "${ref_dir}/results/bt2-hv-index"
        bt2_human_index_path = "${ref_dir}/results/bt2-human-index"
        bt2_other_index_path = "${ref_dir}/results/bt2-other-index"
        bbm_human_index_path = "${ref_dir}/results/bbm-human-index"
        bbm_other_index_path = "${ref_dir}/results/bbm-other-index"
        nodes_path = "${ref_dir}/results/taxonomy-nodes.dmp"
        hv_db_path = "${ref_dir}/results/human-virus-db.tsv.gz"
        viral_taxa_path = "${ref_dir}/results/total-virus-db.tsv.gz"
        // Run initial screen against HV genomes with BBDuk
        bbduk_ch = BBDUK_HITS(reads_ch, hv_ref_path, params.min_kmer_hits, params.k)
        // Carry out stringent adapter removal with Cutadapt and Trimmomatic
        adapt_ch = CUTADAPT(bbduk_ch.fail, adapter_path)
        trim_ch = TRIMMOMATIC(adapt_ch.reads, adapter_path)
        // Run Bowtie2 against an HV database and process output
        bowtie2_ch = BOWTIE2_HV(trim_ch.reads, bt2_hv_index_path, "--no-unal --no-sq --score-min G,1,1")
        bowtie2_sam_ch = PROCESS_BOWTIE2_SAM_PAIRED(bowtie2_ch.sam, genomeid_map_path)
        bowtie2_unconc_ids_ch = EXTRACT_UNCONC_READ_IDS(bowtie2_ch.sam)
        bowtie2_unconc_reads_ch = EXTRACT_UNCONC_READS(bowtie2_ch.reads_unconc.combine(bowtie2_unconc_ids_ch, by: 0))
        bowtie2_reads_combined_ch = COMBINE_MAPPED_BOWTIE2_READS(bowtie2_ch.reads_conc.combine(bowtie2_unconc_reads_ch, by: 0))
        // Filter contaminants
        human_bt2_ch = BOWTIE2_HUMAN(bowtie2_reads_combined_ch, bt2_human_index_path, "")
        other_bt2_ch = BOWTIE2_OTHER(human_bt2_ch.reads_unconc, bt2_other_index_path, "")
        human_bbm_ch = BBMAP_HUMAN(other_bt2_ch.reads_unconc, bbm_human_index_path)
        other_bbm_ch = BBMAP_OTHER(human_bbm_ch.reads_unmapped, bbm_other_index_path)
        // Run Kraken on filtered HV candidates
        tax_ch = TAXONOMY(other_bbm_ch.reads_unmapped, kraken_db_ch)
        // Process Kraken output and merge with Bowtie2 output across samples
        kraken_output_ch = PROCESS_KRAKEN_HV(tax_ch.kraken_output, nodes_path, hv_db_path)
        bowtie2_kraken_merged_ch = MERGE_SAM_KRAKEN(kraken_output_ch.combine(bowtie2_sam_ch, by: 0))
        merged_ch = MERGE_TSVS(bowtie2_kraken_merged_ch.collect().ifEmpty([]))
        // Filter and process putative HV hit TSV
        filtered_ch = FILTER_HV(merged_ch, aln_score_threshold)
        collapsed_ch = COLLAPSE_HV(filtered_ch)
        fasta_ch = MAKE_HV_FASTA(collapsed_ch)
        // Count clades
        count_ch = COUNT_HV_CLADES(collapsed_ch, viral_taxa_path)
    emit:
        tsv = collapsed_ch
        fasta = fasta_ch
        counts = count_ch
        merged_summary = tax_ch.merged_summary
}
