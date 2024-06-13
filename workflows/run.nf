/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { EXTRACT_TARBALL } from "../modules/local/extractTarball"
include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { CLEAN } from "../subworkflows/local/clean" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "cleaned")
include { DEDUP } from "../subworkflows/local/dedup" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "dedup")
include { RIBODEPLETION as RIBO_INITIAL } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_initial", min_kmer_fraction: "0.6", k: "43")
include { RIBODEPLETION as RIBO_SECONDARY } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_secondary", min_kmer_fraction: "0.4", k: "27")
include { TAXONOMY } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D")
include { HV } from "../subworkflows/local/hv"
include { BLAST_HV } from "../subworkflows/local/blast_hv" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")

/************************************
| 0. PREPARE REFERENCES AND INDEXES |
************************************/

// 0.8. Copy sample metadata CSV
process COPY_SAMPLE_METADATA {
    label "BBTools"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    errorStrategy "retry"
    input:
        path(sample_metadata)
    output:
        path("sample-metadata.csv")
    shell:
        '''
        cp !{sample_metadata} sample-metadata.csv
        '''
}

/*****************
| MAIN WORKFLOWS |
*****************/

// Prepare libraries
libraries_ch = Channel
    .fromPath(params.library_tab)
    .splitCsv(header: true)
    .map{row -> [row.sample, row.library]}
    .groupTuple()

// Complete primary workflow
workflow {
    // Prepare references & indexes
    PREPARE_REFERENCES(params.human_index_bt2, params.other_index_bt2, params.human_index_bb, params.other_index_bb, params.hv_index, params.kraken_db, params.sample_tab)
    kraken_db_ch = EXTRACT_TARBALL(params.kraken_db) // TODO: Eliminate this step and just pass reference db path directly
    // Preprocessing
    RAW(libraries_ch, params.raw_dir, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)
    DEDUP(CLEAN_READS.out.reads)
    RIBO_INITIAL(DEDUP.out.reads, params.ribo_db)
    RIBO_SECONDARY(RIBO_INITIAL.out.reads, params.ribo_db)
    // Taxonomic profiling (all ribodepleted)
    TAX = TAXONOMY(RIBO_SECONDARY.out.reads, kraken_db_ch, 1)
    // Taxonomic profiling (pre/post dedup)
    TAX_PRE = TAXONOMY(CLEAN.out.reads, kraken_db_ch, params.classify_dedup_subset)
    TAX_POST = TAXONOMY(DEDUP.out.reads, kraken_db_ch, params.classify_dedup_subset)
    // Extract and count human-viral reads
    HV(RIBO_INITIAL.out.reads, params.genomeid_map, params.hv_index, params.human_index_bt2, params.other_index_bt2,
       params.human_index_bb, params.other_index_bb, kraken_db_ch, params.nodes, params.hv_db,
       params.bt2_score_threshold, params.viral_taxa_db)
    // BLAST validation on human-viral reads (if activated)
    if ( params.blast_hv ) {
        BLAST_HV(HV.out.fasta, params.blast_nt_dir, params.blast_fraction)
    }
    // Process output
    PROCESS_OUTPUT(RAW.out.qc, CLEAN.out.qc, DEDUP.out.qc, RIBO_INITIAL.out.qc, RIBO_SECONDARY.out.qc, TAX.bracken, TAX_PRE.bracken, TAX_POST.bracken, params.classify_dedup_subset)
    // Publish results
}

output {
    directory "${params.pub_dir}"
    mode "copy"
}
