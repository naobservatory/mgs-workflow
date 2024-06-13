/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { EXTRACT_TARBALL } from "../modules/local/extractTarball"
include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { CLEAN } from "../subworkflows/local/clean" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "cleaned")
include { DEDUP } from "../subworkflows/local/dedup" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "dedup")
include { RIBODEPLETION as RIBO_INITIAL } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_initial", min_kmer_fraction: "0.6", k: "43")
include { RIBODEPLETION as RIBO_SECONDARY } from "../subworkflows/local/ribodepletion" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_secondary", min_kmer_fraction: "0.4", k: "27")
include { TAXONOMY as TAXONOMY_FULL } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: 1)
include { TAXONOMY as TAXONOMY_PRE } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: params.classify_dedup_subset)
include { TAXONOMY as TAXONOMY_POST } from "../subworkflows/local/taxonomy" addParams(dedup_rc: false, classification_level: "D", read_fraction: params.classify_dedup_subset)
include { HV } from "../subworkflows/local/hv"
include { BLAST_HV } from "../subworkflows/local/blast_hv" addParams(blast_cpus: "32", blast_mem: "256 GB", blast_filter_mem: "32 GB")
nextflow.preview.output = true

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
    kraken_db_ch = EXTRACT_TARBALL(params.kraken_db) // TODO: Eliminate this step and just pass reference db path directly
    // Preprocessing
    RAW(libraries_ch, params.raw_dir, params.n_reads_trunc)
    CLEAN(RAW.out.reads, params.adapters)
    DEDUP(CLEAN.out.reads)
    RIBO_INITIAL(DEDUP.out.reads, params.ribo_db)
    RIBO_SECONDARY(RIBO_INITIAL.out.reads, params.ribo_db)
    // Taxonomic profiling (all ribodepleted)
    TAXONOMY_FULL(RIBO_SECONDARY.out.reads, kraken_db_ch)
    // Taxonomic profiling (pre/post dedup)
    TAXONOMY_PRE(CLEAN.out.reads, kraken_db_ch)
    TAXONOMY_POST(DEDUP.out.reads, kraken_db_ch)
    // Extract and count human-viral reads
    HV(RIBO_INITIAL.out.reads, params.genomeid_map, params.hv_index, params.human_index_bt2, params.other_index_bt2,
       params.human_index_bb, params.other_index_bb, kraken_db_ch, params.nodes, params.hv_db,
       params.bt2_score_threshold, params.viral_taxa_db)
    // BLAST validation on human-viral reads (if activated)
    if ( params.blast_hv ) {
        BLAST_HV(HV.out.fasta, params.blast_nt_dir, params.blast_fraction)
    }
    // Process output
    PROCESS_OUTPUT(RAW.out.qc, CLEAN.out.qc, DEDUP.out.qc, RIBO_INITIAL.out.qc, RIBO_SECONDARY.out.qc,
                   TAXONOMY_FULL.out.bracken, TAXONOMY_PRE.out.bracken, TAXONOMY_POST.out.bracken, params.classify_dedup_subset)
    // Save parameters
    params_ch = Channel.of(params).collectFile(name: "params.txt", newLine: true) // TODO: Iterate on this
    // Publish results
    publish:
        // Saved inputs
        params_ch >> "input"
        params.sample_tab >> "input"
        params.adapters >> "input"
        // Intermediate files
        CLEAN.out.reads >> "intermediates/cleaned_reads"
        // Final results
        HV.out.tsv >> "results"
        HV.out.counts >> "results"
        PROCESS_OUTPUT.out.composition_full >> "results/taxonomy_final"
        PROCESS_OUTPUT.out.composition_pre >> "results/taxonomy_pre_dedup"
        PROCESS_OUTPUT.out.composition_post >> "results/taxonomy_post_dedup"
}

output {
    directory "${params.pub_dir}"
    mode "copy"
}
