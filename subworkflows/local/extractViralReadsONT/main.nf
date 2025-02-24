/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { MERGE_SAM as MERGE_SAM_VIRUS } from "../../../modules/local/samtools"
include { MERGE_SAM as MERGE_SAM_HUMAN } from "../../../modules/local/samtools"
include { FILTLONG } from "../../../modules/local/filtlong"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
include { BLAST_VIRAL } from "../../../subworkflows/local/blastViral"
include { PULLOUT_FASTQ } from "../../../modules/local/pulloutFastq"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_ONT {
    take:
        reads_ch
        ref_dir
        host_taxon
    main:
        // Get reference_paths
        minimap2_virus_index = "s3://nao-mgs-simon/ont-indices/2024-12-14/minimap2-hv-index"
        minimap2_human_index = "s3://nao-mgs-simon/ont-indices/2024-12-14/minimap2-human-index"
        // minimap2_contam_index = "s3://nao-mgs-simon/ont-indices/2024-12-14/minimap2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        blast_db_path = "${params.ref_dir}/results/${params.blast_db_prefix}"

        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 50, 10000, 90)

        // Mark non-complex read sections
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(filtered_ch)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
            // contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
            // no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_human_ch, minimap2_virus_index, "virus_1", false)

        // Pull out SAM files only
        virus_sam_ch = virus_minimap2_ch.sam.map { it[1] }

        // Pull out FASTQ files only
        // virus_fastq_ch = virus_minimap2_ch.reads_mapped.map { it[1] }

        // BLAST virus reads
        // blast_ch = BLAST_VIRAL(virus_fastq_ch, blast_db_path, params.blast_db_prefix, params.blast_viral_fraction, params.blast_max_rank, params.blast_min_frac, params.random_seed)

        // blast_subset_ch = BLAST_VIRAL.out.blast_subset
        // blast_reads_ch = BLAST_VIRAL.out.subset_reads

        // Merge SAM files
        merged_virus_sam_ch = MERGE_SAM_VIRUS(virus_sam_ch.collect(), "hv")

        // Get original non-masked reads
        // original_reads_ch = masked_ch.input.map { it[1] }
        // original_reads_ch = EXTRACT_FASTQ(original_reads_ch, blast_subset_ch)

        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(merged_virus_sam_ch, genome_meta_path, virus_db_path, host_taxon)

    emit:
        hv_tsv = hv_tsv_ch.output
        // blast_subset = blast_subset_ch
        // blast_reads = blast_reads_ch
}
