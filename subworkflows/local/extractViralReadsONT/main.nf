/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { MERGE_SAM as MERGE_SAM_VIRUS } from "../../../modules/local/samtools"
include { MERGE_SAM as MERGE_SAM_HUMAN } from "../../../modules/local/samtools"
include { CONCATENATE_TSVS as CONCATENATE_HV_TSVS } from "../../../modules/local/concatenateTsvs"
include { FILTLONG } from "../../../modules/local/filtlong"
include { MASK_FASTQ_READS } from "../../../modules/local/maskRead"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
include { BLAST_VIRAL } from "../../../subworkflows/local/blastViral"
include { PULLOUT_FASTQ } from "../../../modules/local/pulloutFastq"
include { CONCATENATE_FASTQ_GZIPPED } from "../../../modules/local/concatenateFastq"
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
        minimap2_virus_index = "${ref_dir}/results/mm2-virus-index"
        minimap2_human_index = "${ref_dir}/results/mm2-human-index"
        minimap2_contam_index = "${ref_dir}/results/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"
        // blast_db_path = "${ref_dir}/results/${blast_db_prefix}"

        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 50, 90)

        // Mask non-complex read sections
        masked_ch = MASK_FASTQ_READS(filtered_ch, 25, 0.55)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
        // contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
        // no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_human_ch, minimap2_virus_index, "hv", false)

        virus_sam_ch = virus_minimap2_ch.sam
        virus_fastq_ch = virus_minimap2_ch.reads_mapped

        // Pull out clean reads from mapped reads
        clean_matched_subset_ch = PULLOUT_FASTQ(virus_fastq_ch.join(filtered_ch.reads))

        // Create common channel for HV SAM and clean HV reads
        sam_and_reads_ch = virus_sam_ch.join(clean_matched_subset_ch.output)

        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(sam_and_reads_ch, genome_meta_path, virus_db_path, host_taxon)

        // Concatenate HV TSVs
        hv_tsvs = hv_tsv_ch.output.map { it[1] }.collect()
        merged_tsv_ch = CONCATENATE_HV_TSVS(hv_tsvs, "hv")

    emit:
        hv_tsv = merged_tsv_ch.output
        // blast_subset = blast_subset_ch
        // blast_reads = blast_reads_ch
}
