/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MINIMAP2 as MINIMAP2_VIRUS } from "../../../modules/local/minimap2"
include { MINIMAP2 as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
include { MINIMAP2_NON_STREAMED as MINIMAP2_CONTAM } from "../../../modules/local/minimap2"
include { CONCATENATE_TSVS as CONCATENATE_HV_TSVS } from "../../../modules/local/concatenateTsvs"
include { ADD_SAMPLE_COLUMN as LABEL_HV_TSVS } from "../../../modules/local/addSampleColumn"
include { FILTLONG } from "../../../modules/local/filtlong"
include { MASK_FASTQ_READS } from "../../../modules/local/maskRead"
include { PROCESS_VIRAL_MINIMAP2_SAM } from "../../../modules/local/processViralMinimap2Sam"
include { BLAST_VIRAL } from "../../../subworkflows/local/blastViral"
include { EXTRACT_SHARED_FASTQ_READS } from "../../../modules/local/extractSharedFastq"
include { CONCATENATE_FILES } from "../../../modules/local/concatenateFiles"
include { CUTADAPT_MANUAL_ADAPTER } from "../../../modules/local/cutadapt"
/***********
| WORKFLOW |
***********/

workflow EXTRACT_VIRAL_READS_ONT {
    take:
        reads_ch
        ref_dir
    main:
        // Get reference_paths
        minimap2_virus_index = "${ref_dir}/results/mm2-virus-index"
        minimap2_human_index = "${ref_dir}/results/mm2-human-index"
        minimap2_contam_index = "${ref_dir}/results/mm2-other-index"
        genome_meta_path = "${ref_dir}/results/virus-genome-metadata-gid.tsv.gz"
        virus_db_path = "${ref_dir}/results/total-virus-db-annotated.tsv.gz"

        // Filter reads by length and quality scores
        filtered_ch = FILTLONG(reads_ch, 50, 15000, 90)

        // Mask non-complex read sections
        masked_ch = MASK_FASTQ_READS(filtered_ch, 25, 0.55)

        // Drop human reads before pathogen identification
        human_minimap2_ch = MINIMAP2_HUMAN(masked_ch.masked, minimap2_human_index, "human", false)
        no_human_ch = human_minimap2_ch.reads_unmapped

        // Identify other contaminants
        contam_minimap2_ch = MINIMAP2_CONTAM(no_human_ch, minimap2_contam_index, "other", false)
        no_contam_ch = contam_minimap2_ch.reads_unmapped

        // Identify virus reads
        virus_minimap2_ch = MINIMAP2_VIRUS(no_human_ch, minimap2_virus_index, "hv", false)

        // Group cleaned reads and sam files by sample
        virus_sam_ch = virus_minimap2_ch.sam
        sam_fastq_ch = virus_sam_ch.join(filtered_ch)

        // Generate HV TSV
        hv_tsv_ch = PROCESS_VIRAL_MINIMAP2_SAM(sam_fastq_ch, genome_meta_path, virus_db_path)
        hv_tsv_labeled_ch = LABEL_HV_TSVS(hv_tsv_ch.output, "sample", "hv_tsv")

        // Concatenate HV TSVs
        hv_tsvs = hv_tsv_labeled_ch.output.map { it[1] }.collect()
        merged_tsv_ch = CONCATENATE_HV_TSVS(hv_tsvs, "hv")

        // Pull out clean reads from mapped reads to feed into BLAST
        virus_fastq_ch = virus_minimap2_ch.reads_mapped
        clean_virus_fastq_ch = EXTRACT_SHARED_FASTQ_READS(virus_fastq_ch.join(filtered_ch.reads))
        trimmed_virus_fastq_ch = CUTADAPT_MANUAL_ADAPTER(clean_virus_fastq_ch.output, "A{50}")
        merged_virus_fastq_ch = CONCATENATE_FILES(trimmed_virus_fastq_ch.reads.map{ it[1] }.collect(), "clean_virus_reads", "fastq.gz")

    emit:
        hits_hv = merged_tsv_ch.output
        hits_fastq = merged_virus_fastq_ch.output
        test_minimap2_virus = virus_sam_ch
        test_fastq_filtered_human = human_minimap2_ch.reads_unmapped
        test_fastq_filtered_contam = contam_minimap2_ch.reads_unmapped
}
