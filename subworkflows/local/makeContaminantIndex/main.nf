/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME as DOWNLOAD_COW } from "../../../modules/local/downloadGenome"
include { DOWNLOAD_GENOME as DOWNLOAD_PIG } from "../../../modules/local/downloadGenome"
include { DOWNLOAD_GENOME as DOWNLOAD_MOUSE } from "../../../modules/local/downloadGenome"
include { DOWNLOAD_GENOME as DOWNLOAD_CARP } from "../../../modules/local/downloadGenome"
include { DOWNLOAD_GENOME as DOWNLOAD_ECOLI } from "../../../modules/local/downloadGenome"
include { CONCATENATE_FASTA_GZIPPED } from "../../../modules/local/concatenateFasta"
include { BBMAP_INDEX } from "../../../modules/local/bbmap"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2"

/***********
| WORKFLOW |
***********/

workflow MAKE_CONTAMINANT_INDEX {
    take:
        cow_url
        pig_url
        mouse_url
        carp_url
        ecoli_url
        contaminants_path
    main:
        // Download reference genomes
        cow_ch = DOWNLOAD_COW(cow_url, "cow_genome")
        pig_ch = DOWNLOAD_PIG(pig_url, "pig_genome")
        mouse_ch = DOWNLOAD_MOUSE(mouse_url, "mouse_genome")
        carp_ch = DOWNLOAD_CARP(carp_url, "carp_genome")
        ecoli_ch = DOWNLOAD_ECOLI(ecoli_url, "ecoli_genome")
        ref_ch = cow_ch.mix(pig_ch, mouse_ch, carp_ch, ecoli_ch, channel.fromPath(contaminants_path)).collect()
        // Concatenate
        genome_ch = CONCATENATE_FASTA_GZIPPED(ref_ch, "ref_concat")
        // Make indexes
        bbmap_ch = BBMAP_INDEX(genome_ch, "bbm-other-index")
        bowtie2_ch = BOWTIE2_INDEX(genome_ch, "bt2-other-index")
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
