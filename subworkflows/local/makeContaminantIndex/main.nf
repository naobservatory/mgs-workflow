/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME as DOWNLOAD_COW } from "../../../modules/local/downloadGenome" addParams(name: "cow_genome")
include { DOWNLOAD_GENOME as DOWNLOAD_PIG } from "../../../modules/local/downloadGenome" addParams(name: "pig_genome")
include { DOWNLOAD_GENOME as DOWNLOAD_MOUSE } from "../../../modules/local/downloadGenome" addParams(name: "mouse_genome")
include { DOWNLOAD_GENOME as DOWNLOAD_CARP } from "../../../modules/local/downloadGenome" addParams(name: "carp_genome")
include { DOWNLOAD_GENOME as DOWNLOAD_ECOLI } from "../../../modules/local/downloadGenome" addParams(name: "ecoli_genome")
include { CONCATENATE_FASTA_GZIPPED } from "../../../modules/local/concatenateFasta" addParams(name: "ref_concat")
include { BBMAP_INDEX } from "../../../modules/local/bbmap" addParams(outdir: "bbm-other-index")
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2" addParams(outdir: "bt2-other-index")

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
        cow_ch = DOWNLOAD_COW(cow_url)
        pig_ch = DOWNLOAD_PIG(pig_url)
        mouse_ch = DOWNLOAD_MOUSE(mouse_url)
        carp_ch = DOWNLOAD_CARP(carp_url)
        ecoli_ch = DOWNLOAD_ECOLI(ecoli_url)
        ref_ch = cow_ch.mix(pig_ch, mouse_ch, carp_ch, ecoli_ch, channel.fromPath(contaminants_path)).collect()
        // Concatenate
        genome_ch = CONCATENATE_FASTA_GZIPPED(ref_ch)
        // Make indexes
        bbmap_ch = BBMAP_INDEX(genome_ch)
        bowtie2_ch = BOWTIE2_INDEX(genome_ch)
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
