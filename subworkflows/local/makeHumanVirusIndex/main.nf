/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_HUMAN_VIRUS_GENOMES } from "../../../modules/local/downloadHumanVirusGenomes"
include { MAKE_HUMAN_VIRUS_ID_MAPPING } from "../../../modules/local/makeHumanVirusIdMapping"
include { CONCATENATE_FASTA_GZIPPED_DIR_DEEP } from "../../../modules/local/concatenateFasta" addParams(name: "human-viral-genomes", suffix: "fna.gz")
include { FILTER_HUMAN_VIRUS_GENOMES } from "../../../modules/local/filterHumanVirusGenomes" 
include { DUSTMASKER_FASTA_GZIPPED } from "../../../modules/local/dustmasker"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2" addParams(outdir: "bt2-hv-index")

/***********
| WORKFLOW |
***********/

workflow MAKE_HUMAN_VIRUS_INDEX {
    take:
        hv_taxid_ch
        hv_patterns_exclude
    main:
        genomes_ch = DOWNLOAD_HUMAN_VIRUS_GENOMES(hv_taxid_ch)
        mapping_ch = MAKE_HUMAN_VIRUS_ID_MAPPING(genomes_ch.metadata, genomes_ch.genomes)
        concat_ch = CONCATENATE_FASTA_GZIPPED_DIR_DEEP(genomes_ch.genomes)
        filtered_ch = FILTER_HUMAN_VIRUS_GENOMES(concat_ch, hv_patterns_exclude)
        masked_ch = DUSTMASKER_FASTA_GZIPPED(filtered_ch)
        bowtie2_ch = BOWTIE2_INDEX(masked_ch)
    emit:
        filtered = filtered_ch
        bt2 = bowtie2_ch
        mapping = mapping_ch
}
