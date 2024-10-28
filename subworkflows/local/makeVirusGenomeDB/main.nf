/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_VIRAL_GENBANK } from "../../../modules/local/downloadViralGenbank"
include { FILTER_VIRAL_GENBANK_METADATA } from "../../../modules/local/filterViralGenbankMetadata" addParams(name: "virus-genome")
include { CONCATENATE_GENOME_FASTA } from "../../../modules/local/concatenateGenomeFasta"
include { FILTER_GENOME_FASTA } from "../../../modules/local/filterGenomeFasta" addParams(name: "virus-genomes-filtered")

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_GENOME_DB {
    take:
        virus_db // TSV giving taxonomic structure and host infection status of virus taxids
        patterns_exclude // File of sequence header patterns to exclude from genome DB
        host_taxa // Tuple of host taxa to include
    main:
        // 1. Download viral Genbank
        dl_ch = DOWNLOAD_VIRAL_GENBANK()
        // 2. Filter genome metadata by taxid to identify genomes to retain
        meta_ch = FILTER_VIRAL_GENBANK_METADATA(dl_ch.metadata, virus_db, host_taxa)
        // 3. Concatenate matching genomes
        concat_ch = CONCATENATE_GENOME_FASTA(dl_ch.genomes, meta_ch.gid)
        // 4. Filter to remove undesired/contaminated genomes
        filter_ch = FILTER_GENOME_FASTA(concat_ch, patterns_exclude)
    emit:
        fasta = filter_ch
        metadata = meta_ch.db
}
