/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_VIRAL_NCBI } from "../../../modules/local/downloadViralNCBI"
include { FILTER_VIRAL_GENBANK_METADATA } from "../../../modules/local/filterViralGenbankMetadata"
include { ADD_GENBANK_GENOME_IDS } from "../../../modules/local/addGenbankGenomeIDs"
include { CONCATENATE_GENOME_FASTA } from "../../../modules/local/concatenateGenomeFasta"
include { FILTER_GENOME_FASTA } from "../../../modules/local/filterGenomeFasta"

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_GENOME_DB {
    take:
        ncbi_viral_params // Parameters for downloading genomic sequences from NCBI's GenBank or RefSeq databases. Any additional sequence-specific options can be provided here, supplementing the default download settings.
        virus_db // TSV giving taxonomic structure and host infection status of virus taxids
        patterns_exclude // File of sequence header patterns to exclude from genome DB
        host_taxa // Tuple of host taxa to include
    main:
        // 1. Download viral Genbank
        dl_ch = DOWNLOAD_VIRAL_NCBI(ncbi_viral_params)
        // 2. Filter genome metadata by taxid to identify genomes to retain
        meta_ch = FILTER_VIRAL_GENBANK_METADATA(dl_ch.metadata, virus_db, host_taxa, "virus-genome")
        // 3. Add genome IDs to Genbank metadata file
        gid_ch = ADD_GENBANK_GENOME_IDS(meta_ch.db, dl_ch.genomes, "virus-genome")
        // 4. Concatenate matching genomes
        concat_ch = CONCATENATE_GENOME_FASTA(dl_ch.genomes, meta_ch.path)
        // 5. Filter to remove undesired/contaminated genomes
        filter_ch = FILTER_GENOME_FASTA(concat_ch, patterns_exclude, "virus-genomes-filtered")
    emit:
        fasta = filter_ch
        metadata = gid_ch
}
