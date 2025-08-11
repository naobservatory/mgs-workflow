/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_VIRAL_NCBI } from "../../../modules/local/downloadViralNCBI"
include { FILTER_VIRAL_GENBANK_METADATA } from "../../../modules/local/filterViralGenbankMetadata"
include { ADD_GENBANK_GENOME_IDS } from "../../../modules/local/addGenbankGenomeIDs"
include { CONCATENATE_GENOME_FASTA } from "../../../modules/local/concatenateGenomeFasta"
include { FILTER_GENOME_FASTA } from "../../../modules/local/filterGenomeFasta"
include { MASK_GENOME_FASTA } from "../../../modules/local/maskGenomeFasta"

/***********
| WORKFLOW |
***********/

workflow MAKE_VIRUS_GENOME_DB {
    take:
        ncbi_viral_params // Parameters for downloading genomic sequences from NCBI's GenBank or RefSeq databases
        virus_db // TSV giving taxonomic structure and host infection status of virus taxids
        params_map // Map containing all parameters
    main:
        // Extract parameters from map
        patterns_exclude = params_map.patterns_exclude
        host_taxa = params_map.host_taxa
        adapters = params_map.adapters
        k = params_map.k
        hdist = params_map.hdist
        entropy = params_map.entropy
        polyx_len = params_map.polyx_len
        
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
	// 6. Mask to remove adapters, low-entropy regions, and polyX
	mask_params = [
		k: k,
		hdist: hdist,
		entropy: entropy,
		polyx_len: polyx_len,
		name_pattern: "virus-genomes"
	]
	mask_ch = MASK_GENOME_FASTA(filter_ch, adapters, mask_params)
    emit:
        fasta = mask_ch.masked
        metadata = gid_ch
}
