# Configurations 

Configurations are specified by `.config` files. These files are used to specify parameters and other configuration options used by Nextflow in executing the pipeline. 

All configurations are stored in the `configs` directory.

# Workflow-specific configurations

## Run workflow (`configs/run.config`)

- `params.mode = "run"`: This instructs the pipeline to execute the [core run workflow](../workflows/run.nf).
- `params.base_dir`: The parent directory for the pipeline working and output directories.
- `params.ref_dir`: The directory containing the outputs of the `index` workflow.
- `params.sample_sheet`: The path to the sample sheet.
- `params.adapters`: The path to the adapter file for adapter trimming.
- `params.grouping`: Whether to group samples by the `group` column in the sample sheet.
- `params.n_reads_profile`: The number of reads per sample to run through taxonomic profiling.
- `params.bt2_score_threshold`: The normalized score threshold for calling a host-infecting virus read (typically 15 or 20).
- `params.blast_viral_fraction`: The fraction of putative host-infecting virus reads to BLAST vs nt (0 = don't run BLAST).
- `params.quality_encoding`: The FASTQ quality encoding (probably phred33, maybe phred64).
- `params.fuzzy_match_alignment_duplicates`: Fuzzy matching the start coordinate of reads for identification of duplicates through alignment (0 = exact matching; options are 0, 1, or 2).
- `params.host_taxon`: The taxon to use for host-infecting virus identification.
- `params.blast_db_prefix`: The prefix for the BLAST database to use for host-infecting virus identification.

## Index workflow (`configs/index.config`)

- `params.mode = "index"`: This instructs the pipeline to execute the [index workflow](../workflows/index.nf).
- `params.base_dir`: The parent directory for the pipeline working and output directories.
- `params.taxonomy_url`: The URL for the NCBI taxonomy dump.
- `params.virus_host_db_url`: The URL for the viral host database.
- `params.human_url`: The URL used for the human genome, which is used to build the human index which we use for screening out human reads.
- `params.genome_urls`: The URLs for the genomes that are common contaminants, which is used to build out the other contaminants index which we use for screening out common contaminants.
- `params.ssu_url`: The URL for the SILVA SSU reference database, which we use for building the ribosomal index to classify reads as ribosomal or not.
- `params.lsu_url`: The URL for the SILVA LSU reference database, which we use for building the ribosomal index to classify reads as ribosomal or not.
- `params.host_taxon_db`: The path to the host taxon database, which we use to establish a host-taxon relationship for the viral genome database.
- `params.contaminants`: The path to the reads that are common contaminants, which we screen out for.
- `params.adapters`: The path to the common adapters which we use to clean the viral reference genomes that we screen for.
- `params.genome_patterns_exclude`: The path to the genome patterns to exclude.
- `params.kraken_db`: The path to the Kraken reference database which we use to taxonomically classify all reads.
- `params.blast_db_name`: The name of the BLAST database to use which is used during BLAST validation.
- `params.ncbi_viral_params`: The parameters to use for the NCBI viral genome database which we use to build out our viral indices for BBDuk and Bowtie2.
- `params.virus_taxid`: The taxid for the virus to use for building the viral genome database.
- `params.viral_taxids_exclude`: The taxids to exclude from the viral genome database.
- `params.host_taxa_screen`: The host taxa to screen for when building the viral genome database.
- `params.kraken_memory`: Initalizes the `run` workflow params to avoid warnings (DO NOT CHANGE)
- `params.classify_dedup_subset`: Initalizes the `run` workflow params to avoid warnings (DO NOT CHANGE)
