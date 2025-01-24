# Configurations 

Configuartions are specified by `.config` files. These files are used to specify parameters and other configuration options used by Nextflow in executing the pipeline. 

All configurations are stored in the `configs` directory.

# Workflow specific configurations

## Run workflow (`configs/run.config`)

- `params.mode = "run"`: This instructs the pipeline to execute the [core run workflow](./run.md).
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
