# v2.9.0.0-dev
- Combined run_dev_se.nf workflow with run.nf: ONT is now an implemented platform!
- Added end-to-end tests for ONT to github actions
- Added a development_mode parameter to LOAD_SAMPLESHEET to allow testing on non-implemented platform/endedness
- Rename `hits_filtered` output of short read pipeline to `hits_final` (so pipeline produces `virus_hits_final.fastq.gz`/`virus_hits_final.tsv.gz`)
    - Also rename `hits_all` output of short read pipeline to `hits_unfiltered`
- Rename `hits_hv` output of ONT pipeline to `hits_final` for consistency with short read pipeline
- Get rid of lingering references to human viruses/HV in comments, variable names, etc.
- Began development of post-hoc validation of putative viral hits in the DOWNSTREAM workflow
    - Split viral hits TSV by assigned species and extract read sequences (SPLIT_VIRAL_TSV_BY_SPECIES)
    - Cluster within species with VSEARCH and obtain representative sequences (CLUSTER_VIRAL_ASSIGNMENTS)
- Updated SUBSET_FASTQ to handle plaintext and FASTA input (and renamed to SUBSET_FASTN)
- Split out merge/join part of TAXONOMY workflow into its own subworkflow (MERGE_JOIN_READS) with associated tests
- Added `test_component_dependencies.py` script to test all modules, subworkflows, and workflows that depend on a given component (e.g., BBDUK, or TAXONOMY)
- Modified EXTRACT_VIRAL_READS_SHORT and EXTRACT_VIRAL_READS_ONT to handle empty input files (which previously caused failures).

# v2.8.3.2
- Modified FASTQ_LABELED to use fixed cpus and memory, and added `--memory` parameter to make full use of available memory.
- Added pass/fail test for FASTQC_LABELED.
- Removed unused QC processes.
- Added rank-raised taxids to viral taxonomy DB output by INDEX workflow.

# v2.8.3.1
- Added `expected-outputs-run.txt` file containing list of expected output files for the `RUN` workflow (excluding BLAST validation).
- Minor updates to logging filenames.

# v2.8.3.0
- **Lowered Bracken read threshold for taxon classification**

# v2.8.2.0
- **Increased runtime Bowtie2 score threshold for viral read identification**
- Updated Github Actions to use NAO secrets to access buckets containing test data
- Removed generate-samplesheet.py, as functionality has moved to internal mgs-metadata repo
- Added ability to set BLAST parameters `qcov_hsp_perc` and `perc_id`
- Added MINIMAP2 classification of ONT reads to PROFILE subworkflow
- Replaced boolean `params.ont` with string `params.platform` and added platform checking to LOAD_SAMPLESHEET
- Fixed bug in running RUN_VALIDATION workflow with a FASTQ file

# v2.8.1.2
- Made Cutadapt mismatch rate parameter configurable
- Fixed issues with BLAST bitscore filtering
- Increased memory allocation for EXTRACT_VIRAL_HITS_TO_FASTQ
- Implemented version compatibility checking between pipeline and index
- Added ONT virus identification support:
    - Created new EXTRACT_VIRAL_READS_ONT subworkflow for processing ONT reads
    - Renamed original EXTRACT_VIRAL_READS workflow to EXTRACT_VIRAL_READS_SHORT to differentiate from ONT processing
    - Added non-streaming version of MINIMAP2 alignment process
    - Added new modules for ONT-specific processing:
        - MASK_FASTQ_READS for masking low complexity regions in reads
        - EXTRACT_SHARED_FASTQ_READS for extracting reads shared between FASTQ files
        - PROCESS_VIRAL_MINIMAP2_SAM for adding reference taxids and clean read information
    - Edited FILTLONG to accept customizable parameters (min_length, max_length, min_mean_q)
    - Added new low-complexity fastq test file.

# v2.8.1.1
- Modified Kraken2 DB handling in index workflow to avoid staging
- Updated defaults in index configs

# v2.8.1.0
- Added downstream duplicate marking functionality via new DOWNSTREAM workflow
    - Fixed JOIN_TSVS to correctly handle many-to-one joins
    - Added strict join mode to JOIN_TSVS
    - Altered PROCESS_VIRAL_BOWTIE2_SAM to make ordering of genome IDs for split alignments predictable (necessary for downstream duplicate marking)
- Updated ANNOTATE_VIRUS_INFECTION to better handle taxa that are missing from Virus-Host DB, and added corresponding tests and documentation.
- Began implementing pipeline components for analyzing ONT data:
    - Added generation of minimap2 indices to INDEX workflow (human, viral, contaminant, and ribosomal).
    - Added LSU and SSU tags to respective small and large ribosomal subunit genomes in the composite ribosomal reference fasta.
    - Added MINIMAP2_INDEX and MINIMAP2 processes for indexing reference genomes and aligning reads to them.
- Added documentation on running the pipeline reproducibly
- Fixed some local unit tests

# v2.8.0.0
- Major changes to many parts of the pipeline as part of a general performance overhaul
    - Modified most processes in the RUN and RUN_VALIDATION workflows to stream data in and out rather than reading whole files
    - As part of the previous change, modified most processes in the RUN and RUN_VALIDATION workflows to work with interleaved rather than paired sequence data
    - Modified BLASTN filtering to take into account bitscore ratio versus best hit for each query
    - Replaced many specific tabular manipulation processes with basic operations: JOIN_TSVS, CONCATENATE_TSVS, ADD_FIXED_COLUMN, etc
    - Removed grouping and group-dependent functionality (in particular, deduplication and clade counting); entire pipeline now operates on a per-sample basis
    - Added unit tests for many processes and workflows
    - Added configurable seeding for testing non-deterministic processes via `params.random_seed`
    - Made Bracken read threshold configurable via `params.bracken_threshold`
    - Removed numerous orphaned modules and processes
- Large changes to outputs:
    - Main output directory no longer contains FASTA files for viral hits (interleaved FASTQ file now saved to intermediates)
    - Clade counts are no longer produced
    - QC and BLAST outputs now show statistics for interleaved files rather than showing forward and reverse reads separately
    - Added new intermediate outputs, including unfiltered viral hits and interleaved FASTQ from EXTRACT_VIRAL_READS
    - Viral hits TSV moved from `virus_hits_db.tsv.gz` to `virus_hits_filtered.tsv.gz`
    - Numerous changes to column names in viral hits TSV, mainly to improve clarity
- Minor changes and fixes:
    - Updated mislabeled processes
    - Fixed bug where multiqc doesn't output sequence length stats if all sequences are the same length
    - Unzipped files in `test-data` directory
    - Added new script, `bin/run_parallel_test.sh`, that allows users to run nf-test tests locally in parallel
    - Assorted updates to documentation
    - Removed some defaults from config files
    - Fixed mislabeled parameter in RUN_VALIDATION workflow

# v2.7.0.3
- Fixing link to configuration file in `README.md`

# v2.7.0.2
- Updated `pipeline-version.txt`

# v2.7.0.1
- Fixed index-related issues from v2.7.0.0:
    - Updated `EXTRACT_VIRAL_READS` to expect updated path to viral genome DB
    - Added `adapters` param to the index config file used to run our tests
    - Updated `RUN` and `RUN_VALIDATION` tests to use up-to-date test index (location: `s3://nao-testing/index/20250130`)

# v2.7.0.0
- Implemented masking of viral genome reference in index workflow with MASK_GENOME_FASTA to remove adapter, low-entropy and repeat sequences.
- Removed TRIMMOMATIC and BBMAP from EXTRACT_VIRAL_READS.
- Restructured subworkflows to take advantage of new viral genome masking:
    - Split PROFILE workflow into SUBSET_TRIM, RUN_QC, and PROFILE subworkflows
    - Moved FASTP read cleaning downstream of BBDUK_HITS (in EXTRACT_VIRAL_READS) and subsetting (in SUBSET_TRIM)
    - Moved FASTQC and MultiQC to after subsetting (in RUN_QC)
    - Removed RAW, CLEAN, and PROCESS_OUTPUT subworkflows
    - Added COUNT_TOTAL_READS subworkflow to count the total number of reads in each sample.
- Replace generate_samplesheet.sh with generate_samplesheet.py
- Fixed bug in extractUnconcReadID that would cause the pipeline to fail if it contained the string 'YT' in the read id.
- Remove `params.quality_encoding` as it was used only by TRIMMOMATIC
- Added length distribution information to QC output
- **Renamed QC output files to reflect the fact that they now only contain QC information on a subset of reads (e.g. `qc_basic_stats.tsv.gz` -> `subset_qc_basic_stats.tsv.gz`)**
- **New QC output files: `read_counts.tsv.gz`, `subset_qc_length_stats.tsv.gz`**

# v2.6.0.0
- Updated version to reflect the new versioning scheme, which is described in `docs/version_schema.md`.

# v2.5.4
- Fixed fatal bug in `configs/run_validation.config` that prevents users from running the `RUN_VALIDATION` workflow.

# v2.5.3
- Added new LOAD_SAMPLESHEET subworkflow to centralize samplesheet processing
- Updated tags to prevent inappropriate S3 auto-cleanup
- Testing infrastructure
  - Split up the tests in `End-to-end MGS workflow test` so that they can be run in parallel on Github Actions.
  - Implemented an end-to-end test that checks if the RUN workflow produces the correct output. The correct output for the test has been saved in `test-data/gold-standard-results` so that the user can diff the output of their test with the correct output to check where their pipeline might be failing.
- Began development of single-end read processing (still in progress)
    - Restructured RAW, CLEAN, QC, TAXONOMY, and PROFILE workflows to handle both single-end and paired-end reads
    - Added new FASTP_SINGLE, TRUNCATE_CONCAT_SINGLE, BBDUK_SINGLE, CONCAT_GROUP_SINGLE, SUBSET_READS_SINGLE and SUBSET_READS_SINGLE_TARGET processes to handle single-end reads
    - Created separate end-to-end test workflow for single-end processing (which will be removed once single-end processing is fully integrated)
    - Modified samplesheet handling to support both single-end and paired-end data
    - Updated generate_samplesheet.sh to handle single-end data with --single_end flag
    - Added read_type.config to handle single-end vs paired-end settings (set automatically based on samplesheet format)
    - Created run_dev_se.config and run_dev_se.nf for single-end development testing (which will be removed once single-end processing is fully integrated)
    - Added single-end samplesheet to test-data

# v2.5.2
- Changes to default read filtering:
    - Relaxed FASTP quality filtering (`--cut_mean_quality` and `--average_qual` reduced from 25 to 20).
    - Relaxed BBDUK viral filtering (switched from 3 21-mers to 1 24-mer).
- Overhauled BLAST validation functionality:
    - BLAST now runs on forward and reverse reads independently
    - BLAST output filtering no longer assumes specific filename suffixes
    - Paired BLAST output includes more information
    - RUN_VALIDATION can now directly take in FASTA files instead of a virus read DB
    - Fixed issues with publishing BLAST output under new Nextflow version
- Implemented nf-test for end-to-end testing of pipeline functionality
    - Implemented test suite in `tests/main.nf.test`
    - Reconfigured INDEX workflow to enable generation of miniature index directories for testing
    - Added Github Actions workflow in `.github/workflows/end-to-end.yml`
    - Pull requests will now fail if any of INDEX, RUN, or RUN_VALIDATION crashes when run on test data.
    - Generated first version of new, curated test dataset for testing RUN workflow. Samplesheet and config file are available in `test-data`. The previous test dataset in `test` has been removed.
- Implemented S3 auto-cleanup:
    - Added tags to published files to facilitate S3 auto-cleanup
    - Added S3 lifecycle configuration file to `ref`, along with a script in `bin` to add it to an S3 bucket
- Minor changes
    - Added logic to check if `grouping` variable in `nextflow.config` matches the input samplesheet, if it doesn't, the code throws an error.
    - Externalized resource specifications to `resources.config`, removing hardcoded CPU/memory values
    - Renamed `index-params.json` to `params-index.json` to avoid clash with Github Actions
    - Removed redundant subsetting statement from TAXONOMY workflow.
    - Added --group_across_illumina_lanes option to generate_samplesheet

# v2.5.1
- Enabled extraction of BBDuk-subset putatively-host-viral raw reads for downstream chimera detection.
- Added back viral read fields accidentally being discarded by COLLAPSE_VIRUS_READS.

# v2.5.0
- Reintroduced user-specified sample grouping and concatenation (e.g. across sequencing lanes) for deduplication in PROFILE and EXTRACT_VIRAL_READS.
- Generalised pipeline to detect viruses infecting arbitrary host taxa (not just human-infecting viruses) as specified by `ref/host-taxa.tsv` and config parameters.
- Configured index workflow to enable hard-exclusion of specific virus taxa (primarily phages) from being marked as infecting ost taxa of interest.
- Updated pipeline output code to match changes made in latest Nextflow update (24.10.0).
- Created a new script `bin/analyze-pipeline.py` to analyze pipeline structure and identify unused workflows and modules.
- Cleaned up unused workflows and modules made obsolete in this and previous updates.
- Moved module scripts from `bin` to module directories.
- Modified trace filepath to be predictable across runs.
- Removed addParams calls when importing dependencies (deprecated in latest Nextflow update).
- Switched from nt to core_nt for BLAST validation.
- Reconfigured QC subworkflow to run FASTQC and MultiQC on each pair of input files separately (fixes bug arising from allowing arbitrary filenames for forward and reverse read files).

# v2.4.0
- Created a new output directory where we put log files called `logging`.
- Added the trace file from Nextflow to the `logging` directory which can be used for understanding cpu, memory usage, and other infromation like runtime. After running the pipeline, `plot-timeline-script.R` can be used to generate a useful summary plot of the runtime for each process in the pipeline.
- Removed CONCAT_GZIPPED.
- Replaced the sample input format with something more similar to nf-core, called `samplesheet.csv`. This new input file can be generated using the script `generate_samplesheet.sh`.
- Now run deduplication on paired-ends reads using clumpify in the taxonomic workflow.
- Fragment length analysis and deduplication analysis.
  - BBtools: Extract the fragment length as well as the number of duplicates from the taxonomic workflow and add them to the `hv_hits_putative_collapsed.tsv.gz`.
  - Bowtie2: Conduct a duplication analysis on the aligned reads, then add the number of duplicates and fragment length to the `hv_hits_putative_collapsed.tsv.gz`.

# v2.3.3
- Added validation workflow for post-hoc BLAST validation of putative HV reads.

# v2.3.2
- Fixed subsetReads to run on all reads when the number of reads per sample is below the set threshold.

# v2.3.1

- Clarifications to documentation (in README and elsewhere)
- Re-added "joined" status marker to reads output by join_fastq.py

# v2.3.0
- Restructured run workflow to improve computational efficiency, especially on large datasets
    - Added preliminary BBDuk masking step to HV identification phase
    - Added read subsampling to profiling phase
    - Deleted ribodepletion and deduplication from preprocessing phase
    - Added riboseparation to profiling phase
    - Restructured profiling phase output
    - Added `addcounts` and `passes` flags to deduplication in HV identification phase
- Parallelized key bottlenecks in index workflow
- Added custom suffix specification for raw read files
- Assorted bug fixes

# v2.2.1
- Added specific container versions to `containers.config`
- Added version & time tracking to workflows
- Added index reference files (params, version) to run output
- Minor changes to default config files

# v2.2.0
- Major refactor
- Start of changelog
