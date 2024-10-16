# v2.4.0
- Created a new output directory where we put log files called `logging`. 
- Added the trace file from Nextflow to the `logging` directory which can be used for understanding cpu, memory usage, and other infromation like runtime. After running the pipeline, `plot-timeline-script.R` can be used to generate a useful summary plot of the runtime for each process in the pipeline.
- Removed CONCAT_GZIPPED.
- Replaced the sample input format with something more similar to nf-core, called `samplesheet.csv`. This new input file can be generated using the script `generate_samplesheet.sh`.
- Now run deduplication on paired-ends reads using clumpify in the taxonomic workflow.
- Fragment length analysis and deduplication analysis.
  - BBtools: Extract the fragment length as well as the number of duplicates from the taxonomic workflow and add them to the `hv_hits_putative_collapsed.tsv.gz`.
  - Bowtie2: Conduct a duplication analysis on the aligned reads, then add the number of duplicates and fragment length to the `hv_hits_putative_collapsed.tsv.gz`.

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
