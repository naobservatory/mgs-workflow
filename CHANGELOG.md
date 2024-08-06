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
