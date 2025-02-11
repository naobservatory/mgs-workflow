# Nucleic Acid Observatory Viral Metagenomics Pipeline

This Nextflow pipeline is designed to process metagenomic sequencing data, characterize overall taxonomic composition, and identify and quantify reads mapping to viruses infecting certain host taxa of interest. It was developed as part of the [Nucleic Acid Observatory](https://naobservatory.org/) project.

The pipeline currently consists of three workflows:

- [`index`](./docs/index.md): Creates indices and reference files used by the `run` and `run_validation` workflows[^1].
- [`run`](./docs/run.md): Performs the main analysis, including QC, viral identification, taxonomic profiling, and optional BLAST validation.
- `run_validation`: Performs part of the `run` workflow dedicated to validation of taxonomic classification with BLAST[^2].

[^1]: The `index` workflow is intended to be run first, after which many instantiations of the `run` workflow can use the same index output files. 
[^2]: The `run_validation` workflow is intended to be run after the `run` workflow if the optional BLAST validation was not selected during the `run` workflow. Typically, this workflow is run on a subset of the host viral reads identified in the `run` workflow, to evaluate the sensitivity and specificity of the viral identification process.

## Documentation

- **Installation and usage:**
    - [Installation instructions](docs/installation.md)
    - [AWS Batch setup](docs/batch.md)
    - [Usage instructions](docs/usage.md)
    - [Troubleshooting](docs/troubleshooting.md)
- **Workflow details:**
    - [Index workflow](docs/index.md)
    - [Run workflow](docs/run.md)
- **Configuration and output:**
    - [Configuration files](docs/config.md)
    - [Pipeline outputs](docs/output.md)
- **Other:**
    - [Version Schema](docs/versioning.md)
