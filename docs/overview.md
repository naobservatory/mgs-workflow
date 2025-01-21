# Overview

The pipeline currently consists of three workflows:

- `run`: Performs the main analysis, including viral identification, taxonomic profiling, and optional BLAST validation.
- `index`: Creates indices and reference files used by the `run` workflow.[^1]
- `run_validation`: Performs BLAST validation on any set of reads to verify their taxonomic classification.[^2]

[^1]: The `index` workflow is intended to be run first, after which many instantiations of the `run` workflow can use the same index output files. 
[^2]: The `run_validation` workflow is intended to be run after the `run` workflow if the optional BLAST validation was not selected during the `run` workflow. Typically, this workflow is run on a subset of the host viral reads identified in the `run` workflow, to evaluate the sensitivity and specificity of the viral identification process.

The `run` workflow consists of four parts:

- **Viral identification**: Reads from viruses infecting specific host taxa of interest (default: vertebrate viruses) are sensititvely and specifically identified using a multi-step pipeline based around k-mer filtering, adapter trimming, read alignment, and taxonomic classification.
- **Taxonomic profile**: A random subset of reads (default: 1M/sample) undergo adapter trimming, ribosomal classification, and taxonomic classification.
- **QC**: The total number of reads are recorded, then a subset of reads (default 1M/sample) undergo quality assessment, both before and after adapter trimming.
- **(Optional) BLAST validation**: Putative host viral reads from part 1 (either all reads or a subset) are checked aginst `core_nt` to evaluate the sensitivity and specificity of the viral identification process.

The following diagram provides a high-level overview of the `run` workflow (each blue box is a subworkflow):

```mermaid
---
title: RUN WORKFLOW
config:
  layout: horizontal
  look: handDrawn
  theme: defaultr
---
flowchart LR
A[Input Files] --> |Raw reads|B[LOAD_SAMPLESHEET]
B --> C[COUNT_TOTAL_READS] & D[EXTRACT_VIRAL_READS] & E[SUBSET_TRIM]
D -.-> |*Optional*| F[BLAST_VIRAL]
E --> |Trim and subset reads| G[RUN_QC] & H[PROFILE]
E --> |Subset reads|G[RUN_QC]
%% Adjust layout by placing subgraphs in specific order
subgraph "Viral identification"
D
end
subgraph "Taxonomic Profile"
H
end
subgraph "QC"
C
G
end
subgraph "BLAST Validation"
F
end
```