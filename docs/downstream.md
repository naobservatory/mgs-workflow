# DOWNSTREAM WORKFLOW

This page describes the structure and function of the `DOWNSTREAM` workflow. This workflow is responsible for downstream analysis of the outputs of the [`RUN` workflow](./run.md), particularly in cases that require comparisons across reads and/or samples[^comp]. Currently, this workflow only performs a single analysis, namely identification and marking of duplicate reads based on their Bowtie2 alignment results.

[^comp]: These are kept to a minimum in the `RUN` workflow to minimize memory demands and maximize parallelization.

## Workflow structure

```mermaid
---
title: DOWNSTREAM WORKFLOW
config:
  layout: horizontal
---
flowchart LR
A(Viral hits tables) & B(Grouping information) --> C[LOAD_DOWNSTREAM_DATA]
C --> D[PREPARE_GROUP_TSVS]
D --> E[MARK_VIRAL_DUPLICATES]
subgraph "Partition and concatenate TSVs by sample group"
D
end
subgraph "Duplicate annotation"
E
end
style A fill:#fff,stroke:#000
style B fill:#fff,stroke:#000
```

## Subworkflows

### Load data into channels (`LOAD_DOWNSTREAM_DATA`)

This subworkflow takes in an input file specifying (1) paths to one or more viral hits tables produced by the `RUN`workflow, and (2) paths to corresponding TSV files specifying the sample groupings to be used for duplicate annotation (see [below](#usage) for more information on this input file). The subworkflow validates that this input file has the required structure, then extracts the filepaths into a channel with the structure expected by the rest of the workflow. (No diagram is provided for this subworkflow.)

### Partition into groups for duplicate annotation (`PREPARE_GROUP_TSVS`)

This subworkflow takes in the channel output by `LOAD_DOWNSTREAM_DATA`, adds sample grouping information to the viral hits tables, then partitions each viral hits table into a separate TSV per sample group. Partitions from different hits tables with matching group annotations are then concatenated together, enabling duplicate annotation across different pipeline runs (e.g. from different data deliveries) as specified by the user.

```mermaid
---
title: PREPARE_GROUP_TSVS
config:
  layout: horizontal
---
flowchart LR
A(Viral hits tables) --> B[SORT_TSV]
C(Grouping information) --> D[SORT_TSV]
B & D --> E[JOIN_TSVS]
E --> F[PARTITION_TSVS]
F --> G[CONCATENATE_TSVS_LABELED]
G --> H(Partitioned sample group TSVs)
subgraph X [Add grouping information to hits tables]
B
D
E
end
subgraph Y [Partition and concatenate TSVs by sample group]
F
G
end
style A fill:#fff,stroke:#000
style C fill:#fff,stroke:#000
style H fill:#000,color:#fff,stroke:#000
```

### Annotate alignment duplicates (`MARK_VIRAL_DUPLICATES`)

This subworkflow takes in partitioned hits tables from `PREPARE_GROUP_TSVS`, then identifies duplicate reads on the basis of their assigned genome ID and alignment coordinates, as determined by Bowtie2 in the `RUN` workflow. In order to be considered duplicates, two read pairs must be mapped to the same genome ID by Bowtie2[^gid], with terminal alignment coordinates that are within a user-specified distance of each other (default 1 nt) at both ends. This fuzzy matching allows for the identification of duplicate reads in the presence of small read errors, alignment errors or overzealous adapter trimming.

For each group of reads identified as duplicates, the algorithm selects the read pair with the highest average quality score to act as the "exemplar" of the group. Each read in the group is annotated with this examplar to identify its duplicate group[^exemplar], enabling downstream deduplication or other duplicate analyses if needed. In addition to an annotated hits TSV containing an additional column for exemplar IDs, the subworkflow also returns a summary TSV giving the number of reads mapped to a given exemplar ID, as well as the fraction of read pairs in the group that are pairwise duplicates[^pairwise].

[^exemplar]: A read with no duplicates will be annotated with itself as the exemplar.
[^pairwise]: Because of the fuzzy matching used to identify duplicates, it is possible for duplicate annotation to be intransitive: i.e. read A is a duplicate of read B, and read B is a duplicate of read C, but read A is not a duplicate of read C. As currently implemented, the algorithm will group a read into a duplicate group if it matches any single read already in that duplicate group, potentially leading to the grouping of reads that would not be considered duplicates of each other in isolation. The reporting of the pairwise duplicate statistic in the summary file allows for quantification of this phenomenon, and potential adjustment of parameters if too high a fraction of non-matching reads are being grouped together in this way.

```mermaid
---
title: MARK_VIRAL_DUPLICATES
config:
  layout: horizontal
---
flowchart LR
A("Partitioned sample group TSVs <br> (PREPARE_GROUP_TSVS)") --> B[MARK ALIGNMENT DUPLICATES]
B --> C[SORT_TSV]
B --> D[SORT_TSV]
C --> E(Annotated hits TSVs)
D --> F(Summary TSVs)
style A fill:#fff,stroke:#000
style E fill:#000,color:#fff,stroke:#000
style F fill:#000,color:#fff,stroke:#000
```

## Usage

> [!IMPORTANT]
> As with the [`RUN` workflow](./usage.md), before following the instructions in this section, make sure you have followed the [installation and setup instructions](./installation.md).

To run the `DOWNSTREAM` workflow, you need:

1. One or more accessible **viral hits tables** produced by the `RUN` workflow. These are [typically saved](./output.md#viral-identification)  in the `RUN` workflow's output directory under `results/virus_hits_filtered.tsv.gz`.
2. For each hit table, an accessible **grouping TSV**, containing the following columns in the specified order:
    - `sample`: Sample ID (must include one row for every value of this column in the viral hits table)
    - `group`: Group IDs to use for duplicate annotatation
3. An accessible **input file CSV** mapping viral hits tables to grouping TSVs, containing the following columns in the specified order:
    - `label`: Arbitrary string label to use for each viral hits table
    - `hits_tsv`: Path to the viral hits table
    - `groups_tsv`: Path to the corresponding grouping TSV
4. A **config file** in a clean launch directory, pointing to:
    - The pipeline mode (`params.mode = "downstream"`);
    - The input file (`params.input_file`);
    - The base directory in which to put the working and output directories (`params.base_dir`);
    - The permitted deviation when identifying alignment duplicates (`params.aln_dup_deviation`).

> [!NOTE]
> Currently, the input file and grouping TSV must be generated manually. We intend to implement programmatic generation of these files in the future.

> [!TIP]
> We recommend starting each pipeline run in a clean launch directory, containing only your input file and config file.

Given these input files, you must choose a run profile as described [here](./usage.md#2-choosing-a-profile). You can then run the pipeline as follows:

```
nextflow run -resume -profile <PROFILE> <PATH/TO/PIPELINE/DIR>
```

where `<PATH/TO/PIPELINE/DIR>` specifies the path to the directory containing the pipeline files from this repository (in particular, `main.nf`) from the launch directory.

Once the pipeline has finished, output and logging files will be available in the `output` subdirectory of the base directory specified in the config file.

> [!IMPORTANT]
> As with the `RUN` workflow, it's highly recommended to clean up your Nextflow working directory after run completion. You can do this manually or with the `nextflow clean` command.
