# Downstream Workflow Column Requirements

## Overview

This document outlines the actual column dependencies in the downstream workflow based on code analysis. This information is crucial for understanding what columns are truly required versus historical artifacts.

## Actually Required Columns

### For Duplicate Marking (`markViralDuplicates`)
- `seq_id` - Primary read identifier for grouping duplicates
- `aligner_genome_id_all` - Groups reads by genome assignment
- `aligner_ref_start` / `aligner_ref_start_rev` - Alignment positions for coordinate-based duplicate detection
- `query_qual` / `query_qual_rev` - Quality scores for selecting exemplar reads

### For Taxonomic Validation (`validateViralAssignments`)
- `aligner_taxid` - Used for taxonomic distance calculations and database joins
- `query_seq` / `query_seq_rev` - Read sequences for FASTQ extraction and clustering
- `sample` - Required for joining with groups TSV

### For Workflow Infrastructure
- `sample` - Essential for joining hits data with sample groups metadata

## Unused Historical Columns

These columns appear in current test data but are **never referenced** by downstream workflow code:

### Kraken2 Columns (Unused)
- `kraken2_classified` - Only used by `processKrakenViral` module (not part of downstream)
- `kraken2_assigned_name` - Only used by `processKrakenViral` module (not part of downstream)
- `kraken2_assigned_taxid` - Only used by `processKrakenViral` module (not part of downstream)
- `kraken2_assigned_host_virus` - Only used by `processKrakenViral` module (not part of downstream)
- `kraken2_length` - Only used by `processKrakenViral` module (not part of downstream)
- `kraken2_encoded_hits` - Only used by `processKrakenViral` module (not part of downstream)

### BBMerge Columns (Unused)
- `bbmerge_frag_length` - Only used by `summarizeBBMerge` module (not part of downstream)

## Analysis Method

This analysis was conducted by:
1. Examining all downstream workflow modules and subworkflows
2. Searching Python scripts for column name references
3. Tracing data flow through processing logic
4. Identifying actual column dependencies vs. test data artifacts

## Implications for Future Development

### For Column Subsetting
When implementing column subsetting in `filterPrimaryAlignments`, these are the minimum required columns for downstream compatibility:
- Core identifier: `seq_id`, `sample`
- Duplicate marking: `aligner_genome_id_all`, `aligner_ref_start`, `aligner_ref_start_rev`, `query_qual`, `query_qual_rev`
- Validation: `aligner_taxid`, `query_seq`, `query_seq_rev`

### For Test Data
- Kraken2 and BBMerge columns can be removed from future test data
- Focus test data on columns that are actually processed
- New LCA-based output format should be compatible with downstream workflow

## Last Updated
2025-06-17 - Initial analysis based on current downstream workflow implementation