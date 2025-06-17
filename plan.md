# Plan: Join LCA TSV with Bowtie2 Viral Reads TSV

## Goal
Create a process to add alignment score and position columns from the Bowtie2 viral reads TSV to the LCA TSV output, keeping only primary alignment data.

## Approach: Modify Upstream + Modular Join + Filter

### Three-Step Process:
1. **Modify `process_viral_bowtie2_sam.py`** to include `is_secondary` column in output TSV
2. **Use existing `JOIN_TSVS` process** to join LCA TSV with Bowtie2 TSV on `seq_id`
3. **Create new filtering process** that processes the joined TSV to keep only primary alignments

### Implementation Tasks:

#### Task 1: Add Secondary Alignment Flag to Bowtie2 Script
- Modify `process_viral_bowtie2_sam.py` to include `is_secondary` column in output headers
- Update output functions to include the `is_secondary` flag in the final TSV
- Ensure compatibility with both paired and unpaired output formats
- Test modifications don't break existing functionality

#### Task 2: Create Primary Alignment Filter Script
- Write Python script (`filter_primary_alignments.py`) that:
  - Reads joined TSV from `JOIN_TSVS` output
  - Filters for rows where `is_secondary == False` (primary alignments only)
  - Outputs clean TSV with LCA data + primary alignment info
- Include proper error handling and logging
- Add type hints and comprehensive docstrings

#### Task 3: Create Nextflow Module
- Create new module `modules/local/filterPrimaryAlignments/main.nf`
- Wrap filter script in Nextflow process
- Define proper input/output channels
- Add resource requirements and labels
- Include input validation

#### Task 4: Integration and Testing
- Create comprehensive nf-test suite
- Test with sample data including edge cases
- Document the new workflow pattern
- Update any relevant documentation

### Benefits:
- ✅ Reuses proven `JOIN_TSVS` functionality  
- ✅ Single responsibility for new process (filtering only)
- ✅ Modular and maintainable design
- ✅ Follows existing pipeline patterns
- ✅ Easy to test individual components

### Files to Create:
- `modules/local/filterPrimaryAlignments/main.nf`
- `modules/local/filterPrimaryAlignments/resources/usr/bin/filter_primary_alignments.py`
- `tests/modules/local/filterPrimaryAlignments/main.nf.test`
- Test data files in `test-data/toy-data/`