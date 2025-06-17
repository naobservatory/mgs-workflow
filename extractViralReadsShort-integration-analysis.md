# extractViralReadsShort Integration Analysis

## Overview

Analysis of compatibility between the new `hits_final` output from extractViralReadsShort (with LCA integration) and the downstream workflow requirements.

## Format Comparison

### Old Format (42 columns)
**Source:** Historical `virus_hits_final.tsv` files
- 1 seq_id column
- 33 Bowtie2 alignment columns
- 7 Kraken2 classification columns (unused by downstream)
- 1 bbmerge_frag_length column (unused by downstream)
- 1 sample column

### New Format (59 columns)
**Source:** extractViralReadsShort with JOIN_TSVS + FILTER_PRIMARY_ALIGNMENTS
- 1 seq_id column ✅
- 33 Bowtie2 alignment columns ✅ (plus `is_secondary`)
- 25 LCA classification columns (new)
- 1 sample column ✅

## Compatibility Analysis

### ✅ Compatible Elements
- **All required Bowtie2 columns present** - Downstream duplicate marking and validation will work
- **Sample column present** - Group joining will work
- **Core identifiers present** - seq_id for read tracking
- **Additional LCA columns** - Will be ignored by downstream (no negative impact)

### ❌ Potentially Incompatible Elements
- **Missing Kraken2 columns** - But analysis shows these are unused by downstream workflow
- **Missing bbmerge_frag_length** - But analysis shows this is unused by downstream workflow
- **Extra `is_secondary` column** - Should be ignored by downstream

### 🔍 Key Finding
**The missing columns (Kraken2, BBMerge) are historical artifacts that are NOT used by the downstream workflow.** They appear in test data but are never referenced by processing code.

## Integration Steps Completed

### 1. extractViralReadsShort Modifications ✅
- Added JOIN_TSVS integration to combine LCA and Bowtie2 data
- Added FILTER_PRIMARY_ALIGNMENTS to keep only primary alignments (`is_secondary=False`)
- Modified workflow to sort data before joining (required by JOIN_TSVS)
- Updated process to handle gzipped input/output

### 2. Column Structure ✅
The new workflow produces this data flow:
1. **PROCESS_VIRAL_BOWTIE2_SAM** → Bowtie2 alignment data (35 columns including `is_secondary`)
2. **LCA_TSV** → Taxonomic classification data (25 columns)
3. **JOIN_TSVS** → Combined data (58 columns)
4. **FILTER_PRIMARY_ALIGNMENTS** → Primary alignments only (58 columns)
5. **ADD_SAMPLE_COLUMN** → Final output (59 columns)

### 3. Test Infrastructure ✅
- Created comprehensive test suite for FILTER_PRIMARY_ALIGNMENTS
- Fixed gzip handling in module and tests
- All tests passing

## Recommendation

**Proceed with direct compatibility testing** rather than building compatibility adapters:

1. **Create test TSV** with new 59-column format
2. **Test downstream workflow** with new format
3. **Expect success** based on column analysis

The downstream workflow should work with the new format because:
- All actually-used columns are present
- Missing columns are unused historical artifacts
- Extra columns will be ignored

## Future Considerations

### Column Subsetting Opportunity
Based on downstream requirements analysis, `filterPrimaryAlignments` could optionally subset to only required columns:
- Reduces file size and processing overhead
- Maintains compatibility
- Focuses on data that's actually used

### Test Data Updates
Once compatibility is confirmed:
- Update downstream test data to use new format
- Remove unused Kraken2/BBMerge columns from test expectations
- Simplify test data maintenance

## Implementation Status

- ✅ extractViralReadsShort integration complete
- ✅ FILTER_PRIMARY_ALIGNMENTS module complete and tested
- ✅ Column analysis complete
- 🔄 Direct compatibility testing pending

## Last Updated
2025-06-17 - Post-integration analysis