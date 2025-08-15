import pytest
import tempfile
import os
from validate_grouping import validate_grouping

def write_tsv(path, header, rows):
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

def test_validate_grouping_superset_and_fail():
    with tempfile.TemporaryDirectory() as tmp:
        grouping_path = os.path.join(tmp, 'grouping.tsv')
        virus_hits_path = os.path.join(tmp, 'virus_hits.tsv')
        validated_output = os.path.join(tmp, 'validated_grouping.tsv')
        samples_without_vv_hits = os.path.join(tmp, 'samples_without_vv_hits.tsv')

        # Case 1: Fail when grouping is not a superset
        write_tsv(grouping_path, ['group', 'sample'], [['g1', 'S1']])
        write_tsv(virus_hits_path, ['sample'], [['S1'], ['S2']])
        with pytest.raises(ValueError):
            validate_grouping(grouping_path, virus_hits_path, samples_without_vv_hits, validated_output)

        # Case 2: Pass when grouping is superset, check outputs
        write_tsv(grouping_path, ['group', 'sample'], [['g1', 'S1'], ['g2', 'S2']])
        write_tsv(virus_hits_path, ['sample'], [['S1']])
        validate_grouping(grouping_path, virus_hits_path, samples_without_vv_hits, validated_output)

        with open(validated_output) as f:
            lines = [l.strip() for l in f if l.strip()]
        assert lines == ['group\tsample', 'g1\tS1']

        with open(samples_without_vv_hits) as f:
            unused_lines = [l.strip() for l in f if l.strip()]
        assert unused_lines == ['S2']

def test_validate_grouping_empty_files():
    with tempfile.TemporaryDirectory() as tmp:
        grouping_path = os.path.join(tmp, 'grouping.tsv')
        virus_hits_path = os.path.join(tmp, 'virus_hits.tsv')
        validated_output = os.path.join(tmp, 'validated_grouping.tsv')
        samples_without_vv_hits = os.path.join(tmp, 'samples_without_vv_hits.tsv')

        # Case: Empty virus hits file with non-empty grouping
        write_tsv(grouping_path, ['group', 'sample'], [['g1', 'S1'], ['g2', 'S2']])
        write_tsv(virus_hits_path, ['sample'], [])  # Empty data, just header
        validate_grouping(grouping_path, virus_hits_path, samples_without_vv_hits, validated_output)

        # Validated output should only have header, no data
        with open(validated_output) as f:
            validated_lines = [l.strip() for l in f if l.strip()]
        assert validated_lines == ['group\tsample']

        # All grouping samples should be in samples_without_vv_hits
        with open(samples_without_vv_hits) as f:
            unused_lines = [l.strip() for l in f if l.strip()]
        assert set(unused_lines) == {'S1', 'S2'}

def test_validate_grouping_perfect_match():
    """Test that when all samples in grouping are present in virus hits, 
    the validated output contains all samples from grouping."""
    with tempfile.TemporaryDirectory() as tmp:
        grouping_path = os.path.join(tmp, 'grouping.tsv')
        virus_hits_path = os.path.join(tmp, 'virus_hits.tsv')
        validated_output = os.path.join(tmp, 'validated_grouping.tsv')
        samples_without_vv_hits = os.path.join(tmp, 'samples_without_vv_hits.tsv')

        # Case: All samples in grouping have virus hits
        write_tsv(grouping_path, ['group', 'sample'], [['g1', 'S1'], ['g1', 'S2'], ['g2', 'S3']])
        write_tsv(virus_hits_path, ['sample'], [['S1'], ['S2'], ['S3']])
        validate_grouping(grouping_path, virus_hits_path, samples_without_vv_hits, validated_output)

        # Read validated output
        with open(validated_output) as f:
            validated_lines = [l.strip() for l in f if l.strip()]
        
        # Should have header plus all 3 samples (sorted by sample name)
        expected_lines = ['group\tsample', 'g1\tS1', 'g1\tS2', 'g2\tS3']
        assert validated_lines == expected_lines

        # No samples should be without viral hits
        with open(samples_without_vv_hits) as f:
            unused_lines = [l.strip() for l in f if l.strip()]
        assert unused_lines == []

