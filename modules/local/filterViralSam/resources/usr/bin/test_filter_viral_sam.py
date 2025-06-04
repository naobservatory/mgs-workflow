#!/usr/bin/env python3

import pytest
import tempfile
import os
import sys
from pathlib import Path

# Add the script directory to path so we can import the module
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from filter_viral_sam import filter_viral_sam_memory_efficient, SamAlignment

class TestSamAlignment:
    def test_parse_basic_sam_line(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGT\tIIIIIIII\tAS:i:45\tYT:Z:CP"
        alignment = SamAlignment.from_sam_line(line)
        
        assert alignment.qname == 'read1'
        assert alignment.flag == 99
        assert alignment.rname == 'chr1'
        assert alignment.pos == 100
        assert alignment.alignment_score == 45
        assert alignment.pair_status == 'CP'
    
    def test_parse_missing_tags(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGT\tIIIIIIII"
        alignment = SamAlignment.from_sam_line(line)
        
        assert alignment.alignment_score is None
        assert alignment.pair_status is None
    
    def test_parse_up_status(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGT\tIIIIIIII\tAS:i:30\tYT:Z:UP"
        alignment = SamAlignment.from_sam_line(line)
        
        assert alignment.pair_status == 'UP'

class TestCalculateNormalizedScore:
    def test_normal_calculation(self):
        import math
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\t" + "A" * 100 + "\t" + "I" * 100 + "\tAS:i:50"
        alignment = SamAlignment.from_sam_line(line)
        alignment.calculate_normalized_score()
        expected = 50 / math.log(100)
        assert abs(alignment.normalized_score - expected) < 1e-10
    
    def test_none_score(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\t" + "A" * 100 + "\t" + "I" * 100
        alignment = SamAlignment.from_sam_line(line)
        alignment.calculate_normalized_score()
        assert alignment.normalized_score == 0
    
    def test_zero_length(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\t\t\tAS:i:50"
        alignment = SamAlignment.from_sam_line(line)
        alignment.calculate_normalized_score()
        assert alignment.normalized_score == 0

class TestCreateUnmappedMate:
    def test_create_mate_from_read1(self):
        line = "read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGT\tIIIIIIII\tAS:i:45\tYT:Z:UP"
        alignment = SamAlignment.from_sam_line(line)
        
        mate = alignment.create_unmapped_mate()
        
        # Should flip from read1 (64) to read2 (128)
        assert mate.flag & 128  # read2
        assert not (mate.flag & 64)  # not read1
        assert mate.flag & 4  # unmapped
        assert mate.cigar == '*'
    
    def test_create_mate_from_read2(self):
        line = "read1\t147\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGT\tIIIIIIII\tAS:i:45\tYT:Z:UP"
        alignment = SamAlignment.from_sam_line(line)
        
        mate = alignment.create_unmapped_mate()
        
        # Should flip from read2 (128) to read1 (64)
        assert mate.flag & 64  # read1
        assert not (mate.flag & 128)  # not read2
        assert mate.flag & 4  # unmapped

class TestFilterViralSam:
    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def teardown_method(self):
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def create_temp_file(self, content, suffix='.txt'):
        fd, path = tempfile.mkstemp(dir=self.temp_dir, suffix=suffix)
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        return path
    
    def test_filtered_read_keeping(self):
        # SAM content with one filtered read that should be kept
        sam_content = """filtered_read\t99\tchr2\t300\t60\t50M\t=\t400\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:45\tYT:Z:CP
filtered_read\t147\tchr2\t400\t60\t50M\t=\t300\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:43\tYT:Z:CP
read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:CP
read1\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:48\tYT:Z:CP"""
        
        filtered_content = """@filtered_read/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@filtered_read/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII
@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 0.1)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should have both pairs: filtered_read (kept) and read1 (passed score threshold)
        assert len(output_lines) == 4
        assert any('filtered_read' in line for line in output_lines)
        assert any('read1' in line for line in output_lines)
    
    def test_score_threshold_pair_logic(self):
        
        # Should keep both because read2 exceeds threshold
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:20\tYT:Z:CP
read1\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:35\tYT:Z:CP
read2\t99\tchr2\t300\t60\t50M\t=\t400\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:15\tYT:Z:CP
read2\t147\tchr2\t400\t60\t50M\t=\t300\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:18\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII
@read2/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read2/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        # Threshold 10 - read1 pair should pass, read2 pair should fail
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 10.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should only have read1 pair (2 lines)
        assert len(output_lines) == 2
        assert all('read1' in line for line in output_lines)
        assert not any('read2' in line for line in output_lines)
    
    def test_single_read_above_threshold(self):
        # Single read above threshold should be kept
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 10.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        assert len(output_lines) == 1
        assert 'read1' in output_lines[0]
    
    def test_single_read_below_threshold(self):
        # Single read below threshold should be filtered out
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:10\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 10.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        assert len(output_lines) == 0
    
    def test_up_read_missing_mate_primary(self):
        # UP read without mate should get an unmapped mate added
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:UP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 1.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should have 2 lines: original read + created unmapped mate
        assert len(output_lines) == 2

        print(output_lines)
        
        # Parse both lines
        read1_line = [line for line in output_lines if '\t99\t' in line][0]
        read2_line = [line for line in output_lines if '\t167\t' in line][0]
        
        read2_parsed = SamAlignment.from_sam_line(read2_line.strip())
        assert read2_parsed.flag & 4  # unmapped
        assert read2_parsed.cigar == '*'
    
    def test_secondary_alignments_with_score_filtering(self):
        # Secondary alignments (flag >= 256) should be processed separately by reference
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:CP
read1\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:48\tYT:Z:CP
read1\t355\tchr2\t300\t60\t50M\t=\t400\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:45\tYT:Z:CP
read1\t403\tchr2\t400\t60\t50M\t=\t300\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:15\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        # Threshold that allows primary but not secondary chr2
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 17.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should only have primary alignments (2 lines)
        assert len(output_lines) == 2
        primary_flags = [SamAlignment.from_sam_line(line.strip()).flag for line in output_lines]
        assert all(flag < 256 for flag in primary_flags)
    
    def test_empty_input(self):
        sam_content = ""
        filtered_content = ""
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 1.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        assert len(output_lines) == 0
    
    def test_skip_header_lines(self):
        # Headers should be skipped (though they shouldn't be present based on workflow)
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:248956422
read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""        

        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 1.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should only have the read line, headers ignored
        assert len(output_lines) == 1
        assert 'read1' in output_lines[0]
        assert all("@" not in line for line in output_lines)

    def test_multiple_secondary_alignment_pairs(self):
        # Test two pairs of secondary alignments to different references
        # with different score patterns
        sam_content = """read1\t99\tchr1\t100\t60\t50M\t=\t200\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:50\tYT:Z:CP
read1\t147\tchr1\t200\t60\t50M\t=\t100\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:48\tYT:Z:CP
read1\t355\tchr2\t300\t60\t50M\t=\t400\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:45\tYT:Z:CP
read1\t403\tchr2\t400\t60\t50M\t=\t300\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:42\tYT:Z:CP
read1\t355\tchr3\t500\t60\t50M\t=\t600\t150\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tAS:i:25\tYT:Z:CP
read1\t403\tchr3\t600\t60\t50M\t=\t500\t-150\tTGCATGCATGCATGCA\tIIIIIIIIIIIIIIII\tAS:i:20\tYT:Z:CP"""
        
        filtered_content = """@read1/1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII"""        
        
        sam_file = self.create_temp_file(sam_content, '.sam')
        filtered_file = self.create_temp_file(filtered_content, '.fastq')
        output_file = os.path.join(self.temp_dir, 'output.sam')
        
        # Threshold 15.0: 
        # - Primary (chr1): scores ~18.0 and ~17.3 - both pass
        # - Secondary chr2: scores ~16.2 and ~15.1 - chr2 pair passes (max > 15)
        # - Secondary chr3: scores ~9.0 and ~7.2 - chr3 pair fails (max < 15)
        filter_viral_sam_memory_efficient(sam_file, filtered_file, output_file, 15.0)
        
        with open(output_file, 'r') as f:
            output_lines = f.readlines()
        
        # Should have primary pair (2) + chr2 secondary pair (2) = 4 lines
        assert len(output_lines) == 4
        
        # Check that we have the right references
        references = [SamAlignment.from_sam_line(line.strip()).rname for line in output_lines]
        assert 'chr1' in references  # Primary
        assert 'chr2' in references  # Secondary that passed
        assert 'chr3' not in references  # Secondary that failed
        
        # Check flag distribution: 2 primary (< 256), 2 secondary (>= 256)
        flags = [SamAlignment.from_sam_line(line.strip()).flag for line in output_lines]
        primary_flags = [f for f in flags if f < 256]
        secondary_flags = [f for f in flags if f >= 256]
        assert len(primary_flags) == 2
        assert len(secondary_flags) == 2

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
