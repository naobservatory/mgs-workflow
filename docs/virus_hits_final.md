# Columns of `virus_hits_final.tsv.gz`:

## Present in both single and paired end:
- `seq_id`: Name of read that was identified as a viral hit
- `aligner_genome_id`: GenBank ID for best viral genome match to read, as identified by our aligner (bowtie2 or minimap2). "Best" means having the highest length normalized alignment score. For paired-end data, this is the highest-scoring match across forward and reverse reads. 
- `aligner_taxid`: NCBI taxon ID for taxon best matching read, as identified by our aligner (bowtie2 or minimap2). For paired-end data, this is the highest-scoring match across forward and reverse reads. 
- `aligner_best_alignment_score`: Alignment score (directly from aligner) of best-scoring alignment (for paired-end data, score for forward read's best alignment). 
- `aligner_edit_distance`: Edit distance between read and aligned genome (for paired-end data, edit distance for forward read)
- `aligner_ref_start`: Location of start of alignment on reference (for paired-end data, location of forward read alignment on reference)
- `aligner_map_qual`: Mapping quality (MAPQ) as returned by bowtie2/minimap2 (for paired-end data, mapping quality of forward read)
- `aligner_cigar`: CIGAR string representing alignment of read to aligned genome (for paired-end data, CIGAR string for forward read). Note that this is the CIGAR string as returned by the aligner. If `query_rc_by_aligner` is true, you should reverse-complement the query sequence before comparing it against the CIGAR string. 
- `query_len`: Length of read, after trimming (for paired-end data, length of forward read)
- `query_seq`: Sequence of read (for paired-end data, sequence of forward read). Not reverse-complemented (we undo any reverse-complement performed by aligner).
- `query_rc_by_aligner`: Was query reverse-complemented by aligner? True/False  
- `query_qual`: PHRED quality scores for read (for paired-end data, quality scores for forward read). Like `query_seq`, not reverse-complemented. 
- `aligner_length_normalized_score`: Length-normalized alignment score: `aligner_best_alignment_score` divided by natural log of `query_len`. (For paired-end data, max of scores for forward and reverse reads.)
- `sample`: Sample name corresponding to read (from our samplesheet)

## Extra columns, for paired end only
- `aligner_genome_id_all`: GenBank IDs for viral genomes matching forward and reverse reads, joined by "/" 
- `aligner_taxid_all` NCBI taxon ID for taxons matching forward and reverse reads, joined by "/" 
- `aligner_fragment_length`: Inferred fragment length, as calculated by the aligner (NA if forward and reverse reads align to different genome IDs/align discordantly)
- `aligner_genome_id_fwd`: GenBank ID for viral genome matching forward read
- `aligner_genome_id_rev`: GenBank ID for viral genome matching reverse read
- `aligner_taxid_fwd`: NCBI taxon ID for taxon matching forward read 
- `aligner_taxid_rev`: NCBI taxon ID for taxon matching reverse read
- `aligner_best_alignment_score_rev`: Alignment score of best-scoring alignment of reverse read 
- `aligner_next_alignment_score`: Score of second-best alignment of forward read (currently absent for single-end data as we do not have minimap2 return multiple alignments)
- `aligner_next_alignment_score_rev`: Score of second-best alignment of reverse read
- `aligner_edit_distance_rev`: Edit distance between reverse read and aligned genome
- `aligner_ref_start_rev`: Location of start of alignment of reverse read on reference 
- `aligner_map_qual_rev`: Mapping quality (MAPQ) of reverse read
- `aligner_cigar_rev`: CIGAR string representing alignment of reverse read to aligned genome 
- `query_len_rev`: Length of reverse read
- `query_seq_rev`: Sequence of reverse read. Like query_seq, not reverse-complemented.
- `query_rc_by_aligner_rev`: Was (reverse-read) query reverse-complemented by aligner? True/False  
- `query_qual_rev`: PHRED quality scores for reverse read 
- `aligner_length_normalized_score_fwd`: Length-normalized alignment score for forward read
- `aligner_length_normalized_score_rev`: Length-normalized alignment score for reverse read
- `aligner_pair_status`: Pair status (UU for unpaired, CP for concordant pair, DP for discordant pair, UP for unaligned pair) 
- `kraken2_classified`: was the read classified by kraken2 (True/False)? (all kraken2 columns are present in paired-end only as only our paired-short-read pipeline runs kraken2)
- `kraken2_assigned_name`: Name of taxon to which read was assigned (e.g. "Human astrovirus"); "unclassified" if unassigned
- `kraken2_assigned_taxid`: NCBI taxon ID for taxon to which the read was assigned; 0 if unassigned
- `kraken2_assigned_host_virus`: Infection status annotation 0-3 (see https://github.com/naobservatory/mgs-workflow/blob/dev/docs/annotation.md) 
- `kraken2_length`: Length of merged read passed to kraken2
- `kraken2_encoded_hits`: space-delimited list indicating LCA mapping of each minimizer in the sequence (from kraken raw output)
- `bbmerge_frag_length`: Inferred fragment length, as calculated by bbmerge (NA if bbmerge could not calculate fragment length due to merge failure)

