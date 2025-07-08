## Aligner LCA Output Columns (excluding `seq_id`)

The aligner's `LCA_TSV` output is structured around three distinct groups of sequence alignments, identified by seq_id: All, Natural, and Artificial. Each of these groups includes the same set of eight statistical columns that describe the alignments within that category.

### Data Categories & Definitions

First, let's define the key terms used to categorize the alignment data:

- Natural Alignments: These are alignments to taxa that do not descend from the NCBI taxid 81077. This group represents naturally occurring sequences.

- Artificial Alignments: These are alignments to taxa that do descend from the NCBI taxid 81077, which is designated as the parent for artificial or synthetic sequences.

- All Alignments: This category is a superset, encompassing both Natural and Artificial alignments for a given sequence.

- Classified Taxa: A key distinction within each group is whether a taxon is "classified." A taxon is considered classified if its NCBI taxon name does not contain the terms "sp." or "unclassified." This helps filter out ambiguous or low-confidence assignments.

### Data Processing Logic

The calculation of the Lowest Common Ancestor (LCA) and the selection of the top-scoring taxon follow specific rules:

- Top Taxid Selection: When multiple alignments have the same top score, a classified taxon will always be chosen over an unclassified one. If there's a tie between multiple classified taxa, the one that appears first is selected.
- LCA Calculation:
  - If the top-scoring alignment is a classified taxon, the LCA is determined using only the classified taxa from the alignment set.
  - If the top-scoring alignment is an unclassified taxon, the LCA calculation will include both classified taxa and any unclassified taxa that share the top score.

### Output Column Descriptions

There are eight core statistical columns, which are repeated for each of the three data categories (`All`, `Natural`, and `Artificial`). The suffix of each column name (e.g., `_natural`, `_artificial`) indicates which category it describes.

#### Core Column Definitions:

* `_taxid_lca`: The NCBI taxon ID for the Lowest Common Ancestor (LCA) of all alignments in the group.
* `_n_assignments`: The total number of alignments in the group for that read.
* `_n_assignments_classified`: The subset of `_n_assignments` that are to classified taxa.
* `_taxid_top`: The taxon ID of the top-scoring alignment in the group.
* `_taxid_top_classified`: A boolean (`True`/`False`) value indicating if the top-scoring taxon is classified.
* `_length_normalized_score_min`: The minimum alignment score in the group.
* `_length_normalized_score_max`: The maximum alignment score in the group.
* `_length_normalized_score_mean`: The average alignment score for the group.

#### Column Groups:

* **All** (e.g., `aligner_taxid_lca`, `aligner_n_assignments`): These columns provide statistics for *all* alignments (Natural + Artificial).
* **Natural** (e.g., `aligner_taxid_lca_natural`, `aligner_n_assignments_natural`): These columns provide the same statistics but *only* for the Natural alignments.
* **Artificial** (e.g., `aligner_taxid_lca_artificial`, `aligner_n_assignments_artificial`): These columns provide the same statistics but *only* for the Artificial alignments.


## Processed Bowtie2/Minimap2 TSV Output Columns (excluding `seq_id`)

- `genome_id`: GenBank ID for best viral genome match to read, as identified by our aligner (bowtie2 or minimap2). "Best" means having the highest length normalized alignment score. For paired-end data, this is the highest-scoring match across forward and reverse reads. 
- `genome_id_all`: GenBank IDs for viral genomes matching forward and reverse reads, joined by "/" 
- `taxid`: NCBI taxon ID for taxon best matching read, as identified by our aligner (bowtie2 or minimap2). For paired-end data, this is the highest-scoring match across forward and reverse reads. 
- `taxid_all`: NCBI taxon ID for taxons matching forward and reverse reads, joined by "/" 
- `fragment_length`: Inferred fragment length, as calculated by the aligner (NA if forward and reverse reads align to different genome IDs/align discordantly)
- `best_alignment_score`: Alignment score (directly from aligner) of best-scoring alignment (for paired-end data, score for forward read's best alignment). 
- `next_alignment_score`: Alignment score of second-best alignment (for paired-end data, score for forward read's second-best alignment; NA for minimap2)
- `genome_id_fwd`: GenBank ID for viral genome matching forward read
- `genome_id_rev`: GenBank ID for viral genome matching reverse read
- `taxid_fwd`: NCBI taxon ID for taxon matching forward read 
- `taxid_rev`: NCBI taxon ID for taxon matching reverse read
- `best_alignment_score_rev`: Alignment score of best-scoring alignment of reverse read 
- `next_alignment_score_rev`: Score of second-best alignment of reverse read
- `edit_distance`: Edit distance between read and aligned genome (for paired-end data, edit distance for forward read)
- `edit_distance_rev`: Edit distance between reverse read and aligned genome
- `ref_start`: Location of start of alignment on reference (for paired-end data, location of forward read alignment on reference)
- `ref_start_rev`: Location of start of alignment of reverse read on reference 
- `map_qual`: Mapping quality (MAPQ) as returned by bowtie2/minimap2 (for paired-end data, mapping quality of forward read)
- `map_qual_rev`: Mapping quality (MAPQ) of reverse read
- `cigar`: CIGAR string representing alignment of read to aligned genome (for paired-end data, CIGAR string for forward read). Note that this is the CIGAR string as returned by the aligner. If `query_rc` is true, you should reverse-complement the query sequence before comparing it against the CIGAR string. 
- `cigar_rev`: CIGAR string representing alignment of reverse read to aligned genome 
- `query_len`: Length of read, after trimming (for paired-end data, length of forward read)
- `query_len_rev`: Length of reverse read
- `query_seq`: Sequence of read (for paired-end data, sequence of forward read). Not reverse-complemented (we undo any reverse-complement performed by aligner).
- `query_seq_rev`: Sequence of reverse read. Like query_seq, not reverse-complemented.
- `query_rc`: Was query reverse-complemented by aligner? True/False  
- `query_rc_by_rev`: Was (reverse-read) query reverse-complemented by aligner? True/False  
- `query_qual`: PHRED quality scores for read (for paired-end data, quality scores for forward read). Like `query_seq`, not reverse-complemented. 
- `query_qual_rev`: PHRED quality scores for reverse read 
- `length_normalized_score_fwd`: Length-normalized alignment score for forward read
- `length_normalized_score_rev`: Length-normalized alignment score for reverse read
- `length_normalized_score`: Length-normalized alignment score: `best_alignment_score` divided by natural log of `query_len`. (For paired-end data, max of scores for forward and reverse reads.)
- `pair_status`: Pair status (UU for unpaired, CP for concordant pair, DP for discordant pair, UP for unaligned pair) 
- `classification`: Alignment type (primary, secondary or supplementary)

