# Explanation of our custom lowest common ancestor algorithm (LCA_TSV)

In deciding which alignments to include in our LCA calculation, we distinguish between taxids with a fully-resolved place in the taxonomy tree (which we call "classified") and those with ambiguity in their position (which we call "unclassified"). A taxon was designated as unclassified if its name, or that of any of its ancestors, contained the strings "unclassified" or "sp.". Unclassified taxa were only included in the LCA calculation if they had the highest score across all taxid assignments for a given species and were not tied for highest score with any classified taxa[^unclassified]. 

[^unclassified]: This adjustment to naive LCA was necessary because many sequences in Genbank are initially assigned to such "unclassified" taxa, even if further investigation would result in a more informative classification. As a result, a read from a given virus often aligns to a mixture of classified and unclassified taxa; in this situation, naive LCA results in a very high-level assignment (e.g. at the family level or even higher) that didn't accurately reflect the sum of our information about the read. Discarding alignments to unclassified taxa unless they are clearly the best alignment for a sequence avoids this issue.

In addition, we distinguish in our LCA calculation between artificial sequences (those descended from [taxid 81077](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=81077&lvl=3&lin=f&keep=1&srchmode=1&unlock), "artificial sequences") and natural sequences (all others). We calculate LCA assignments (and associated statistics) separately for natural sequences, artificial sequences, and all sequences together; by default, only the natural LCA is used for downstream analysis[^artificial].

[^artificial]: Artificial sequences are often generated from natural sequences, and can easily interfere with the LCA process. For example, an adenovirus read will match both adenovirus and a range of adenovirus-derived engineered vectors; taking the LCA of these two groups of sequences will result in an output assignment of root (taxid 1) even if there is no reason to believe the sequence actually came from a vector.

The algorithm operates as follows:

- Top taxid selection: Among all alignments within the same seq_id with the highest score:
  - Classified taxa are always preferred over unclassified ones
  - If multiple classified taxa tie for the top score, the first encountered is selected
  - If only unclassified taxa have the top score, the first encountered is selected
- LCA calculation rules:
  - If the top taxid is classified: Only classified taxa are included in the LCA calculation (all unclassified alignments are excluded)
  - If the top taxid is unclassified: The LCA includes all classified taxa PLUS any unclassified taxa that share the top score
- Output categories: The algorithm produces three separate LCA calculations:
  - Natural: Excludes all artificial alignments 
  - Artificial: Only artificial alignments
  - Combined: Includes both natural and artificial alignments

More information on the output of `LCA_TSV` for the `RUN` workflow can be found in `docs/lca_intermediates.md`.

