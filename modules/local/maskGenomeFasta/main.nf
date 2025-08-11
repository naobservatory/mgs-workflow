// Filter genomes to exclude specific patterns in sequence headers
process MASK_GENOME_FASTA {
    label "large"
    label "BBTools"
    input:
        path(filtered_genomes)
        path(adapters)
        val(params_map)
    output:
        path("${name_pattern}-masked.fasta.gz"), emit: masked
		path("${name_pattern}-mask-adapters-entropy.stats.txt"), emit: log1
		path("${name_pattern}-mask-polyx.stats.txt"), emit: log2
	script:
	// Extract parameters from map
	k = params_map.k
	hdist = params_map.hdist
	entropy = params_map.entropy
	polyx_len = params_map.polyx_len
	name_pattern = params_map.name_pattern
	// Groovy runs first – build the poly-X literal once
	def polyx = ['A','C','G','T'].collect { it * (polyx_len as int) }.join(',')
	"""
	# first pass – adapters + entropy
	bbduk.sh \
		in=${filtered_genomes} \
		out=intermediate-masking.fasta.gz \
		ref=${adapters} \
		stats=${name_pattern}-mask-adapters-entropy.stats.txt \
		k=${k} hdist=${hdist} mm=f mask=N rcomp=t \
		entropy=${entropy} entropymask=t mink=8 hdist2=1

	# second pass – poly-X masking
	bbduk.sh \
		in=intermediate-masking.fasta.gz \
		out=${name_pattern}-masked.fasta.gz \
		literal=${polyx} \
		stats=${name_pattern}-mask-polyx.stats.txt \
		k=${polyx_len} hdist=0 mm=f mask=N rcomp=F
	"""
}