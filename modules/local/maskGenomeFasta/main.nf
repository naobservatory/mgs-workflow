// Filter genomes to exclude specific patterns in sequence headers
process MASK_GENOME_FASTA {
    label "large"
    label "BBTools"
    input:
        path(filtered_genomes)
        path(adapters)
        val(k)
        val(hdist)
        val(entropy)
        val(polyx_len)
        val(name_pattern)
    output:
        path("${name_pattern}-masked.fasta.gz"), emit: masked
	path("${name_pattern}-mask-adapters-entropy.stats.txt"), emit: log1
	path("${name_pattern}-mask-polyx.stats.txt"), emit: log2
    shell:
 	// Simplest way to mask polyX regions is just to pass them as literals,
	// e.g. "AAAAA,CCCCC,GGGGG,TTTTT" for polyx_len=5
	polyx = ['A', 'C', 'G', 'T'].collect { it * (polyx_len as int) }.join(',')
        '''
	# Define input/output
	in=!{filtered_genomes}
	out1=intermediate-masking.fasta.gz
	out2=!{name_pattern}-masked.fasta.gz
	ref=!{adapters}
	stats1=!{name_pattern}-mask-adapters-entropy.stats.txt
	stats2=!{name_pattern}-mask-polyx.stats.txt
	par1="k=!{k} hdist=!{hdist} mm=f mask=N rcomp=t entropy=!{entropy} entropymask=t mink=8 hdist2=1"
	par2="k=!{polyx_len} hdist=0 mm=f mask=N rcomp=F"
	# Execute masking in sequence: first adapter/entropy masking, then polyX masking
	bbduk.sh in=${in} out=${out1} ref=${ref} stats=${stats1} ${par1}
	bbduk.sh in=${out1} out=${out2} literal=!{polyx} stats=${stats2} ${par2}
	'''
}
