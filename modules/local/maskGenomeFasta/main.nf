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
	path("${name_pattern}-mask-adapters-entropy.stats.txt), emit: log1
	path("${name_pattern}-mask-polyx.stats.txt), emit: log2
    shell:
        '''
        zcat !{collated_genomes} | grep "^>" | grep -vif !{patterns_exclude} | sed 's/>//' > names.txt
        seqtk subseq !{collated_genomes} names.txt | gzip -c > !{name_pattern}.fasta.gz
	# Define input/output
	in=!{filtered_genomes}
	out1=intermediate-masking.fasta.gz
	out2=!{name_pattern}-masked.fasta.gz
	ref=!{adapters}
	stats1=!{name_pattern}-mask-adapters-entropy.stats.txt
	stats2=!{name_pattern}-mask-polyx.txt
	par1="k=!{k} hdist=!{hdist} mm=f mask=N rcomp=t entropy=!{entropy} entropymask=t mink=8 hdist2=1"
	par2="k=!{polyx_len} hdist=0 mm=f mask=N rcomp=F"
 	# Simplest way to mask polyX regions is just to pass them as literals, 
	# e.g. "AAAAA,CCCCC,GGGGG,TTTTT" for polyx_len=5
	polyx='A' * polyx_len + ',' + 'C' * polyx_len + ',' + 'G' * polyx_len + ',' + 'T' * polyx_len 
	# Execute masking in sequence: first adapter/entropy masking, then polyX masking
	bbduk.sh in=${in} out=${out1} ref=${ref} stats=${stats1} ${par1}
	bbduk.sh in=${out1} out=${out2} literal=${polyx} stats=${stats2} ${par2}         
	'''
}
