/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED_MERGED } from "../../../modules/local/subsetReads"
include { BLAST_LOCAL as BLAST_LOCAL_FWD } from "../../../modules/local/blast"
include { BLAST_LOCAL as BLAST_LOCAL_REV } from "../../../modules/local/blast"
include { FILTER_BLAST as FILTER_BLAST_FWD } from "../../../modules/local/filterBlast"
include { FILTER_BLAST as FILTER_BLAST_REV } from "../../../modules/local/filterBlast"
include { PAIR_BLAST } from "../../../modules/local/pairBlast"

/***********
| WORKFLOW |
***********/

workflow BLAST_VIRAL {
    take:
        viral_fasta
        blast_db_dir
        blast_db_prefix
        read_fraction
        blast_cpus
        blast_mem
        blast_filter_mem
    main:
        // Subset viral reads for BLAST
        if ( read_fraction < 1 ) {
            subset_ch = SUBSET_READS_PAIRED_MERGED(viral_fasta, read_fraction, "fasta")
        } else {
            subset_ch = viral_fasta
        }
        // BLAST forward and reverse reads separately against prepared DB
        blast_ch_1 = BLAST_LOCAL_FWD(subset_ch.map{it[0]}, blast_db_dir, blast_db_prefix, blast_cpus, blast_mem)
        blast_ch_2 = BLAST_LOCAL_REV(subset_ch.map{it[1]}, blast_db_dir, blast_db_prefix, blast_cpus, blast_mem)
        // Filter BLAST output (keeping forward and reverse separate)
        filter_ch_1 = FILTER_BLAST_FWD(blast_ch_1, blast_filter_mem, 5)
        filter_ch_2 = FILTER_BLAST_REV(blast_ch_2, blast_filter_mem, 5)
        // Combine alignment information across pairs
        filter_out_1 = filter_ch_1.map{ file ->
            def newFile = file.parent / "fwd_${file.name}"
            file.copyTo(newFile)
            return newFile
            }
        filter_out_2 = filter_ch_2.map{ file ->
            def newFile = file.parent / "rev_${file.name}"
            file.copyTo(newFile)
            return newFile
            }
        pair_ch = PAIR_BLAST(filter_out_1, filter_out_2)
    emit:
        blast_subset = subset_ch
        blast_paired = pair_ch
}
