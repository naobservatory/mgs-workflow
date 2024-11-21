/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME } from "../../../modules/local/downloadGenome"
include { CONCATENATE_FASTA_GZIPPED } from "../../../modules/local/concatenateFasta"
include { BBMAP_INDEX } from "../../../modules/local/bbmap"
include { BOWTIE2_INDEX } from "../../../modules/local/bowtie2"

/***********
| WORKFLOW |
***********/

workflow MAKE_CONTAMINANT_INDEX {
    take:
        genome_urls
        contaminants_path
    main:
        // Download reference genomes
        ref_ch = channel
            .fromList(genome_urls.entrySet())
            .map { entry -> 
                tuple(entry.value, entry.key)  // (url, name)
            }

        downloaded_ch = DOWNLOAD_GENOME(ref_ch)

        combined_ch = downloaded_ch
            .mix(channel.fromPath(contaminants_path))
            .collect()

        // Then use combined_ch for the rest of your workflow
        genome_ch = CONCATENATE_FASTA_GZIPPED(combined_ch, "ref_concat")

        // Make indexes
        bbmap_ch = BBMAP_INDEX(genome_ch, "bbm-other-index")
        bowtie2_ch = BOWTIE2_INDEX(genome_ch, "bt2-other-index")
    emit:
        bbm = bbmap_ch
        bt2 = bowtie2_ch
}
