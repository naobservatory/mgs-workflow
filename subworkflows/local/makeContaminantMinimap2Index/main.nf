/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DOWNLOAD_GENOME } from "../../../modules/local/downloadGenome"
include { CONCATENATE_FASTA_GZIPPED } from "../../../modules/local/concatenateFasta"
include { MINIMAP2_INDEX } from "../../../modules/local/minimap2"

/***********
| WORKFLOW |
***********/

workflow MAKE_CONTAMINANT_MINIMAP2_INDEX {
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
        minimap2_ch = MINIMAP2_INDEX(genome_ch, "mm2-other-index")
    emit:
        mm2 = minimap2_ch
}
