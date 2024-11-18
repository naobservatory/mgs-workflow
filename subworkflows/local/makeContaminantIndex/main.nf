/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

//  include { DOWNLOAD_GENOME as DOWNLOAD_COW } from "../../../modules/local/downloadGenome"
//  include { DOWNLOAD_GENOME as DOWNLOAD_PIG } from "../../../modules/local/downloadGenome"
//  include { DOWNLOAD_GENOME as DOWNLOAD_MOUSE } from "../../../modules/local/downloadGenome"
//  include { DOWNLOAD_GENOME as DOWNLOAD_CARP } from "../../../modules/local/downloadGenome"
//  include { DOWNLOAD_GENOME as DOWNLOAD_ECOLI } from "../../../modules/local/downloadGenome"
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
            .of(genome_urls)
            .splitCsv(sep: ' ')  // Split on spaces
            .flatten()  // Flatten the resulting list
            .map { url -> 
                def uuid = UUID.randomUUID().toString().take(8)  // Takes first 8 chars of UUID
                tuple(url, "genome_${uuid}")
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
