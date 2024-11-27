// Copy a file to a new location with a custom path
process CONCAT_GROUP {
    label "coreutils"
    label "single"
    input:
        tuple val(samples), path(fastq_1_list), path(fastq_2_list), val(group)

    output:
        tuple val(group), path("${group}_R{1,2}.fastq.gz")

    script:
        """
        cat ${fastq_1_list.join(' ')} > ${group}_R1.fastq.gz
        cat ${fastq_2_list.join(' ')} > ${group}_R2.fastq.gz
        """
}
