// Labeled version - concatenate FASTA/FASTQ files preserving extension
process CONCATENATE_FASTN_LABELED {
    label "single"
    label "coreutils"
    input:
        tuple val(label), path(files)
        val(name)
    output:
        tuple val(label), path("${label}_${name}.*"), emit: output
        tuple val(label), path("${label}_input_${files[0]}"), emit: input
    script:
        // Extract extension from first file, handling double extensions like .fasta.gz
        def first_file = files[0].toString()
        if (!first_file.contains('.')) {
            throw new Exception("Input file ${first_file} has no extension")
        }
        def extension = first_file.substring(first_file.indexOf('.') + 1)
        """
        cat ${files} > ${label}_${name}.${extension}
        ln -s ${files[0]} ${label}_input_${files[0]} # Link input to output for testing
        """
}