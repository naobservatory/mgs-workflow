/***********
| WORKFLOW |
***********/

workflow LOAD_SAMPLESHEET_STREAMED {
    take:
        sample_sheet
        single_end
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Validate headers
        def expected_headers_se = ['sample', 'fastq']
        def expected_headers_pe = ['sample', 'fastq_1', 'fastq_2']
        def required_headers = single_end ? expected_headers_se : expected_headers_pe
        def headers = file(sample_sheet).readLines().first().tokenize(',')*.trim()
        if (headers != required_headers) {
            throw new Exception("""Invalid samplesheet header. 
                Expected: ${required_headers.join(', ')}
                Found: ${headers.join(', ')}
                Please ensure the samplesheet has the correct columns in the specified order.""".stripIndent())
        }

        // Construct samplesheet channel
        if (single_end) {
            samplesheet = Channel
                .fromPath(sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq)) }
            samplesheet_ch = samplesheet.map { sample, read -> tuple(sample, [read]) }
        } else {
            samplesheet = Channel
                .fromPath(sample_sheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
            samplesheet_ch = samplesheet.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
        }

    emit:
        samplesheet = samplesheet_ch
        start_time_str = start_time_str
}
