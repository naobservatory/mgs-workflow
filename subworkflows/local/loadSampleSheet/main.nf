/***********
| WORKFLOW |
***********/

workflow LOAD_SAMPLESHEET {
    take:
        sample_sheet
        grouping
        single_end
    main:
        // Start time
        start_time = new Date()
        start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

        // Define expected headers
        def expected_headers_se = ['sample', 'fastq']
        def expected_headers_pe = ['sample', 'fastq_1', 'fastq_2']
        
        // Read actual headers
        def headers = file(sample_sheet).readLines().first().tokenize(',')*.trim()
        
        // Validate headers based on single_end and grouping parameters
        def required_headers = single_end ? expected_headers_se : expected_headers_pe
        if (grouping) required_headers = required_headers + ['group']
        
        // Check if headers match exactly
        if (headers != required_headers) {
            throw new Exception("""Invalid samplesheet header. 
                Expected: ${required_headers.join(', ')}
                Found: ${headers.join(', ')}
                Please ensure the samplesheet has the correct columns in the specified order.""".stripIndent())
        }

        // Check if grouping column exists in samplesheet
        check_grouping = file(sample_sheet).readLines()[0].contains('group')
        if (grouping != check_grouping) {
            if (grouping && !check_grouping) {
                throw new Exception("Grouping enabled in config file, but group column absent from samplesheet.")
            } else if (!grouping && check_grouping) {
                throw new Exception("Grouping is not enabled in config file, but group column is present in the samplesheet.")
            }
        }


        if (single_end) {
            if (grouping) {
                samplesheet = Channel
                    .fromPath(sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq), row.group) }
                samplesheet_ch = samplesheet.map { sample, read, group -> tuple(sample, [read]) }
                group_ch = samplesheet.map { sample, read, group -> tuple(sample, group) }
            } else {
                samplesheet = Channel
                    .fromPath(sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq)) }
                samplesheet_ch = samplesheet.map { sample, read -> tuple(sample, [read]) }
                group_ch = Channel.empty()
            }
        } else {
            if (grouping) {
                samplesheet = Channel
                    .fromPath(sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2), row.group) }
                samplesheet_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, [read1, read2]) }
                group_ch = samplesheet.map { sample, read1, read2, group -> tuple(sample, group) }
            } else {
                samplesheet = Channel
                    .fromPath(sample_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
                samplesheet_ch = samplesheet.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
                group_ch = Channel.empty()
                }
            }
    emit:
        samplesheet = samplesheet_ch
        group = group_ch
        start_time_str = start_time_str
}
