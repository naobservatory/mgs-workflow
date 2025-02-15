process MARK_DUPLICATE_FASTA_HEADERS {
    label "biopython"
    label "single"

    input:
        tuple path(reads)

    output:
        path("uniq_headers.fasta.gz"), emit: reads
    shell:
        '''
        #!/usr/bin/env python
        from Bio import SeqIO
        import gzip
        from collections import defaultdict

        # Track header occurrences
        header_counts = defaultdict(int)
        records = []

        # First pass - count headers
        with gzip.open("!{reads}", 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header_counts[record.id] += 1
                records.append(record)

        # Second pass - modify duplicate headers
        with gzip.open("uniq_headers.fasta.gz", 'wt') as outf:
            for record in records:
                if header_counts[record.id] > 1:
                    header_counts[record.id] -= 1
                    suffix = chr(96 + header_counts[record.id])  # a, b, c, etc.
                    record.id = f"{record.id}-{suffix}"
                    record.description = record.id
                SeqIO.write(record, outf, "fasta")
        '''
}