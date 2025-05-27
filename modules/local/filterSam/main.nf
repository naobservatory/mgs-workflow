process FILTER_SAM {
    label "pysam"
    input:
        tuple val(sample), path(sam)
        val(threshold)
    output:
        tuple val(sample), path("${sample}_filtered.sam"), emit: sam
    script:
        """
        #!/usr/bin/env python3

        import pysam
        import math
        import sys

        input_file = "${sam}"
        output_file = "${sample}_filtered.sam"
        threshold = ${threshold}

        def extract_alignment_score(read):
            try:
                return read.get_tag("AS")
            except KeyError:
                for tag_name, tag_value in read.get_tags():
                    if tag_name == "AS":
                        return int(tag_value)
            return None


        with pysam.AlignmentFile(input_file, "r") as infile:
            with pysam.AlignmentFile(output_file, "w", template=infile) as outfile:
                for read in infile:
                    query_len = len(str(read.query_sequence))
                    alignment_score = extract_alignment_score(read)
                    if alignment_score is None:
                        continue
                    if alignment_score is not None and (alignment_score / math.log(query_len)) >= threshold:
                        outfile.write(read)
        """
}
