/*
Given a TSV with two taxid columns, compute the vertical taxonomic distance
(number of parent-child steps) between the two taxids for each row. Rows
for which the two taxids are the same are given a distance of 0; rows for which
one taxid is not an ancestor of the other are given a distance of NA.
*/

process COMPUTE_TAXID_DISTANCE {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv) // Sorted TSV with two taxid columns
        val(taxid_field_1) // Column header for first taxid field
        val(taxid_field_2) // Column header for second taxid field
        val(distance_field) // Column header for new distance field
        path(nodes_db) // TSV containing taxonomic structure (mapping taxids to parent taxids)
    output:
        tuple val(sample), path("distance_${tsv}"), emit: output // Distance-computed TSV
        tuple val(sample), path("input_${tsv}"), emit: input // Input file for testing
    shell:
        '''
        # Set up and run Python script
        io="-i !{tsv} -o distance_!{tsv} -n !{nodes_db}"
        par="-t1 !{taxid_field_1} -t2 !{taxid_field_2} -d !{distance_field}"
        compute_taxid_distance.py ${io} ${par}
        # Link input file to output for testing
        ln -s !{tsv} input_!{tsv}
        '''
}