// Validate that all samples in the grouping file are present in the viral hits table
// Filter out the samples from the grouping file that have 0 viral hits
process VALIDATE_GROUPING {
    label "python"
    label "single"
    input:
        tuple val(label), path(input_file), path(groups_file)
        val(join_field)
        val(output_label)
    output:
        tuple val(label), path("${label}_${output_label}_validated_${join_field}.tsv.gz"), emit: output
        tuple val(label), path("${label}_${output_label}_samples_with_zero_viral_hits.tsv"), emit: zero_vv_log
        tuple val(label), path("input_${input_file}"), path("input_${groups_file}"), emit: input
    shell:
        '''
        out_validated=!{label}_!{output_label}_validated_!{join_field}.tsv.gz
        out_zero_vv=!{label}_!{output_label}_zero_vv_samples.tsv
        validate_grouping.py !{input_file} !{groups_file} !{join_field} ${out_validated} ${out_zero_vv}
        # Link input files for testing
        ln -s !{input_file} input_!{input_file}
        ln -s !{groups_file} input_!{groups_file}
        '''
}
