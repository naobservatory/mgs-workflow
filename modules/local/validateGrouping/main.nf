// Filter out the samples from the grouping file that have 0 viral hits
process VALIDATE_GROUPING {
    label "python"
    label "single"
    input:
        tuple val(label), path(input_file), path(groups_file)
    output:
        tuple val(label), path(input_file), path("validated_${groups_file}"), emit: output
        tuple val(label), path("zero_vv_${label}_${input_file}"), emit: zero_vv_log
        tuple val(label), path("input_${input_file}"), path("input_${groups_file}"), emit: input
    shell:
        '''
        out_validated=validated_!{groups_file}
        out_zero_vv=zero_vv_!{label}_!{input_file}
        validate_grouping.py !{input_file} !{groups_file} ${out_validated} ${out_zero_vv}
        # Link input files for testing
        ln -s !{input_file} input_!{input_file}
        ln -s !{groups_file} input_!{groups_file}
        '''
}
