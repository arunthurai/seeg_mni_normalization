#from Mau's previous work see https://github.com/mcespedes99/MNI_Registration

rule synthstrip_T1w:
    input:
        t1=inputs.path["base_T1w"],
    output:
        mask=temp(
            bids(
                root='work',
                **inputs.wildcards['base_T1w'],
                desc="temp",
                suffix="mask.nii.gz"
            )
        ),
    group:
        "anat"
    threads: 8
    shell:
        "mri_synthstrip -i {input.t1} -m {output.mask} --no-csf"


rule fixheader_synthstrip:
    input:
        t1=inputs.path["base_T1w"],
        mask=rules.synthstrip_T1w.output.mask,
    output:
        mask=bids(
            root='work',
            **inputs.wildcards['base_T1w'],
            desc="brain",
            suffix="mask.nii.gz"
        ),
    group:
        "anat"
    shell:
        "c3d {input.t1} {input.mask} -copy-transform -o {output.mask}"


rule n4_bias_correction:
    input:
        ref_t1 = inputs.path['base_T1w'],
        mask=rules.fixheader_synthstrip.output.mask,

    output:
        n4 = bids(
                root = 'work',
                suffix = 'T1w.nii.gz',
                desc = 'preproc',
                **inputs.wildcards['base_T1w'],
            ),
    group:
        "anat"
    threads: 8
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input.ref_t1} -x {input.mask} -o {output.n4}"

rule mask_subject_T1w:
    input:
        T1w=rules.n4_bias_correction.output.n4,
        mask=rules.fixheader_synthstrip.output.mask,
    output:
        final_T1w=bids(
            root='work',
            **inputs.wildcards['base_T1w'],
            suffix="T1w.nii.gz",
            desc="masked"
        ),
    group:
        "anat"
    shell:
        "c3d {input} -multiply -o {output}"