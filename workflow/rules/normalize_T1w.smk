#from Mau's previous work see https://github.com/mcespedes99/MNI_Registration


def get_template_prefix(root, subj_wildcards, template):
    """creates prefix for template files, including subject/session wildcards
    so that DAGs for each subject/session are kept independent.
        e.g.: sub-001/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym"""

    path_entities = bids(root=root, **subj_wildcards).split("/")[
        :-1
    ]  # leave out the file prefix

    path_entities.append(f"tpl-{template}")  # sub-folder name
    path_entities.append(f"tpl-{template}")  # file prefix

    return "/".join(path_entities)


rule mask_template_t1w:
    input:
        t1=str(Path(workflow.basedir).parent /config["template_t1w"]),
        mask=str(Path(workflow.basedir).parent /config["template_mask"]),
    output:
        t1=temp(
        get_template_prefix(
            root='work', subj_wildcards=inputs.wildcards['base_T1w'], template=config["template"]
        )
        )
        + "_desc-masked_T1w.nii.gz",
    group:
        "anat"
    shell:
        "c3d {input} -multiply -o {output}"

if config['greedy']:
    rule greedy_t1_to_template:
        input:
            flo=[
                bids(
                    root='work',
                    **inputs.wildcards['base_T1w'],
                    suffix="T1w.nii.gz",
                    desc="masked",
                )
            ],
                
            ref=[
                get_template_prefix(
                    root='work', subj_wildcards=inputs.wildcards['base_T1w'], template=config["template"]
                )
                + "_desc-masked_T1w.nii.gz",
            ],
            #seg=[
                #get_template_prefix(
                #root='work', subj_wildcards=inputs.wildcards['base_T1w'], template=config["template"]
                #),
            #]
        params:
            input_fixed_moving=lambda wildcards, input: [
                f"-i {fixed} {moving}" for fixed, moving in zip(input.ref, input.flo)
            ],
            input_moving_warped=lambda wildcards, input, output: [
                f"-rm {moving} {warped}"
                for moving, warped in zip(input.flo, output.warped_flo)
            ],
            affine_iterations="100x50x10",
            fluid_iterations="100x50x10",  #default 100x50x10
            gradient_sigma="1.732vox",  #default 1.732vox
            warp_sigma="0.707vox",  #default 0.707vox
            timestep="1.0",  #default 1.0
        output:
            warp=bids(
                root=config['output_dir'],
                datatype="warps",
                suffix="warp.nii.gz",
                from_="subject",
                to=config["template"],
                **inputs.wildcards['base_T1w']
            ),
            invwarp=bids(
                root=config['output_dir'],
                datatype="warps",
                suffix="invwarp.nii.gz",
                from_="subject",
                to=config["template"],
                **inputs.wildcards['base_T1w']
            ),
            warped_flo=[
                bids(
                    root=config['output_dir'],
                    datatype="anat",
                    suffix="T1w.nii.gz",
                    space=config["template"],
                    desc="greedy",
                    **inputs.wildcards['base_T1w']
                )
            ],
            affine=bids(
                root=config['output_dir'],
                datatype="warps",
                suffix="affine.txt",
                from_="subject",
                to=config["template"],
                desc="itk",
                **inputs.wildcards['base_T1w']
            ),
            affine_xfm_ras=bids(
                root=config['output_dir'],
                datatype="warps",
                suffix="affine.txt",
                from_="subject",
                to=config["template"],
                desc="ras",
                **inputs.wildcards['base_T1w']
            ),
        threads: 8
        resources:
            mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
            time=60,  # 1 hrs
        group:
            "anat"
        log:
            bids(
                root="logs",
                suffix="greedy.log",
                template=config["template"],
                **inputs.wildcards['base_T1w']
            ),
        shell:
            #affine first
            "greedy -d 3 -threads {threads} -a -m NCC 2x2x2 {params.input_fixed_moving} -o {output.affine_xfm_ras} -ia-image-centers -n {params.affine_iterations} &> {log} && "
    
            "greedy -d 3 -threads {threads} -m NCC 2x2x2 {params.input_fixed_moving} -it {output.affine_xfm_ras} -o {output.warp} -oinv {output.invwarp} -n {params.fluid_iterations} -s {params.gradient_sigma} {params.warp_sigma} -e {params.timestep} &>> {log} && "
    
            "c3d_affine_tool {output.affine_xfm_ras} -oitk {output.affine} &>> {log} && "
    
            "greedy -d 3 -threads {threads} -rf {input.ref[0]} {params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} &>> {log}"
            #then deformable:
            #then convert affine to itk format that ants uses
            #and finally warp the moving image

if config['ants']:
    rule ants_T1_template:
        input:
            flo= bids(
                root='work',
                **inputs.wildcards['base_T1w'],
                suffix="T1w.nii.gz",
                desc="masked"
                ),
            ref = get_template_prefix(
                root='work', 
                subj_wildcards=inputs.wildcards['base_T1w'], 
                template=config["template"]
                ) + "_desc-masked_T1w.nii.gz"
        #params:
           
        output:
            prefix = bids(
                root = config['output_dir'],
                from_='sub',
                to = config['template'],
                **inputs.wildcards['base_T1w']
            ),
            warp=bids(
                root=config['output_dir'],
                datatype="warps_ants",
                suffix="warp.nii.gz",
                from_="subject",
                to=config["template"],
                **inputs.wildcards['base_T1w']
            ),
            invwarp=bids(
                root=config['output_dir'],
                datatype="warps_ants",
                suffix="invwarp.nii.gz",
                from_="subject",
                to=config["template"],
                **inputs.wildcards['base_T1w']
            ),
            warped_flo=bids(
                    root=config['output_dir'],
                    datatype="anat",
                    suffix="T1w.nii.gz",
                    space=config["template"],
                    desc="ants",
                    **inputs.wildcards['base_T1w']
                ),
            affine=bids(
                root=config['output_dir'],
                datatype="warps_ants",
                suffix="affine.txt",
                from_="subject",
                to=config["template"],
                desc="itk",
                **inputs.wildcards['base_T1w']
            ),
            affine_xfm_ras=bids(
                root=config['output_dir'],
                datatype="warps_ants",
                suffix="affine.txt",
                from_="subject",
                to=config["template"],
                desc="ras",
                **inputs.wildcards['base_T1w']
            ),
        threads: 8
        resources:
            mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
            time=60,  # 1 hrs
        group:
            "anat"
        log:
            bids(
                root="logs",
                suffix="ants.log",
                template=config["template"],
                **inputs.wildcards['base_T1w']
            ),
        shell:
            'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} antsRegistration '
            '--write-composite-transform '
            '-d 3 --float 1 -u 1 -w [0.01,0.99] -z 1 '
            '-o [{output.prefix}, {output.warped_flo}] '
            '-t Rigid[0.1] '
            '-m CC[{input.ref}, {input.flo}, 1, 32, Regular, 0.25] '
            '-c [1000x500x250x0,1e-6,10] '
            '-f 6x4x2x1 '
            '-s 4x2x1x0 '
            '-t Affine[0.1] '
            '-m CC[{input.ref}, {input.flo}, 1, 32, Regular, 0.25] '
            '-c [1000x500x250x0,1e-6,10] '
            '-f 6x4x2x1 '
            '-s 4x2x1x0 '
            '-t SyN[0.1,3,0] '
            '-m CC[{input.ref}, {input.flo}, 1, 4] '
            '-c [100x100x70x50x10,1e-9,10] '
            '-f 12x6x4x2x1 '
            '-s 6x3x2x1x0vox'

#rule warp_dseg_from_template:

