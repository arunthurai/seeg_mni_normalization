
rule transform_contacts:
    input:
        contacts = inputs.path['seeg_contacts'],
        actual = inputs.path['seeg_actual'],
        planned = inputs.path['seeg_planned'],
        deformable_warp = bids(
            root=config['output_dir'],
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="subject",
            to=config["template"],
            **inputs.wildcards['base_T1w']
            ),
        affine = bids(
            root=config['output_dir'],
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="ras",
            **inputs.wildcards['seeg_contacts']
        ),

    output:
        deformable_contacts = bids(
         root = config['output_dir'],
         desc = 'warp',
         space=config["template"],
         suffix = "SEEGA.fcsv",
         **inputs.wildcards['seeg_contacts']
        ),
        deformable_actual = bids(
         root = config['output_dir'],
         desc = 'warp',
         space=config["template"],
         suffix = "actual.fcsv",
         **inputs.wildcards['seeg_contacts']
        ),
        deformable_planned = bids(
         root = config['output_dir'],
         desc = 'warp',
         space=config["template"],
         suffix = "planned.fcsv",
         **inputs.wildcards['seeg_contacts']
        ),

        affine_contacts = bids(
         root = config['output_dir'],
         desc = 'affine',
         space=config["template"],
         suffix = "SEEGA.fcsv",
         **inputs.wildcards['seeg_contacts']
        ),
        affine_actual = bids(
         root = config['output_dir'],
         desc = 'affine',
         space=config["template"],
         suffix = "actual.fcsv",
         **inputs.wildcards['seeg_contacts']
        ),
        affine_planned = bids(
         root = config['output_dir'],
         desc = 'affine',
         space=config["template"],
         suffix = "planned.fcsv",
         **inputs.wildcards['seeg_contacts']
        )
    script:
        "../scripts/warp_coords.py"