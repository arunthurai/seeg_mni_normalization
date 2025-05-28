#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 19:04:41 2025

@author: arun
"""
#see https://github.com/afids/afids-pred/blob/main/notebooks/2.x_LeadDBS_apply_transforms_stn.ipynb
import numpy as np
import pandas as pd
import scipy.io
import os
import SimpleITK as sitk
import re

def determineFCSVCoordSystem(input_fcsv,overwrite_fcsv=False):
	# need to determine if file is in RAS or LPS
	# loop through header to find coordinate system
	coordFlag = re.compile('# CoordinateSystem')
	verFlag = re.compile('# Markups fiducial file version')
	headFlag = re.compile('# columns')
	coord_sys=None
	headFin=None
	ver_fin=None

	with open(input_fcsv, 'r') as myfile:
		firstNlines=myfile.readlines()[0:3]

	for row in firstNlines:
		row=re.sub("[\s\,]+[\,]","",row).replace("\n","")
		cleaned_dict={row.split('=')[0].strip():row.split('=')[1].strip()}
		if None in list(cleaned_dict):
			cleaned_dict['# columns'] = cleaned_dict.pop(None)
		if any(coordFlag.match(x) for x in list(cleaned_dict)):
			coord_sys = list(cleaned_dict.values())[0]
		if any(verFlag.match(x) for x in list(cleaned_dict)):
			verString = list(filter(verFlag.match,  list(cleaned_dict)))
			assert len(verString)==1
			ver_fin = verString[0].split('=')[-1].strip()
		if any(headFlag.match(x) for x in list(cleaned_dict)):
			headFin=list(cleaned_dict.values())[0].split(',')


	if any(x in coord_sys for x in {'LPS','1'}):
		df = pd.read_csv(input_fcsv, skiprows=3, header=None)

		if df.shape[1] != 13:
			df=df.iloc[:,:14]

		df[1] = -1 * df[1] # flip orientation in x
		df[2] = -1 * df[2] # flip orientation in y

		if overwrite_fcsv:
			with open(input_fcsv, 'w') as fid:
				fid.write("# Markups fiducial file version = 4.11\n")
				fid.write("# CoordinateSystem = 0\n")
				fid.write("# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")

			df.rename(columns={0:'node_id', 1:'x', 2:'y', 3:'z', 4:'ow', 5:'ox',
								6:'oy', 7:'oz', 8:'vis', 9:'sel', 10:'lock',
								11:'label', 12:'description', 13:'associatedNodeID'}, inplace=True)

			df['associatedNodeID']= pd.Series(np.repeat('',df.shape[0]))
			df.round(6).to_csv(input_fcsv, sep=',', index=False, lineterminator="", mode='a', header=False, float_format='%.6f')

			print(f"Converted LPS to RAS: {os.path.dirname(input_fcsv)}/{os.path.basename(input_fcsv)}")
	return coord_sys,headFin

def transform_points(fcsv_path, transform):
    fiducial_points = []
    labels = []
    with open(fcsv_path, 'r') as file:
        for line in file.readlines():
            # Skip comment lines starting with '#'
            if not line.startswith('#'):
                # Extract the properties and coordinates from each line
                fields = line.strip().split(',')
                x, y, z = float(fields[1]), float(fields[2]), float(fields[3])
                fiducial_points.append([x, y, z])
                labels.append(fields[11])
    coords = np.asarray(fiducial_points)
    transformed_coords = np.zeros(coords.shape)
    transform = np.loadtxt(transform)
    transform = np.linalg.inv(transform)

    M = transform[:3, :3]
    abc = transform[:3, 3]

    for i in range(len(coords)):
        vec = coords[i, :]
        tvec = M.dot(vec) + abc
        transformed_coords[i, :] = tvec[:3]

    transformed_coords = np.round(transformed_coords, 1).astype(float)
    labels = np.array(labels).reshape(-1,1)
    
    final_array = np.hstack((transformed_coords, labels))
    
    # print(final_array.shape)
    # print(final_array)
    return final_array


def apply_warp_deformation(transform_path, fiducial_file):
    """
    Transforms fiducial points from a fiducial file using a transformation matrix.

    Parameters:
    - transform_path: str, path to the transformation matrix file
    - fiducial_file: str, path to the fiducial file

    Returns:
    - transformed_fiducial_points: list of transformed fiducial points
    - fiducial_properties: list of properties corresponding to each fiducial point
    """
    
    # Reads the transform and casts the output to a compatible format
    transform_image = sitk.ReadImage(transform_path)
    transform_image = sitk.Cast(transform_image, sitk.sitkVectorFloat64)

    # Load it as a transform
    transform = sitk.Transform(transform_image)

    # Loop through the file and extract the fiducial points
    fiducial_points = []
    labels = []
    with open(fiducial_file, 'r') as file:
        for line in file.readlines():
            # Skip comment lines starting with '#'
            if not line.startswith('#'):
                # Extract the properties and coordinates from each line
                fields = line.strip().split(',')
                x, y, z = float(fields[1]), float(fields[2]), float(fields[3])
                fiducial_points.append([x, y, z])
                labels.append(fields[11])

    fiducial_points = np.array(fiducial_points) * np.array([-1, -1, 1])

    # Apply the transform to the fiducial points
    transformed_fiducial_points = []
    for point in fiducial_points:
        transformed_point = transform.TransformPoint(point.tolist())
        transformed_fiducial_points.append(transformed_point)

    transformed_fiducial_points = np.array(transformed_fiducial_points) * np.array([-1, -1, 1])
    labels = np.array(labels).reshape(-1,1)
    # print(transformed_fiducial_points)
    # print(labels.reshape(-1,1))
    final_array = np.hstack((transformed_fiducial_points, labels))
    
    # print(final_array.shape)
    # print(final_array)
    return final_array

def array_to_fcsv(coord_array, output_fcsv):
    with open(output_fcsv, "w") as fid:
        fid.write("# Markups fiducial file version = 4.11\n")
        fid.write("# CoordinateSystem = RAS\n")
        fid.write(
            "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n"
        )

    out_df = {
        "node_id": [],
        "x": [],
        "y": [],
        "z": [],
        "ow": [],
        "ox": [],
        "oy": [],
        "oz": [],
        "vis": [],
        "sel": [],
        "lock": [],
        "label": [],
        "description": [],
        "associatedNodeID": [],
    }

    for i, coord in enumerate(coord_array):
        out_df["node_id"].append(i + 1)
        out_df["x"].append(coord[0])
        out_df["y"].append(coord[1])
        out_df["z"].append(coord[2])
        out_df["ow"].append(0)
        out_df["ox"].append(0)
        out_df["oy"].append(0)
        out_df["oz"].append(0)
        out_df["vis"].append(1)
        out_df["sel"].append(1)
        out_df["lock"].append(1)
        out_df["label"].append(coord[3])
        out_df["description"].append("elec_type")
        out_df["associatedNodeID"].append("vtkMRMLScalarVolumeNode2")

    out_df = pd.DataFrame(out_df)
    out_df.to_csv(
        output_fcsv,
        sep=",",
        index=False,
        lineterminator="",
        mode="a",
        header=False,
        float_format="%.3f",
    )

if __name__ == "__main__":
    determineFCSVCoordSystem(snakemake.input['contacts'], overwrite_fcsv=True)
    determineFCSVCoordSystem(snakemake.input['planned'], overwrite_fcsv=True)
    determineFCSVCoordSystem(snakemake.input['actual'], overwrite_fcsv=True)
    
    warp_SEEGA = apply_warp_deformation(snakemake.input['deformable_warp'], snakemake.input['contacts'])
    warp_planned = apply_warp_deformation(snakemake.input['deformable_warp'], snakemake.input['planned'])
    warp_actual = apply_warp_deformation(snakemake.input['deformable_warp'], snakemake.input['actual'])
    
    array_to_fcsv(warp_SEEGA, snakemake.output['deformable_contacts'])
    array_to_fcsv(warp_planned, snakemake.output['deformable_planned'])
    array_to_fcsv(warp_actual, snakemake.output['deformable_actual'])

    affine_SEEGA = transform_points(snakemake.input['contacts'], snakemake.input['affine'])
    affine_planned = transform_points(snakemake.input['planned'], snakemake.input['affine'])
    affine_actual = transform_points(snakemake.input['actual'], snakemake.input['affine'])
    
    array_to_fcsv(affine_SEEGA, snakemake.output['affine_contacts'])
    array_to_fcsv(affine_planned, snakemake.output['affine_planned'])
    array_to_fcsv(affine_actual, snakemake.output['affine_actual'])