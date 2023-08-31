# To use once create_all_gammas is done with chunks.
# Create the large dataframe with E/Pebble vs step.
# Equivalent to compile_gammas, but reapplies matrices and helps for other energy grids.

# Import required libraries
import struct
import os
from glob import glob
import re
import math
import numpy as np
import subprocess
from time import time
import pandas as pd
import matplotlib.pyplot as plt # not needed if you don't want to plot
import pickle
import warnings
import sys
from scipy.interpolate import interp2d

# Function to sort filenames or other strings containing numbers in a human-friendly way
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# Function to read gamma emission data for a particular material
def read_gamma_file(gamma_file_path, mat_name):
    with open(gamma_file_path, 'r') as file:
        lines = file.read().splitlines()

    nlines = len(lines)
    line_idx = 0
    materials = {}
    print(f'\tReading gamma file at: {gamma_file_path}')
    while line_idx < nlines:
        while f'{mat_name}z' not in lines[line_idx]:
            line_idx +=1
            if line_idx >= nlines:
                break
        if line_idx >= nlines:
            break
        zone = int(lines[line_idx].split(f'{mat_name}z')[1].split(' =')[0].split('_')[0])
        if f'{mat_name}z{zone}' not in materials:
            line_idx +=1
            table = []
            while '];' not in lines[line_idx]:
                if lines[line_idx] and lines[line_idx][0] != '%':
                    table.append([int(lines[line_idx].split()[0])] + [float(i) for i in lines[line_idx].split()[1:]])
                line_idx +=1
            materials[f'{mat_name}z{zone}'] = pd.DataFrame(table, columns = ['ZAI', 'Specific intensity', 'Total emission rate', 'Cumulative material fraction', 'Energy', 'Relative intensity', 'Cumulative nuclide total'])
            materials[f'{mat_name}z{zone}']['Energy emission rate'] = materials[f'{mat_name}z{zone}']['Relative intensity']/materials[f'{mat_name}z{zone}']['Specific intensity']*materials[f'{mat_name}z{zone}']['Total emission rate']
            materials[f'{mat_name}z{zone}'] = materials[f'{mat_name}z{zone}'].loc[:, ['ZAI', 'Energy', 'Energy emission rate']]
            materials[f'{mat_name}z{zone}'].name = f'{mat_name}z{zone}'
        line_idx +=1

    print(f'\t\t{len(materials)} materials treated')
    return materials

# Function to aggregate gamma data, optionally using an energy resolution for binning
def cumulate_gamma(material_df, resolution=None):
    warnings.filterwarnings("ignore", message=".*Index.ravel returning ndarray is deprecated.*") # Will have to find a solution at some point (will be deprecated)
    #print('\tCumulating gamma counts for same energy')
    # Do not bin
    if isinstance(resolution, type(None)):
        new_table = material_df.groupby('Energy').agg({"ZAI": list, "Energy emission rate": "sum"})
    # Bin
    else:
        try:
            _ = len(resolution)
            bins = resolution
        except:
            bins = np.arange(material_df['Energy'].min(), material_df['Energy'].max()+resolution/2, resolution)
        cut_data = pd.cut(material_df['Energy'], bins=bins).apply(lambda x: x.mid)
        new_table = material_df.groupby(cut_data).agg({"ZAI": list, "Energy emission rate": "sum"})
    return new_table.sort_index()

# Function to reprocess gamma data using interpolation and a treatment matrix
def treat_gammas(df_gammas, new_Egrid, treatment_matrix):
    warnings.filterwarnings("ignore", message=".*`interp2d` is deprecated*") # Will have to find a solution at some point (will be deprecated)
    Egrid_centers = np.array(df_gammas.index)
    new_Egrid_centers = (new_Egrid[:-1] + new_Egrid[1:])/2
    treatment_interp = interp2d(new_Egrid_centers, new_Egrid_centers, treatment_matrix, kind='linear')
    treatment = treatment_interp(Egrid_centers, Egrid_centers)
    treated_gammas = treatment @ df_gammas
    treated_gammas.index = Egrid_centers
    return treated_gammas

path_folder = '/global/scratch/users/co_nuclear/gamma_maker/chunked_gammas/' # directory containing the chunked gamma data
decay_index = 0 # index of the decay step (suffix of the files)
matrices_to_apply = [{'name': 'incident', 'Egrid_path': './detector_data/geometry_Egrid.txt', 'Egrid_scale':1, 'matrix_path': './detector_data/geometry.txt'},
                     {'name': 'absorbed', 'Egrid_path': './detector_data/NaI_Egrid.txt', 'Egrid_scale':1e-3, 'matrix_path': './detector_data/NaI_absorption.txt'},
                     {'name': 'detected', 'Egrid_path': './detector_data/NaI_Egrid.txt', 'Egrid_scale':1e-3, 'matrix_path': './detector_data/NaI_resolution.txt'}]

print('Reading energy grids and matrices')
for i in range(len(matrices_to_apply)):
    print(f'\t{matrices_to_apply[i]["name"].title()}')
    matrices_to_apply[i]['Egrid'] = np.loadtxt(matrices_to_apply[i]['Egrid_path'])*matrices_to_apply[i]['Egrid_scale']
    matrices_to_apply[i]['matrix'] = np.loadtxt(matrices_to_apply[i]['matrix_path'])
names = ['cumulated'] + [matrices_to_apply[i]['name'] for i in range(len(matrices_to_apply))]

print('Reading created gamma files and processing')
files = glob(f'{path_folder}/*/working_directory/gammas_*_{decay_index}.pkl')
files = natural_sort([f for f in files if re.search(r'gammas_([0-9]+)_' + str(decay_index) + r'\.pkl', f)])
steps = [int(f.split('gammas_')[-1]. split('_')[0]) for f in files]
print(f'\tFound {len(files)} files from step {steps[0]} to step {steps[-1]}')

# Loop over each gamma emission file, process it, and save the results
data = [{} for i in range(len(matrices_to_apply)+1)]
for i, file in enumerate(files):
    print(f'\tStep {steps[i]}')
    with open(file, 'rb') as f:
        # Source to cumulated (binned) gammas
        materials = pickle.load(f)[0]
        print('\t\tCumulated')
        data[0][steps[i]] = pd.DataFrame({i:cumulate_gamma(materials[list(materials.keys())[i]], resolution=matrices_to_apply[0]['Egrid'])['Energy emission rate'] for i in range(len(materials))})

        # Apply matrices
        for j in range(len(matrices_to_apply)):
            print(f'\t\t{matrices_to_apply[j]["name"].title()}')
            data[j+1][steps[i]] = treat_gammas(data[j][steps[i]], matrices_to_apply[j]['Egrid'], matrices_to_apply[j]['matrix'])
    
    # Process
    for j in range(len(data)):
        data[j][steps[i]] = data[j][steps[i]].stack()

# Gather all processed data into a single data frame and save it
print('Gathering data and saving...')
compiled_data = {}
for i in range(len(data)):
    print(f'\t{names[i].title()}')
    compiled_data[names[i]] = pd.DataFrame(data[i]).T
    compiled_data[names[i]].index.name = 'Step'
    compiled_data[names[i]].columns.names = ['Energy', 'Pebble']
    
    print('\t\tSaving')
    compiled_data[names[i]].to_pickle(f'{path_folder}/{prefix}{names[i]}_{decay_index}.pkl', 'wb')