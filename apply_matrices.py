# In case detector matrix has changed, reapply matrices.

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
gamma_type_source = 'cumulated'
matrices_to_apply = [{'name': 'incident', 'Egrid_path': './detector_data/geometry_Egrid.txt', 'Egrid_scale':1, 'matrix_path': './detector_data/geometry.txt'}, # scale is important, for instance if the geometry energy grid is in keV and the detector energy grid is in Mev, a 1e-3 scaling is necessary.
                     {'name': 'absorbed', 'Egrid_path': './detector_data/NaI_Egrid.txt', 'Egrid_scale':1e-3, 'matrix_path': './detector_data/NaI_absorption.txt'},
                     {'name': 'detected', 'Egrid_path': './detector_data/NaI_Egrid.txt', 'Egrid_scale':1e-3, 'matrix_path': './detector_data/NaI_resolution.txt'}]

print('Reading energy grids and matrices')
for i in range(len(matrices_to_apply)):
    print(f'\t{matrices_to_apply[i]["name"].title()}')
    matrices_to_apply[i]['Egrid'] = np.loadtxt(matrices_to_apply[i]['Egrid_path'])*matrices_to_apply[i]['Egrid_scale']
    matrices_to_apply[i]['matrix'] = np.loadtxt(matrices_to_apply[i]['matrix_path'])
names = ['cumulated'] + [matrices_to_apply[i]['name'] for i in range(len(matrices_to_apply))]

print('Reading created source files and processing')
files = natural_sort(glob(f'{path_folder}/*/working_directory/gammas_*_{gamma_type_source}_{decay_index}.pkl'))
steps = [int(f.split('gammas_')[-1]. split('_')[0]) for f in files]
print(f'\tFound {len(files)} files from step {steps[0]} to step {steps[-1]}')

# Loop over each gamma emission file, process it, and save the results
data = [{} for i in range(len(matrices_to_apply)+1)]
for i, file in enumerate(files):
    print(f'\tStep {steps[i]}')
    with open(file, 'rb') as f:
        # Source to cumulated (binned) gammas
        data[0][steps[i]] = pd.read_pickle(f)
        print('\t\tCumulated')
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