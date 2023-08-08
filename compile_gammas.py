# To use once create_all_gammas is done with chunks.
# Create the large dataframe with E/Pebble vs step

import pickle
import pandas as pd
import numpy as np
from glob import glob
import re
import sys

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

path_folder = './chunked_gammas'
prefix = 'gammas_'
decay_index = 0

for gamma_type in ['cumulated', 'incident', 'absorbed', 'detected']:
    files = natural_sort(glob(f'{path_folder}/{prefix}*/*/*_{gamma_type}_{decay_index}.pkl'))
    if len(files)==0:
        raise Exception(f'No file found as: {path_folder}/{prefix}*/*_{gamma_type}_{decay_index}.pkl')
    data = None
    for file in files:
        step = int(file.split(f'_{gamma_type}_')[0].split('_')[-1])
        print(step, file)
        with open(file, "rb") as f:
            materials = pickle.load(f)
        if isinstance(data, type(None)):
            series = materials.stack()
            data = pd.DataFrame(series, columns=[step])
        else:
            data[step] = materials.stack()
    
    print('Gathering data...')
    data = pd.DataFrame(data).T
    data.index.name = 'Step'
    data.columns.names = ['Energy', 'Pebble']
    
    print('Storing data')
    data.to_pickle(f'{path_folder}/{prefix}{gamma_type}_{decay_index}.pkl')
    print('Done')