# To get discharged data from HxF
# Create dataframe with Pebble vs step for each field
# Create dataframe per step

# Import modules
import os
from glob import glob
import re
import pandas as pd
import sys
import shutil 
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

#
path_output = './chunked_gammas'
path_folder_data = '/global/scratch/users/yvesrobert/HxF_dev/Cases/HTR10_restart/Data/'
step_start = 1200
step_end = 1202
saving = True
if len(sys.argv) == 3:
    step_start = int(sys.argv[1])
    step_end = int(sys.argv[2])

# Script
os.makedirs(f'{path_output}/Discharged_data', exist_ok=True)
data_files = natural_sort(glob(f'{path_folder_data}/discharged_fuel*.csv'))
data_files = [f for f in data_files if int(f.split('_')[-1].split('.csv')[0])>=step_start and int(f.split('_')[-1].split('.csv')[0])<=step_end]

data ={}
for file in data_files:
    print(file)
    step = int(file.split('_')[-1].split('.csv')[0])
    df = pd.read_csv(file, index_col=0).reset_index()
    df = df.rename(columns={df.columns[0]: "Fuel zone index"}) 
    if saving:  
        df.to_csv(f'{path_output}/Discharged_data/discharged_{step}.csv')
    
    if len(data)==0:
        for field in df.columns.drop(['Fuel zone index', 'id', 'r', 'uni', 'mat_name', 'initial', 'insertion_step', 'discarded', 'pass_nsteps', 'integrated_power_pebbles', 'pass_integrated_power_pebbles'] + [col for col in df if '_unc' in col]):
            data[field] = {}
    for field in df.columns.drop(['Fuel zone index', 'id', 'r', 'uni', 'mat_name', 'initial', 'insertion_step', 'discarded', 'pass_nsteps', 'integrated_power_pebbles', 'pass_integrated_power_pebbles'] + [col for col in df if '_unc' in col]):
        data[field][step] = df[field]

print('Compiling data...')
for field in data:
    data[field] = pd.concat(data[field], axis=1).T
    if saving:
        data[field].to_csv(f'{path_output}/Discharged_data/discharged_{field}.csv')

shutil.copy2(f'{path_folder_data}/cycle_{step_end}.csv', f'{path_output}/Discharged_data/')
print('Done')