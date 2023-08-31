# Use create_all_gammas.sh instead of that directly

# Import modules
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

# Function for natural sorting of a list
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# Function to execute a bash command and print the output
def execute_bash_command(command, printing=True, infile=None):
    print(f'\tRunning {command}')
    if printing and infile:
        print(f'\t\tOutput will be written in {infile}')
        with open(infile, 'w') as f:
            pass
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    ## Input
    # Read the command output line by line and print it live
    print('\t\tWaiting for the job to end')
    while True:
        output = process.stdout.readline().decode('utf-8')
        if output and printing:
            if not infile:
                print(output.strip())
            else:
                with open(infile, 'a') as f:
                    f.write(output)
        if output == '' and process.poll() is not None:
            break    
    # Wait for the command to finish and get the return code
    return_code = process.wait()
    
    return return_code

# Function to read a restart file and extract material data. Can put any ZAI in ZAI_fields, or one/multiple of the following in nonZAI_fiels: n, name, bu_global, bu_days, nnuc, adens, mdens, burnup 
def read_restart(list_paths, step=0, mat_parent='fuel', nonZAI_fields=['burnup'], ZAI_fields=[], df=True):
    if isinstance(list_paths, str):
        list_paths = [list_paths]
    elif len(list_paths)==0:
        raise Exception('No restart file given')
    materials = dict()
    ZAI_indices = []
    for path_to_file in list_paths:
        print('\t',path_to_file)
        current_step = -1
        burnups = dict()
        # Read restart file
        with open(path_to_file, mode='rb') as file:  # b is important -> binary
            while True:
                # Read material block
                s = file.read(8)
                if not s:
                    break
                material = {}
                material['n'] = struct.unpack("q", s)[0]  # length of material name
                material['name'] = struct.unpack("{}s".format(material['n']), file.read(material['n']))[0].decode('UTF-8') # material name
                material['bu_global'] = struct.unpack("d", file.read(8))[0] # BU of snapshot
                material['bu_days'] = struct.unpack("d", file.read(8))[0] # time of snapshot
                material['nnuc'] = struct.unpack("q", file.read(8))[0] # Number of nuclides in material
                material['adens'] = struct.unpack("d", file.read(8))[0] # Atomic density of material
                material['mdens'] = struct.unpack("d", file.read(8))[0] # Mass density of material
                material['burnup'] = struct.unpack("d", file.read(8))[0] # Burnup of material
                if len(burnups) == 0 or material['bu_global'] != burnups[list(burnups.keys())[-1]]:
                    current_step += 1
                    burnups[current_step] = material['bu_global']

                # Check if material name matches
                if material['name'][:min(len(mat_parent), len(material['name']))+1] != f'{mat_parent}z' or current_step != step:
                    # Seek to the next block by calculating the number of bytes to skip
                    bytes_to_skip = 16 * material['nnuc']  # Size of the data block (16 bytes for each nuclide)
                    file.seek(bytes_to_skip, 1)  # Move the file pointer forward by the calculated number of bytes
                    continue

                materials[material['name']] = {field: material[field] for field in nonZAI_fields}
                if len(ZAI_fields)==0:
                    # Seek to the next block by calculating the number of bytes to skip
                    file.seek(16 * material['nnuc'], 1)  # Move the file pointer forward by the calculated number of bytes
                    continue

                # Just once
                if len(ZAI_indices) == 0:
                    adens_list = []
                    ZAI_list = []
                    for i in range(material['nnuc']):
                        ZAI, adens = struct.unpack("qd", file.read(16))
                        ZAI_list.append(ZAI)
                        adens_list.append(adens)
                    ZAI_indices = [ZAI_list.index(int(ZAI)) for ZAI in ZAI_fields]
                    ZAI_indices, ZAI_fields = zip(*sorted(zip(ZAI_indices, ZAI_fields)))
                    for i, index in enumerate(ZAI_indices):
                        materials[material['name']][int(ZAI_fields[i])] = adens_list[index]
                # The rest of the cases
                else:
                    last_index = 0
                    for i in range(len(ZAI_fields)):
                        index = ZAI_indices[i]
                        file.seek(16*(index-last_index), 1)
                        ZAI, adens = struct.unpack("qd", file.read(16))
                        materials[material['name']][int(ZAI)] = adens
                        last_index = int(index)+1
                    file.seek(16*(material['nnuc']-last_index), 1)
    # Sort and convert materials information
    if df:
        materials = pd.DataFrame.from_dict(materials, orient='index')
        materials = materials.loc[natural_sort(materials.index)]
    else:
        materials = {key: materials[key] for key in natural_sort(materials)}
    return materials

# Function to read a restart file and extract one material data.
def read_restart_material(list_paths, step=0, mat_name='fuel', df=True):
    found = False
    if isinstance(list_paths, str):
        list_paths = [list_paths]
    for path_to_file in list_paths:
        if found:
            break
        print('\t', path_to_file)
        current_step = -1
        burnups = dict()
        # Read restart file
        with open(path_to_file, mode='rb') as file:  # b is important -> binary
            while True:
                # Read material block
                s = file.read(8)
                if not s:
                    break
                material = {}
                material['n'] = struct.unpack("q", s)[0]  # length of material name
                material['name'] = struct.unpack("{}s".format(material['n']), file.read(material['n']))[0].decode('UTF-8') # material name
                material['bu_global'] = struct.unpack("d", file.read(8))[0] # BU of snapshot
                material['bu_days'] = struct.unpack("d", file.read(8))[0] # time of snapshot
                material['nnuc'] = struct.unpack("q", file.read(8))[0] # Number of nuclides in material
                material['adens'] = struct.unpack("d", file.read(8))[0] # Atomic density of material
                material['mdens'] = struct.unpack("d", file.read(8))[0] # Mass density of material
                material['burnup'] = struct.unpack("d", file.read(8))[0] # Burnup of material
                if len(burnups) == 0 or material['bu_global'] != burnups[list(burnups.keys())[-1]]:
                    current_step += 1
                    burnups[current_step] = material['bu_global']

                # Check if material name matches
                if material['name'] != mat_name or current_step != step:
                    # Seek to the next block by calculating the number of bytes to skip
                    bytes_to_skip = 16 * material['nnuc']  # Size of the data block (16 bytes for each nuclide)
                    file.seek(bytes_to_skip, 1)  # Move the file pointer forward by the calculated number of bytes
                else:
                    found = True
                    for i in range(material['nnuc']):
                        ZAI, adens = struct.unpack("qd", file.read(16))
                        material[ZAI] = adens
                    break
    if found:
        # Sort and convert materials information
        if df:
            material = pd.DataFrame(material, index=[mat_name])
        else:
            material = {key: material[key] for key in natural_sort(material)}
        return material
    else:
        raise Exception(f'Did not find material {mat_name}')


# Function to read a gamma file and extract material data
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

# Function to cumulate gamma data based on energy. If a resolution is given, bin the energies accordingly
def cumulate_gamma(material_df, resolution=None):
    warnings.filterwarnings("ignore", message=".*Index.ravel returning ndarray is deprecated.*") # Will have to find a solution at some point (will be deprecated)
    print('\tCumulating gamma counts for same energy')
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


def decay_materials(decay_steps, prefix_restart, divided_mat_name, r_trisos, N_trisos, sss_exe, acefile, decfile, nfyfile, step=0, path_output='./', sample_pebbles=np.inf, ncores=None, print_Serpent=True):
    # Process input
    restart_files = natural_sort(glob(f'{prefix_restart}*'))
    N_files = len(restart_files)
    V_per_material = 4/3*math.pi*r_trisos**3*N_trisos
    decay_steps = " ".join(np.atleast_1d(decay_steps).astype(str))
    if not ncores:
        import multiprocessing
        ncores = multiprocessing.cpu_count()
    
    # Read restart 
    materials = read_restart(restart_files, step, divided_mat_name, ['name'], [])
    if len(materials) > sample_pebbles:
        materials = materials.sample(n=sample_pebbles)
    materials = materials.sort_index()
    N_materials = len(materials)
    print(f'\tNumber of materials: {N_materials}')
    indices = materials['name'].str.split(f'{divided_mat_name}z').str[1].astype(int).values - 1 # remove 1 because first index is 1
    N_max_index = max(indices) # get the maximum material zone index

    # Create pebbles bed file
    mat_uni = f'0.0 0.0 0.0 0.0 u_mat'
    else_uni = f'0.0 0.0 0.0 0.0 u_else'
    pbed = np.full(N_max_index+1, else_uni)
    pbed[indices] = mat_uni
    np.savetxt(f"{path_output}/pebbles.inp", pbed, fmt="%s")

    path_restart = f'decay.wrk'

    # Write Serpent input for decay if needed
    # Input
    with open(f'{path_output}/input_decay.txt', 'w') as f:
        f.write(f'''
mat {divided_mat_name}   -1 burn 1 vol {N_materials*V_per_material}
92235.12c 1
92238.12c 1 
6000.12c 1
8016.12c 1

% Dummy geometry, does not matter as long as u_mat and u_else are in it
surf s_mat sph 0 0 0 1
cell c_mat u_mat {divided_mat_name} -s_mat
cell c_else u_else void -s_mat

surf infinity inf
cell c geom void -infinity

% Import pebble bed with the right indices
pbed u_pb geom "pebbles.inp"
surf s1 sqc 0.0 0.0 0.0 100.0
cell int 0 fill u_pb -s1
cell ext 0 outside s1

% Divide pebbles in zones
div {divided_mat_name} peb u_pb 1 u_mat

% Nuclear libraries
set acelib "{acefile}"
set declib "{decfile}"
set nfylib "{nfyfile}"

% Read restart
set rfr idx {step} "{prefix_restart}" {N_files} % indicate how many files there are (for decomposed cases)

% Write new restart
set rfw 1 "{path_output}/{path_restart}"

% Decay calculation
dep decstep {decay_steps}

% Depletion/transport (does not matter, no-run)
set pop 1000 20 20
set power 1
set opti 1''')

    # Run serpent
    t0 = time()
    command = f'{sss_exe} {path_output}/input_decay.txt -omp {ncores} -norun'
    _ = execute_bash_command(command, printing=print_Serpent, infile=f'{path_output}/decay_log.txt')
    t = time()
    print(f'\t{t-t0:.0f} seconds elapsed (~{(t-t0)/N_materials:.2f} seconds/material)')    
    return f"{path_output}/{path_restart}"
    
def make_gammas(prefix_restart, divided_mat_name, r_trisos, N_trisos, sss_exe, acefile, acefile_photons, decfile, nfyfile, pdatafile, step=0, path_output='./', prefix_output='gammas', sample_pebbles=np.inf, ncores=None, print_Serpent=True):
    # Write Serpent input for photon source
    # Process input
    restart_files = natural_sort(glob(f'{prefix_restart}*'))
    N_files = len(restart_files)
    V_per_material = 4/3*math.pi*r_trisos**3*N_trisos
    if not ncores:
        import multiprocessing
        ncores = multiprocessing.cpu_count()
    
    # Read restart 
    materials = read_restart(restart_files, 0, divided_mat_name, ['name'], [])
    if len(materials) > sample_pebbles:
        materials = materials.sample(n=sample_pebbles)
    materials = materials.sort_index()
    N_materials = len(materials)
    print(f'\tNumber of materials: {N_materials}')
    indices = materials['name'].str.split(f'{divided_mat_name}z').str[1].astype(int).values - 1 # remove 1 because first index is 1
    N_max_index = max(indices) # get the maximum material zone index

    # Create pebbles bed file
    mat_uni = f'0.0 0.0 0.0 0.0 u_mat'
    else_uni = f'0.0 0.0 0.0 0.0 u_else'
    pbed = np.full(N_max_index+1, else_uni)
    pbed[indices] = mat_uni
    np.savetxt(f"{path_output}/pebbles.inp", pbed, fmt="%s")

    # Input
    with open(f'{path_output}/input_gamma_{step}.txt', 'w') as f:
        f.write(f'''
mat {divided_mat_name}   -1 vol {N_materials*V_per_material}
92235.12c 1
92238.12c 1 
6000.12c 1
8016.12 1

% Dummy geometry, does not matter as long as u_mat and u_else are in it
surf s_mat sph 0 0 0 1
cell c_mat u_mat {divided_mat_name} -s_mat
cell c_else u_else void -s_mat

surf infinity inf
cell c geom void -infinity

% Import pebble bed with the right indices
pbed u_pb geom "pebbles.inp"
surf s1 sqc 0.0 0.0 0.0 100.0
cell int 0 fill u_pb -s1
cell ext 0 outside s1

% Divide pebbles in zones
div {divided_mat_name} peb u_pb 1 u_mat

% Nuclear libraries
set acelib "{acefile}"  "{acefile_photons}"
set declib "{decfile}"
set nfylib "{nfyfile}"
set pdatadir "{pdatafile}"

% Read restart
set rfr idx {step} "{prefix_restart}" {N_files} % indicate how many files there are (for decomposed cases)

% Photon calculation (does not matter, no-run)
src 1 g sg -1 1
set nps 1000
set opti 1''')

    # Run serpent
    t0 = time()
    command = f'{sss_exe} {path_output}/input_gamma_{step}.txt -omp {ncores} -norun'
    _ = execute_bash_command(command, printing=print_Serpent, infile=f'{path_output}/gamma_{step}_log.txt')
    t = time()
    print(f'\t{t-t0:.0f} seconds elapsed (~{(t-t0)/N_materials:.2f} seconds/material)')
    
    # Treating
    print('Reading created gamma file and processing')
    t0 = time()
    materials = read_gamma_file(f'input_gamma_{step}.txt_gsrc.m', divided_mat_name)
    cumulated_gammas = pd.DataFrame({i:cumulate_gamma(materials[list(materials.keys())[i]], resolution=Egrid)['Energy emission rate'] for i in range(len(materials))})    
    
    # Write binary file for faster data recovery
    print('Storing data')
    with open(f'{path_output}/{prefix_output}_{step}.pkl', "wb") as f:
        pickle.dump((materials, ), f)
    with open(f'{path_output}/{prefix_output}_cumulated_{step}.pkl', "wb") as f:
        pickle.dump(cumulated_gammas, f)
    t = time()
    print(f'\t{t-t0:.0f} seconds elapsed (~{(t-t0)/N_materials:.2f} seconds/material)')
    return materials, cumulated_gammas

def treat_gammas(df_gammas, new_Egrid, treatment_matrix, treatment_suffix, path_output='./', prefix_output='gammas'):
    warnings.filterwarnings("ignore", message=".*`interp2d` is deprecated*") # Will have to find a solution at some point (will be deprecated)
    Egrid_centers = np.array(df_gammas.index)
    new_Egrid_centers = (new_Egrid[:-1] + new_Egrid[1:])/2
    treatment_interp = interp2d(new_Egrid_centers, new_Egrid_centers, treatment_matrix, kind='linear')
    treatment = treatment_interp(Egrid_centers, Egrid_centers)
    treated_gammas = treatment @ df_gammas
    treated_gammas.index = Egrid_centers
    with open(f'{path_output}/{prefix_output}_{treatment_suffix}.pkl', "wb") as f:
        pickle.dump(treated_gammas, f)    
    return treated_gammas

if __name__=='__main__':
    ############################### INPUT ###############################

    path_folder_restart = '/global/scratch/users/yvesrobert/HxF_dev/Cases/HTR10_restart_P1T09/wrk_Serpent/'
    step_start = 1436
    step_end = 2040

    input_file_name = 'input.inp'


    # If using argument, change step start and step end based on argument
    if len(sys.argv) == 4:
        path_folder_restart = sys.argv[1]
        step_start = int(sys.argv[2])
        step_end = int(sys.argv[3])
    print(f'From {step_start} to {step_end}')
    steps = range(step_start, step_end+1)
    path_output = './'

    # To import restart data
    divided_mat_name = 'fuel' # name of the fuel material used in HxF

    # Fuel volume calculation
    r_trisos = 0.025
    N_trisos = 8335

    # Nuclear data paths
    acefile = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7u.xsdata"
    acefile_photons = "/global/home/groups/co_nuclear/serpent_photon_data/mcplib.xsdata"
    decfile = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7.dec"
    nfyfile = "/global/home/groups/co_nuclear/serpent/xsdata/endfb7/sss_endfb7.nfy"
    pdatafile = "/global/home/groups/co_nuclear/serpent_photon_data/photon_data"

    # Serpent executable
    sss_exe = '/global/home/groups/co_nuclear/HxF_tools/serpent2.2.0_HxF/sss2' # Do not change version for now

    # Others
    N_pebbles_max = np.inf # as calculations are heavy, indicate the maximum number of pebbles to extract
    decay_steps = [3] # days of decay before gamma emission
    decays_diff = np.diff([0]  + list(np.atleast_1d(decay_steps)))
    print_Serpent = True

    # Gamma processing
    Egrid = np.loadtxt('./detector_data/geometry_Egrid.txt')
    geometry_matrix = np.loadtxt('./detector_data/geometry.txt')
    Egrid_detector = np.loadtxt('./detector_data/NaI_Egrid.txt')*1e-3 # keV -> MeV
    absorption_detector_fine = np.loadtxt('./detector_data/NaI_absorption.txt'),
    resolution_detector_fine = np.loadtxt('./detector_data/NaI_resolution.txt')

    #####################################################################

    os.makedirs('working_directory', exist_ok=True)
    os.chdir('working_directory')

    # Set up Serpent library paths
    os.environ['LD_LIBRARY_PATH'] = "/global/home/groups/co_nuclear/gd_library/code/lib"
    os.environ['LIBRARY_PATH'] = "/global/home/groups/co_nuclear/gd_library/code/lib"
    os.environ['CPATH'] = "/global/home/groups/co_nuclear/gd_library/code/include"

    # Script
    for step in steps:
        print(step)
        prefix_restart = f'{path_folder_restart}/{input_file_name}.wrk_{2000000+step}' # path to the restart files (if using dd, just put prefix, 2000000 is the prefix for discharged data)
        if not isinstance(decays_diff, type(None)) and len(decays_diff)>0:
            print('\tDecaying materials')
            t0 = time()
            path_restart = decay_materials(decays_diff, prefix_restart, divided_mat_name, r_trisos, N_trisos, sss_exe, acefile, decfile, nfyfile, step=0, path_output=path_output, sample_pebbles=N_pebbles_max, print_Serpent=print_Serpent)
        else:
            path_restart = prefix_restart

        for i in range(len(np.atleast_1d(decay_steps))):
            print(f'\tMaking gammas for step {i} ({decay_steps[i]} days)')
            materials, cumulated_gammas = make_gammas(path_restart, divided_mat_name, r_trisos, N_trisos, sss_exe, acefile, acefile_photons, decfile, nfyfile, pdatafile, step=i, path_output=path_output, prefix_output=f'gammas_{step}', sample_pebbles=N_pebbles_max, print_Serpent=print_Serpent)
            
            print(f'\tTreating gammas for step {i}')
            incident_gammas = treat_gammas(cumulated_gammas, Egrid, geometry_matrix, treatment_suffix=f'{step}_incident_{i}', path_output='./', prefix_output='gammas')
            absorbed_gammas = treat_gammas(incident_gammas, Egrid_detector, absorption_detector_fine,   treatment_suffix=f'{step}_absorbed_{i}', path_output='./', prefix_output='gammas')
            detected_gammas = treat_gammas(absorbed_gammas, Egrid_detector, resolution_detector_fine,   treatment_suffix=f'{step}_detected_{i}', path_output='./', prefix_output='gammas')

        #print(f'{time()-t0:.0f} seconds elapsed (~{(time()-t0)/cumulated_gammas.shape[1]/len(decays_diff):.2f} seconds/material.decay_step)')
