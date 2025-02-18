## This is the params.py file for the plotting script

import os
investigation = 'ref_D_coeff'
T = 800

path = f'results/{investigation}/{T}K'
output_path = f'plots/{investigation}/{T}K'
if not os.path.exists(output_path):
    os.makedirs(output_path)

dump_name_list = []

for file in os.listdir(path):
    if file.endswith(".lammpstrj"):
        dump_name_list.append(f'{path}/{file}')



#this is the number of non-overlapping timelag samples that is allowed to contribute
#to a mean squared displacement calculation at that time lag
min_sample_number = 60
#this is the minimum number of dump frames to start sampling at, should be roughly twice the
#mean jump wait time.
min_sample_rate = 6

# this is how often the HeDump files were written in the simulation
dump_freq=10

# this is the size of the cell
n=16