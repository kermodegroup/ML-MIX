import os
import sys
T = sys.argv[1]
sim = sys.argv[2]

print(f'T: {T}')
print(f'sim: {sim}')

num_traj = 5
path = f'results/{T}/{sim}'

#get the list of dump files at path
dumps = os.listdir(path)
#create an array with the full path to each dump file
dump_name_list = [f'{path}/{dump_name}' for dump_name in dumps]
#dump_name_list = [f'results/{pot}/{run}/dump{i}.lammpstrj' for i in range(1,num_traj+1)]
#delete the second entry
# del dump_name_list[1]

dump_step = 2000
min_sample_number = 10
min_sample_rate = 2
n = 16

output_path = f'plots/many_single_Ds/{T}/{sim}'
#output_path = f'plots/espresso_avg_over_all_Ds/{T}/{sim}'
if not os.path.exists(output_path):
    os.makedirs(output_path)