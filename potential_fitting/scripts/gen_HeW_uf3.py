# ----------------------------------------------------------------------
# This script is part of the ML-MIX repository.
# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------

import json
import numpy as np
import argparse
#open json file as dict

#use argparse to get path to json file
parser = argparse.ArgumentParser()
parser.add_argument('--json_path',type=str)
parser.add_argument('--in_file_name',type=str)
parser.add_argument('--out_file_name',type=str)

args = parser.parse_args()

json_path = args.json_path
in_file_name = args.in_file_name
out_file_name = args.out_file_name

f = open(f'{json_path}/{in_file_name}.json')

data = json.load(f)

print(data.keys())
print(data['coefficients'].keys())
zeros_1b = 0.0
zeros_2b = np.zeros_like(data['coefficients']['W-W']).tolist()
zeros_3b = np.zeros_like(data['coefficients']['W-W-W']).tolist()


knots_2b = data['knots']['W-W']
knots_3b = data['knots']['W-W-W']

knots_map_2b = data['knots_map']['W-W']
knots_map_3b = data['knots_map']['W-W-W']

new_element_list = ['W','He']

He_2b = ['He-He','W-He']
He_3b = ['He-He-He','W-He-He','W-W-He','He-W-W','He-He-W']

data['coefficients']['He'] = zeros_1b

for key in He_2b:
    data['coefficients'][key] = zeros_2b
    data['knots'][key] = knots_2b
    data['knots_map'][key] = knots_map_2b

for key in He_3b:
    data['coefficients'][key] = zeros_3b
    data['knots'][key] = knots_3b
    data['knots_map'][key] = knots_map_3b

data['element_list'] = new_element_list

#delete the data coverage entry
del data['data_coverage']


with open(f'{json_path}/{out_file_name}.json','w') as f:
    json.dump(data,f,indent=4)
f.close()


