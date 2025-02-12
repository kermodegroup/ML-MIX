#!/bin/bash

# script path
script_path="../../scripts"

# Model parameters
alpha="1e-7"
elements="Si"
e_ref="-158.54496821"
order="3"
totaldegree="15"
rcut="6.0"
model_name1="model_unconstrained_alpha_${alpha}_${order}_${totaldegree}"
model_name2="model_constrained_alpha_${alpha}_${order}_${totaldegree}"
nprocs=40

# Name of folder for this set of potentials
folder_name="energy_diff_constraints_alpha_${alpha}_${order}_${totaldegree}"

# Create the folder 'potentials' if it doesn't exist
mkdir -p cheap_potentials

# Create the main folder
mkdir "cheap_potentials/$folder_name"

# Create the subfolders
mkdir "cheap_potentials/$folder_name/raw_pot"
mkdir "cheap_potentials/$folder_name/lammps_pot"
mkdir "cheap_potentials/$folder_name/analysis"

# echo command
echo "calling Julia for ACE fit..."

julia "$script_path/fit_ACE.jl" \
    --raw_output "cheap_potentials/$folder_name/raw_pot" \
    --lammps_output "cheap_potentials/$folder_name/lammps_pot" \
    --analysis_path "cheap_potentials/$folder_name/analysis" \
    --model_name_unconstrained "$model_name1" \
    --model_name_constrained "$model_name2" \
    --alpha "$alpha" \
    --nprocs "$nprocs" \
    --elements "$elements" \
    --e_ref "$e_ref" \
    --order "$order" \
    --totaldegree "$totaldegree" \
    --rcut "$rcut"

echo Testing elastic constants and lattice parameter...

python "$script_path/test_elastic_constants.py" "cheap_potentials/$folder_name/lammps_pot" "cheap_potentials/$folder_name/analysis" "cheap_potentials/$folder_name" "$model_name1" "$model_name2"