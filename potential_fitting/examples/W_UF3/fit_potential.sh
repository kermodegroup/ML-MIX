#!/bin/bash

# script path
script_path="../../scripts"

# User name - to associate with fitted potential
user_name="Joe_Bloggs"

# constrained fit parameters
alpha="1e-10"
max_strain="0.005"
lambda_max=1000000
lambda_points=1000

# model names
folder_name="energy_diff_constraints_e_${max_strain}_alpha_${alpha}_lambda_max_${lambda_max}"
model_name1="model_unconstrained"
model_name2="model_constrained_alpha_${alpha}"

# Model parameters
element="W"
rcut2b=6.0
rcut3b=3.5
resolution2b=15
resolution3b=6
rcutmin2b=0.001
rcutmin3b=1.5
energy_key="energy"
ridge1b=0.0
ridge2b=0.0
ridge3b=1e-8
curvature2b=1e-8
curvature3b=0.0

# Create the folder 'potentials' if it doesn't exist
mkdir -p cheap_potentials

# Create the main folder
mkdir "cheap_potentials/$folder_name"

# Create the subfolders
mkdir "cheap_potentials/$folder_name/raw_pot"
mkdir "cheap_potentials/$folder_name/lammps_pot"
mkdir "cheap_potentials/$folder_name/analysis"

# echo command
echo "calling python for UF3 fit..."

# Run fit_UF3_to_ACE.py with user-specified model names and additional parameters
python $script_path/fit_UF3.py \
    --raw_output "cheap_potentials/$folder_name/raw_pot" \
    --lammps_output "cheap_potentials/$folder_name/lammps_pot" \
    --analysis_path "cheap_potentials/$folder_name/analysis" \
    --model_name_unconstrained "$model_name1" \
    --model_name_constrained "$model_name2" \
    --alpha "$alpha" \
    --element "$element" \
    --rcut2b "$rcut2b" \
    --rcut3b "$rcut3b" \
    --resolution2b "$resolution2b" \
    --resolution3b "$resolution3b" \
    --rcutmin2b "$rcutmin2b" \
    --rcutmin3b "$rcutmin3b" \
    --energy_key "$energy_key" \
    --lambda_max "$lambda_max" \
    --ridge1b "$ridge1b" \
    --ridge2b "$ridge2b" \
    --ridge3b "$ridge3b" \
    --curvature2b "$curvature2b" \
    --curvature3b "$curvature3b" \
    --lambda_points "$lambda_points"

# Get all files matching the user-specified model names in UF3_raw
files=$(find "cheap_potentials/$folder_name/raw_pot" -name "*$model_name1*" -o -name "*$model_name2*")

# Run make_lammps.py for each file and rename the output file
model_paths=()
for file in $files; do
    filename=$(basename "$file")
    model_name="${filename%.*}"
    #echo command
    echo Generating LAMMPS potentials...

    python ../../../external/uf3/lammps_plugin/scripts/generate_uf3_lammps_pots.py \
    -a $user_name \
    -u metal \
    -m "$file" \
    -d "cheap_potentials/$folder_name/lammps_pot"
    mv "cheap_potentials/$folder_name/lammps_pot/W.uf3" "cheap_potentials/$folder_name/lammps_pot/$model_name.uf3"
done

# now call the script test_elastic_constants.py, 
# passing the path to UF3_lammps, the path to the elastic_constants folder,
# and the model names as arguments

echo testing elastic constants and lattice parameter...

python "$script_path/test_elastic_constants.py" "cheap_potentials/$folder_name/lammps_pot" "cheap_potentials/$folder_name/analysis" "cheap_potentials/$folder_name" "$model_name1" "$model_name2"

