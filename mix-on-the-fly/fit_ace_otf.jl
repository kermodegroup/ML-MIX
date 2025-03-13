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

println("########Entering Julia########")
using Pkg
# path to ML-MIX is first sys arg
# path to config.yaml is second sys arg
ml_mix_path = ARGS[1]
Pkg.activate(ml_mix_path) # activate the environment
using YAML
config_path = ARGS[2]
config = YAML.load_file(config_path)
using ACEpotentials
jsonpath = "./"
lammpspath = "./"

@info("Loading data...")
data = read_extxyz("otf-fit-data.xyz")
if length(data) == 0
    println("No data found in file")
    exit(1)
end

model_name = config["output_name"]
elements = Symbol.(el for el in config["fit-elements"])
e_ref = Dict(el => config["e_refs"][el] for el in config["fit-elements"])
order = config["order"]
totaldegree = config["totaldegree"]
rcut = config["rcut"]
smoothness_prior_strength = config["smoothness_prior_strength"]
energy_key = nothing
force_key = "fit_forces"
virial_key = nothing

if "pae" in keys(data[1].data)
    pae_key = "pae"
else
    pae_key = nothing
end

if "mask" in keys(data[1].data)
    mask_key = "mask"
else
    mask_key = nothing
end


@info("Setting up basis...")

model = acemodel(
    elements = elements,
    Eref = e_ref,
    order = order,
    totaldegree = totaldegree,
    rcut = rcut,
)
@show length(model.basis)


@info("Fitting model...")
# fit reference unconstrained ACE model to As data
solver = ACEfit.BLR(factorization=:svd)
P = smoothness_prior(model; p = smoothness_prior_strength) #this is our Tikhonov regularizer

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0 , "V" => 1.0, "PAE" => 1.0))
acefit!(model, 
        data; 
        solver=solver, 
        prior = P, 
        weights=weights,
        energy_key = energy_key,
        force_key = force_key,
        virial_key = virial_key,
        pae_key = pae_key,
        mask_key = mask_key)

export2json("$jsonpath/$model_name.json", model)
export2lammps("$lammpspath/$model_name.yace", model)

@info("Evaluating training error...")
rmse, mae = ACEpotentials.linear_errors(data, 
                                        model, 
                                        energy_key=energy_key, 
                                        force_key=force_key,
                                        virial_key=virial_key,
                                        pae_key=pae_key,
                                        mask_key=mask_key)

rmses = rmse.second["set"]["E"]
force_RMSE = rmse.second["set"]["F"]
if pae_key != nothing
    pae_RMSE = rmse.second["set"]["PAE"]
end


#write errors to file
open("otf-RMSE.txt", "w") do f
    write(f, "Force train RMSE: $force_RMSE\n")
    if pae_key != nothing
        write(f, "PAE RMSE: $unconstrained_pae_err\n")
    end
end

println("########Exiting Julia########")

