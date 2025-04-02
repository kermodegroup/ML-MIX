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
import ACE1x:_set_params!
jsonpath = "./"
lammpspath = "./"


model_name = config["output_name"]
@info("Generating placeholder ACE model $model_name.yace")
elements = Symbol.(el for el in config["fit-elements"])
e_ref = Dict(el => config["e_refs"][el] for el in config["fit-elements"])
order = config["order"]
totaldegree = config["totaldegree"]
rcut = config["rcut"]
smoothness_prior_strength = config["smoothness_prior_strength"]


@info("Setting up basis...")

model = acemodel(
    elements = elements,
    Eref = e_ref,
    order = order,
    totaldegree = totaldegree,
    rcut = rcut,
)


export2json("$jsonpath/$model_name.json", model)
export2lammps("$lammpspath/$model_name.yace", model)
println("########Exiting Julia########")