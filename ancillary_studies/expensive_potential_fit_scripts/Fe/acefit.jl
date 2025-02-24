# ----------------------------------------------------------------------
# This script is part of the ML-MIX repository.
# 
# Copyright (2025)  Lakshmi Shenoy, Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------


using Pkg
Pkg.activate("../../../")
using Distributed
addprocs(30, exeflags="--project=$(Base.active_project())")
@everywhere using ACEpotentials
# Read data
data_file = "espresso_iron_2023.xyz"
data = read_extxyz(data_file)
@show length(data)
# Test-Train split
idx_test = 1:100:length(data) 
idx_train = [ i for i in 1:length(data) if !(i in idx_test) ]
data_test = data[idx_test]
data_train = data[idx_train]
# Define model
model = acemodel(elements = [:Fe], #, :P],
                order = 3,
                totaldegree = 20,
                rcut = 5.5,
                Eref = [:Fe => -3456.00113622806]) #, :P => -216.0991804514])
@show length(model.basis)
solver = ACEfit.BLR(factorization=:svd)
P = smoothness_prior(model; p = 4) 
weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0 , "V" => 1.0 ),
               "slice_sample" => Dict("E" => 80.0, "F" => 0.1 , "V" => 1.0),
               "prim_random" => Dict("E" => 80.0, "F" => 0.1 , "V" => 1.0))
# ACE fit
acefit!(model, data_train; solver=solver, prior = P, weights=weights)
@info("Training Error Table")
ACEpotentials.linear_errors(data_train, model) 
@info("Test Error Table")
ACEpotentials.linear_errors(data_test, model) 
export2json("jace_espresso.json", model)
export2lammps("jace_espresso.yace", model)
