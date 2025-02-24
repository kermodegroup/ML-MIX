# ----------------------------------------------------------------------
# This script is part of the ML-MIX repository.
# 
# Copyright (2025)  Matthew Nutter, Fraser Birks
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
addprocs(24, exeflags="--project=$(Base.active_project())")
using ACEpotentials
import Random
using LinearAlgebra: norm, Diagonal
dataset = read_extxyz("W_He_ACE_dataset_relabelled.xyz")
config_types = [at.data["config_type"].data for at in dataset]
function count_configs(config_types)
    config_counts = [sum(config_types.==ct) for ct in unique(config_types)]
    config_dict = Dict([ct=>cc for (ct,cc) in zip(unique(config_types), config_counts)])
end;
println("There are ", length(unique(config_types)), " unique config_types "*
        "in the dataset:")
display(count_configs(config_types))
#setup weights and solver
weights = Dict("test1" => Dict("E" => 4000.0, "F" => 100.0 , "V" => 200.0 ),
"test12" => Dict("E" => 288.6751345948129, "F" => 10.0 , "V" => 0.0 ),
"test47" => Dict("E" => 145.86499149789455, "F" => 10.0 , "V" => 0.0 ),
"test53" => Dict("E" => 137.36056394868902, "F" => 10.0 , "V" => 0.0 ),
"test54" => Dict("E" => 136.08276348795434, "F" => 10.0 , "V" => 0.0 ),
"test55" => Dict("E" => 134.83997249264843, "F" => 10.0 , "V" => 0.0 ),
"test56" => Dict("E" => 133.6306209562122, "F" => 10.0 , "V" => 0.0 ),
"test57" => Dict("E" => 132.4532357065044, "F" => 10.0 , "V" => 0.0 ),
"test58" => Dict("E" => 131.30643285972255, "F" => 10.0 , "V" => 0.0 ),
"test59" => Dict("E" => 130.18891098082386, "F" => 10.0 , "V" => 0.0 ),
"test88" => Dict("E" => 106.60035817780522, "F" => 10.0 , "V" => 0.0 ),
"test92" => Dict("E" => 104.2572070285374, "F" => 10.0 , "V" => 0.0 ),
"test94" => Dict("E" => 103.14212462587933, "F" => 10.0 , "V" => 0.0 ),
"test105" => Dict("E" => 97.59000729485332, "F" => 10.0 , "V" => 0.0 ),
"test107" => Dict("E" => 96.67364890456635, "F" => 10.0 , "V" => 0.0 ),
"test120" => Dict("E" => 91.28709291752769, "F" => 10.0 , "V" => 0.0 ),
"test121" => Dict("E" => 90.9090909090909, "F" => 10.0 , "V" => 0.0 ),
"test122" => Dict("E" => 90.53574604251853, "F" => 10.0 , "V" => 0.0 ),
"test127" => Dict("E" => 88.73565094161138, "F" => 10.0 , "V" => 0.0 ),
"test128" => Dict("E" => 88.38834764831843, "F" => 10.0 , "V" => 0.0 ),
"test133" => Dict("E" => 86.71099695241199, "F" => 10.0 , "V" => 0.0 ),
"test135" => Dict("E" => 86.06629658238704, "F" => 10.0 , "V" => 0.0 ),
"test137" => Dict("E" => 85.43576577167609, "F" => 10.0 , "V" => 0.0 ),
"test139" => Dict("E" => 84.8188929679971, "F" => 10.0 , "V" => 0.0 ),
"test141" => Dict("E" => 84.2151921066519, "F" => 10.0 , "V" => 0.0 ),
"test150" => Dict("E" => 81.64965809277261, "F" => 10.0 , "V" => 0.0 ),
"test152" => Dict("E" => 81.11071056538127, "F" => 10.0 , "V" => 0.0 ),
"test165" => Dict("E" => 77.8498944161523, "F" => 10.0 , "V" => 0.0 ),
"test167" => Dict("E" => 77.38232325341369, "F" => 10.0 , "V" => 0.0 ),
"test178" => Dict("E" => 74.95316889958615, "F" => 10.0 , "V" => 0.0 ),
"test180" => Dict("E" => 74.53559924999298, "F" => 10.0 , "V" => 0.0 ),
"test182" => Dict("E" => 74.12493166611011, "F" => 10.0 , "V" => 0.0 ),
"test195" => Dict("E" => 71.61148740394329, "F" => 10.0 , "V" => 0.0 ),
"test197" => Dict("E" => 71.24704998790965, "F" => 10.0 , "V" => 0.0 ),
"test210" => Dict("E" => 69.00655593423542, "F" => 10.0 , "V" => 0.0 ),
"test212" => Dict("E" => 68.68028197434451, "F" => 10.0 , "V" => 0.0 ),
"test225" => Dict("E" => 66.66666666666667, "F" => 10.0 , "V" => 0.0 ),
"test227" => Dict("E" => 66.3723311599972, "F" => 10.0 , "V" => 0.0 ),
"test255" => Dict("E" => 62.62242910851495, "F" => 10.0 , "V" => 0.0 ),
"test257" => Dict("E" => 62.37828615518053, "F" => 10.0 , "V" => 0.0 ),
"test272" => Dict("E" => 60.63390625908324, "F" => 10.0 , "V" => 0.0 ),
"test285" => Dict("E" => 59.23488777590923, "F" => 10.0 , "V" => 0.0 ),
"test287" => Dict("E" => 59.02813361009553, "F" => 10.0 , "V" => 0.0 ),
"test360" => Dict("E" => 52.70462766947299, "F" => 10.0 , "V" => 0.0 ),
"test362" => Dict("E" => 52.55883312276367, "F" => 10.0 , "V" => 0.0 ),
"test422" => Dict("E" => 48.679238351123544, "F" => 10.0 , "V" => 0.0 ),
"test437" => Dict("E" => 47.836487323493984, "F" => 10.0 , "V" => 0.0 ),
"test452" => Dict("E" => 47.03604341917986, "F" => 10.0 , "V" => 0.0 ))
model = acemodel(elements = [:W, :He],
                Eref = [:W => -1908.566156526, :He => -78.51555127829],
                order = 3,
                rcut = 5.0,
                totaldegree = 20)#,
                # pair_envelope = (:r, 3),
                # pair_transform = (:agnesi, 1, 4),
                # envelope = (:x, 2, 3),
                # transform = (:agnesi, 2, 4))
@show length(model.basis);
acefit!(model, dataset;
weights=weights);
ACEpotentials.linear_errors(dataset, model);
export2json("jace.json", model)
export2lammps("jace.yace", model)