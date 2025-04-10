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
using LinearAlgebra: Diagonal
using Plots
jsonpath = "./"
lammpspath = "./"

@info("Loading data...")
data = read_extxyz("otf-fit-data.xyz")
if length(data) == 0
    println("No data found in file")
    exit(1)
end

# now need to read in the extra As and Ah data
key_types = ["energy-key", "force-key", "virial-key", "pae-key", "mask-key"]

if "extra_hard_constraints" in keys(config)
    extra_hard_constraints = config["extra_hard_constraints"]
    for key in key_types
        if key in keys(extra_hard_constraints)
            extra_hard_constraints[key] = extra_hard_constraints[key]
        else
            extra_hard_constraints[key] = nothing
        end
    end
    Ah_data = read_extxyz(extra_hard_constraints["path"])
    alpha_val = extra_hard_constraints["alpha"]

    total_atoms = 0
    for i in 1:length(Ah_data)
        global total_atoms
        total_atoms += length(Ah_data[i])
    end
else
    extra_hard_constraints = nothing
end

if "extra_soft_constraints" in keys(config)
    extra_soft_constraints = config["extra_soft_constraints"]
    for key in key_types
        if key in keys(extra_soft_constraints)
            extra_soft_constraints[key] = extra_soft_constraints[key]
        else
            extra_soft_constraints[key] = nothing
        end
    end
    As_data = read_extxyz(extra_soft_constraints["path"])
else
    extra_soft_constraints = nothing
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

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0 , "V" => 1.0, "PAE" => 1.0))

# fit reference unconstrained ACE model to As data
solver = ACEfit.BLR(factorization=:svd)
P = smoothness_prior(model; p = smoothness_prior_strength) #this is our Tikhonov regularizer


if extra_hard_constraints != nothing
    Ah = map( Ah_data ) do data_point
        AtomsData(
            data_point;
            energy_key = extra_hard_constraints["energy-key"],
            force_key  = extra_hard_constraints["force-key"],
            virial_key = extra_hard_constraints["virial-key"],
            pae_key    = extra_hard_constraints["pae-key"],
            mask_key   = extra_hard_constraints["mask-key"],
            weights    = weights,
            v_ref      = model.Vref,
        )
    end
    Ah, Yh, Wh = ACEfit.assemble(Ah, model.basis)
end



if extra_soft_constraints != nothing
    As = map( As_data ) do data_point
        AtomsData(
            data_point;
            energy_key = extra_hard_constraints["energy-key"],
            force_key  = extra_soft_constraints["force-key"],
            virial_key = extra_soft_constraints["virial-key"],
            pae_key    = extra_soft_constraints["pae-key"],
            mask_key   = extra_soft_constraints["mask-key"],
            weights    = weights,
            v_ref      = model.Vref,
        )
    end
    As, Ys, Ws = ACEfit.assemble(As, model.basis)
end


As_otf = map( data ) do data_point
    AtomsData(
        data_point;
        energy_key = energy_key,
        force_key  = force_key,
        virial_key = virial_key,
        pae_key    = pae_key,
        mask_key   = mask_key,
        weights    = weights,
        v_ref      = model.Vref,
    )
 end
As_otf, Y_otf, W_otf = ACEfit.assemble(As_otf, 
                                        model.basis)

if extra_soft_constraints != nothing
    @info("Adding extra soft constraints...")
    As = vcat(As_otf, As)
    Ys = vcat(Y_otf, Ys)
    Ws = vcat(W_otf, Ws)
else
    @info("No extra soft constraints found...")
    As = As_otf
    Ys = Y_otf
    Ws = W_otf
end


if extra_hard_constraints == nothing
    @info("No hard constraints found, solving for unconstrained model...")
    Ap = (Diagonal(Ws)*As)/P
    Y = Ws.*Ys
    res = ACEfit.solve(solver, Ap, Y)
    c = P \ res["C"]
else
    @info("Fitting with hard constraints...")

    function h(lambda_::Float64)
        Ap = vcat(Diagonal(Ws)*As, sqrt(lambda_)*Diagonal(Wh)*Ah)/P
        Y = vcat(Ws.*Ys, sqrt(lambda_)*Wh.*Yh)
        res = ACEfit.solve(solver, Ap, Y)
        c = P \ res["C"]
        #alpha is in RMS energy per atom
        val = sqrt((1/total_atoms)*sum((Ah*c - Yh).^2)) - alpha_val
        return val, c
    end

    function objective(lambda_)
        val, c = h(lambda_)
        return val
    end

    function interval_bisection(f, a, b, tol, lambda_vals, h_vals, maxnum=10)
        it = 0
        while (((b-a)/2 > tol) && (it <= maxnum))
            c = (a+b)/2
            tmp_c = f(c)
            push!(lambda_vals, c)
            push!(h_vals, tmp_c)
            tmp_a = f(a)
            push!(lambda_vals, a)
            push!(h_vals, tmp_a)
            if abs(tmp_c) < tol
                return c
            elseif tmp_a*tmp_c < 0
                b = c
            else
                a = c
            end
            it += 1
        end
        if (it >= maxnum)
            println("Max number of iterations reached!")
        end
        return (a+b)/2
    end

    #first get h(0.0)
    val, c = h(0.0)
    println("value is, ", val)
    if (val<alpha_val)
        print("value is less than alpha!")
    end

    #Define a vector of lambda_ values
    #first take logarithmic steps through objective until we find a negative value

    lambda_vals = Float64[]
    h_vals = Float64[]

    trial_lambda = 1.0
    for i in 1:10
        global trial_lambda, lambda_vals, h_vals
        println("trial_lambda is, ", trial_lambda)
        push!(lambda_vals, trial_lambda)
        value = objective(trial_lambda)
        push!(h_vals, value)
        println("value is, ", value)
        if (value<0)
            print("value is less than 0")
            break
        end
        trial_lambda = trial_lambda*10.0
    end

    #now do interval bisection with lambda_ and lambda_/10.0
    lambda_root = interval_bisection(objective, trial_lambda/10.0, trial_lambda, alpha_val/10, lambda_vals, h_vals)
    # Apply h to each element of lambda_vec using broadcasting
    #vals = objective.(lambda_vec)

    val, c = h(lambda_root)
    push!(lambda_vals, lambda_root)
    push!(h_vals, val)

    println("all lambda vals tried", lambda_vals)
    println("corresponding h vals", h_vals)
    # plot vals vs lambda_vec
    scatter(lambda_vals, h_vals.+alpha_val)
    scatter!([lambda_root], [val+alpha_val], label="selected root (within tol)", color="red")
    #plot horizontal line at alpha
    hline!([alpha_val], label="alpha")
    xlabel!("lambda")
    ylabel!("h(lambda) + alpha")
    title!("lambda search")
    #set axes to logarithmic
    xaxis!(:log)
    yaxis!(:log)
    #save plot to file
    savefig("lambda_search.png")
end

_set_params!(model, c)


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
        write(f, "PAE RMSE: $pae_RMSE\n")
    end
end

println("########Exiting Julia########")

