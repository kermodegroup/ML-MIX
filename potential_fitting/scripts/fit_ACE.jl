using Pkg
Pkg.activate("../../../.") # activate the environment
using Distributed
using DelimitedFiles
using Plots
import ACE1x:_set_params!
using LinearAlgebra: Diagonal
using ArgParse
function parse_arguments()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--raw_output"
            help = "Path to output file for raw JSON potential"
            arg_type = String
            required = true
        "--lammps_output"
            help = "Path to output file for LAMMPS potential"
            arg_type = String
            required = true
        "--analysis_path"
            help = "Path to output file for analysis"
            arg_type = String
            required = true
        "--model_name_unconstrained"
            help = "Output name of unconstrained model"
            arg_type = String
            required = true
        "--model_name_constrained"
            help = "Output name of constrained model"
            arg_type = String
            required = true
        "--alpha"
            help = "Hard constraint strength"
            arg_type = Float64
            required = true
        "--nprocs"
            help = "Number of processors to use"
            arg_type = Int64
            required = true
        "--elements"
            help = "Input element"
            arg_type = String
            required = true
        "--e_ref"
            help = "Reference energies for each element"
            arg_type = Float64
            required = true
        "--order"
            help = "Correlation order of the model"
            arg_type = Int64
            required = true
        "--totaldegree"
            help = "Total degree of the model"
            arg_type = Int64
            required = true
        "--rcut"
            help = "Cutoff radius for the model"
            arg_type = Float64
            required = true
        "--energy_key"
            help = "Key for energy in the data"
            arg_type = String
            default = "energy"
        "--force_key"
            help = "Key for forces in the data"
            arg_type = String
            default = "forces"
        "--virial_key"
            help = "Key for virial stresses in the data"
            arg_type = String
            default = nothing
    end

    return parse_args(s)
end


# Parse arguments
args = parse_arguments()
nprocs = args["nprocs"]

addprocs(nprocs-1, exeflags="--project=$(Base.active_project())")
@everywhere using ACEpotentials
@everywhere using ACEfit
@everywhere import ACEfit: assemble, solve
@everywhere using JuLIP.MLIPs


jsonpath = args["raw_output"]
lammpspath = args["lammps_output"]
analysis_path = args["analysis_path"]
model_name_1 = args["model_name_unconstrained"]
model_name_2 = args["model_name_constrained"]
alpha = args["alpha"]
elements = Symbol.(split(args["elements"], ","))
e_ref = Dict(elements[i] => args["e_ref"][i] for i in 1:length(elements))
order = args["order"]
totaldegree = args["totaldegree"]
rcut = args["rcut"]
energy_key = args["energy_key"]
force_key = args["force_key"]
virial_key = args["virial_key"]


#set up model hyperparameters and basis

println("reading data...")
# Read data
As_data = read_extxyz("As_data/data_for_As.xyz")
Ah_data = read_extxyz("Ah_data/data_for_Ah.xyz")
#get total number of atoms in Ah_data
total_atoms = 0
for i in 1:length(Ah_data)
    global total_atoms
    total_atoms += length(Ah_data[i])
end
println("total atoms in hard constraint: ", total_atoms)
base_bulk_struct = read_extxyz("Ah_data/just_bulk.xyz")


println("setting up basis...")

model = acemodel(
    elements = elements,
    Eref = e_ref,
    order = order,
    totaldegree = totaldegree,
    rcut = rcut,
)
@show length(model.basis)



# fit reference unconstrained ACE model to As data
println("Fitting unconstrained reference model...")
solver = ACEfit.BLR(factorization=:svd)
P = smoothness_prior(model; p = 4) #this is our Tikhonov regularizer

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0 , "V" => 1.0 ))

#TODO write this out explicitly to stop assembling As twice.
acefit!(model, 
        As_data; 
        solver=solver, 
        prior = P, 
        weights=weights,
        energy_key = energy_key,
        force_key = force_key,
        virial_key = virial_key)

export2json("$jsonpath/$model_name_1.json", model)
export2lammps("$lammpspath/$model_name_1.yace", model)

@info("Training Error Table")
rmse, mae = ACEpotentials.linear_errors(As_data, model, energy_key=energy_key, force_key=force_key)
println(rmse)
println(rmse.second)
println(rmse.second["set"])
println(rmse.second["set"]["E"])
unconstrained_energy_err = rmse.second["set"]["E"]
unconstrained_force_err = rmse.second["set"]["F"]

# Now build Ah and As design matrices
Ah = map( Ah_data ) do data_point
    AtomsData(
        data_point;
        energy_key = energy_key,
        force_key  = nothing,
        virial_key = virial_key,
        weights    = weights,
        v_ref      = model.Vref,
    )
 end

As = map( As_data ) do data_point
    AtomsData(
        data_point;
        energy_key = energy_key,
        force_key  = "forces",
        virial_key = virial_key,
        weights    = weights,
        v_ref      = model.Vref,
    )
 end

Ah_base = map( base_bulk_struct ) do data_point
    AtomsData(
        data_point;
        energy_key = energy_key,
        force_key  = nothing,
        virial_key = nothing,
        weights    = weights,
        v_ref      = model.Vref,
    )
 end

Ah, Yh, Wh = ACEfit.assemble(Ah, 
                            model.basis)
                            #energy_ref= model.Vref)
As, Ys, Ws = ACEfit.assemble(As, 
                            model.basis)
                            #energy_ref= model.Vref)
                            
Ah_base, Yh_base, Wh_base = ACEfit.assemble(Ah_base, 
                            model.basis)
                            #energy_ref= model.Vref)
println(size(Ah))
println(size(Ah_base))

#subtract Ah_base from each row of Ah
println(Ah_base[1,1], ' ', Ah_base[1,2])
println(Ah[1,1], ' ',Ah[1,2])
Ah = Ah .- Ah_base
println(Ah[1,1], ' ',Ah[1,2])

Yh = Yh .- Yh_base


function h(lambda_::Float64)
    Ap = vcat(Diagonal(Ws)*As, sqrt(lambda_)*Diagonal(Wh)*Ah)/P
    Y = vcat(Ws.*Ys, sqrt(lambda_)*Wh.*Yh)
    res = ACEfit.solve(solver, Ap, Y)
    c = P \ res["C"]
    #alpha is in RMS energy per atom
    val = sqrt((1/total_atoms)*sum((Ah*c - Yh).^2)) - alpha
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
if (val<alpha)
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
lambda_root = interval_bisection(objective, trial_lambda/10.0, trial_lambda, alpha/10, lambda_vals, h_vals)
# Apply h to each element of lambda_vec using broadcasting
#vals = objective.(lambda_vec)

val, c = h(lambda_root)
push!(lambda_vals, lambda_root)
push!(h_vals, val)

_set_params!(model, c)
rmse, mae = ACEpotentials.linear_errors(As_data, 
                            model,
                            energy_key=energy_key,
                            force_key=force_key)

constrained_energy_err = rmse.second["set"]["E"]
constrained_force_err = rmse.second["set"]["F"]


export2json("$jsonpath/$model_name_2.json", model)
export2lammps("$lammpspath/$model_name_2.yace", model)

println("all lambda vals tried", lambda_vals)
println("corresponding h vals", h_vals)
# plot vals vs lambda_vec
scatter(lambda_vals, h_vals.+alpha)
scatter!([lambda_root], [val+alpha], label="selected root (within tol)", color="red")
#plot horizontal line at alpha
hline!([alpha], label="alpha")
xlabel!("lambda")
ylabel!("h(lambda) + alpha")
title!("lambda search")
#set axes to logarithmic
xaxis!(:log)
yaxis!(:log)
#save plot to file
savefig("$analysis_path/lambda_search.png")

#write errors to file
open("$analysis_path/errors.txt", "w") do f
    write(f, "alpha: $alpha\n")
    write(f, "lambda_root: $lambda_root\n")
    write(f, "hval at root: $val\n")
    write(f, "unconstrained_energy_err: $unconstrained_energy_err\n")
    write(f, "unconstrained_force_err: $unconstrained_force_err\n")
    write(f, "constrained_energy_err: $constrained_energy_err\n")
    write(f, "constrained_force_err: $constrained_force_err\n")
end
