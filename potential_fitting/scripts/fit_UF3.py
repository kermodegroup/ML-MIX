import numpy as np
import os

from uf3.data import composition
from uf3.representation import bspline
from uf3.data import io
from uf3.representation import process
from uf3.regression import least_squares


import ase.io

import jax.numpy as jnp
from jax import vmap, jit

import argparse

import matplotlib.pyplot as plt
#get arguments from command line in the form of
#python fit_matching_UF3.py path model_name_1 model_name_2 alpha
def parse_args():
    parser = argparse.ArgumentParser(description='Fit UF3 model with constraints.')
    parser.add_argument('--raw_output', help="Path to output file for raw JSON potential", type=str, required=True)
    parser.add_argument('--lammps_output', help="Path to output file for LAMMPS potential", type=str, required=True)
    parser.add_argument('--analysis_path', help="Path to output file for analysis", type=str, required=True)
    parser.add_argument('--model_name_unconstrained', help="Output name of unconstrained model", type=str, required=True)
    parser.add_argument('--model_name_constrained', help="Output name of constrained model", type=str, required=True)
    parser.add_argument('--alpha', help="Hard constraint strength", type=float, required=True)
    parser.add_argument('--element', help="Input element", type=str, required=True)
    parser.add_argument('--rcut2b', help="Cutoff radius for the 2 body term", type=float, default=6.0)
    parser.add_argument('--rcut3b', help="Cutoff radius for the 3 body term", type=float, default=3.5)
    parser.add_argument('--resolution2b', help="Resolution of 2 body UF3 model", type=int, default=15)
    parser.add_argument('--resolution3b', help="Resolution of 3 body UF3 model", type=int, default=6)
    parser.add_argument('--rcutmin2b', help="Minimum cutoff radius for the 2 body term", type=float, default=0.001)
    parser.add_argument('--rcutmin3b', help="Minimum cutoff radius for the 3 body term", type=float, default=1.5)
    parser.add_argument('--energy_key', help="Key for energy in the data", type=str, default="energy")
    parser.add_argument('--skip_featurization', help="Skip featurization step", action='store_true')
    parser.add_argument('--lambda_max', help="Maximum lambda value for inequality constrained regression", type=float, default=10**6)
    parser.add_argument('--ridge1b', help="Ridge regularization for 1 body term", type=float, default=0.0)
    parser.add_argument('--ridge2b', help="Ridge regularization for 2 body term", type=float, default=0.0)
    parser.add_argument('--ridge3b', help="Ridge regularization for 3 body term", type=float, default=1e-8)
    parser.add_argument('--curvature2b', help="Curvature regularization for 2 body term", type=float, default=1e-8)
    parser.add_argument('--curvature3b', help="Curvature regularization for 3 body term", type=float, default=0.0)
    parser.add_argument('--lambda_points', help="Number of lambda points for grid search", type=int, default=1000)
    return parser.parse_args()

args = parse_args()
path = args.raw_output
model_name_1 = args.model_name_unconstrained
model_name_2 = args.model_name_constrained
alpha = args.alpha
skip_featurization = args.skip_featurization
rcut2b = args.rcut2b
rcut3b = args.rcut3b
rcutmin2b = args.rcutmin2b
rcutmin3b = args.rcutmin3b
el = args.element
lambda_max = args.lambda_max
resolution2b = args.resolution2b
resolution3b = args.resolution3b
ridge1b = args.ridge1b
ridge2b = args.ridge2b
ridge3b = args.ridge3b
curvature2b = args.curvature2b
curvature3b = args.curvature3b
lambda_points = args.lambda_points
analysis_path = args.analysis_path


element_list = [el.upper()]
degree = 3

chemical_system = composition.ChemicalSystem(element_list=element_list,
                                             degree=degree)
print(chemical_system)

print("Trios:", chemical_system.interactions_map[3])

#set hyperparameters

r_min_map = {("W", "W"): rcutmin2b,
             ("W", "W", "W"): [rcutmin3b, rcutmin3b, rcutmin3b],
            }
r_max_map = {("W", "W"): rcut2b,
             ("W", "W", "W"): [rcut3b, rcut3b, 2*rcut3b],
            }
resolution_map = {("W", "W"): resolution2b,
                  ("W", "W", "W"): [resolution3b, resolution3b, 2*resolution3b],
                 }
trailing_trim = 3
leading_trim = 0
bspline_config = bspline.BSplineBasis(chemical_system,
                                    r_min_map=r_min_map,
                                    r_max_map=r_max_map,
                                    resolution_map=resolution_map,
                                    leading_trim=leading_trim,
                                    trailing_trim=trailing_trim)
                                    #offset_1b=False)


#load in eq constraint configs
eq_constraint_filename = "Ah_data/data_for_Ah.xyz"
#traj_c = ase.io.read(eq_constraint_filename, index=":")
# for t in traj_c:
#     t.calc.results['energy'] -= len(t)*E0_energy

# modified_fname = f"{path}/modified_hard_constraint"
# ase.io.write(modified_fname,traj_c,parallel=False)
data_coordinator = io.DataCoordinator()
data_coordinator.dataframe_from_trajectory(eq_constraint_filename,
                                        prefix='EOS')
df_constraint = data_coordinator.consolidate()

bulk_reference_filename = "Ah_data/just_bulk.xyz"
data_coordinator = io.DataCoordinator()
data_coordinator.dataframe_from_trajectory(bulk_reference_filename,
                                        prefix='bulk_ref')
df_bulk_ref = data_coordinator.consolidate()

data_filename = "As_data/data_for_As.xyz"
traj_d = ase.io.read(data_filename, index=":")
print('TRAJ LEN',len(traj_d))
data_coordinator = io.DataCoordinator()

data_coordinator.dataframe_from_trajectory(data_filename,
                                        prefix='ace')
df_data = data_coordinator.consolidate()


if not skip_featurization:
    #delete the file 'test_uff.h5' if it exists


    #do both models
    if os.path.exists(f'{path}/fit_data.h5'):
        os.remove(f'{path}/fit_data.h5')
    if os.path.exists(f'{path}/hard_constraint.h5'):
        os.remove(f'{path}/hard_constraint.h5')
    if os.path.exists(f'{path}/bulk_reference.h5'):
        os.remove(f'{path}/bulk_reference.h5')
    representation = process.BasisFeaturizer(bspline_config)
    representation.batched_to_hdf_uff(f'{path}/hard_constraint.h5',df_data=df_constraint,n_jobs=1)
    representation.batched_to_hdf_uff(f'{path}/fit_data.h5',df_data=df_data,n_jobs=1)
    representation.batched_to_hdf_uff(f'{path}/bulk_reference.h5',df_data=df_bulk_ref,n_jobs=1)



# regularizer = bspline_config.get_regularization_matrix(ridge_1b=0.0,
#                                                     ridge_2b=1e-2,
#                                                     ridge_3b=1e-3,
#                                                     curvature_2b=1e-2,
#                                                     curvature_3b=1e-3)
regularizer = bspline_config.get_regularization_matrix(ridge_1b=ridge1b,
                                                       ridge_2b=ridge2b,
                                                       ridge_3b=ridge3b,
                                                       curvature_2b=curvature2b,
                                                       curvature_3b=curvature3b)




model_uff = least_squares.WeightedLinearModel(bspline_config,
                                            regularizer=regularizer)



print('fitting unconstrained model reference...')
model_uff.fit_from_file(f'{path}/fit_data.h5', list(df_data.index), 
                        batch_size=2500, UFF=True,
                        energy_key="energy")


c_unconstrained = model_uff.coefficients.copy()
print('zero_coeff', model_uff.coefficients[0])
#model_uff.coefficients[0] += E0_energy
model_uff.to_json(f"{path}/{model_name_1}.json")

#now build constraint matrices

print('building design matrix for constraint 1')
Ah,Yh = model_uff.build_design_matrix(f'{path}/hard_constraint.h5', list(df_constraint.index), 
                                            UFF=True,
                                            energy_key="energy",
                                            just_energies=True)


bulk_ref, Y_bulk = model_uff.build_design_matrix(f'{path}/bulk_reference.h5', list(df_bulk_ref.index),
                                                UFF=True,
                                                energy_key="energy",
                                                just_energies=True)

#subtract bulk_ref from all rows of Ah

for i in range(np.shape(Ah)[0]):
    Ah[i,:] -= bulk_ref.flatten()

# here, print the rank of the design matrix
print('rank of Ah:',np.linalg.matrix_rank(Ah))

# print the number of free parameters
print('number of free parameters:',np.shape(Ah)[1])

Yh = Yh - Y_bulk


print('building design matrix for data 1')

As,Ys = model_uff.build_design_matrix(f'{path}/fit_data.h5', list(df_data.index),
                                                UFF=True,
                                                energy_key="energy")


Omega = least_squares.freeze_regularizer(model_uff.regularizer, model_uff.mask)

if alpha == 0:
    print('alpha=0 is equality constrained regression, solving exactly...')

    #num basis funcs
    n_B = np.shape(Ah)[1]

    print('performing SVD for constraint...')
    n_C = np.shape(Ah)[0]
    U, S, VT = np.linalg.svd(Ah.T,full_matrices=True)
    print('Done!')
    S = np.diag(S)
    #print how many non zero singular values there are
    Yh_trans = VT@Yh
    c1 = np.linalg.solve(S.T,Yh_trans)

    #now get A1 and A2
    A = As @ U

    A1 = A[:,:n_C]
    A2 = A[:,n_C:]

    #A2 is our new design matrix
    #y-A1ytilde is our new target

    b = Ys - (A1 @ c1)

    #fit the model
    print('performing solve with regularisation...')
    # res = np.linalg.lstsq(A2,b,rcond=None)
    #build gram matrix with A2
    G = A2.T @ A2
    #get ordinate
    d = A2.T @ b

    regularizer = least_squares.freeze_regularizer(model_uff.regularizer, model_uff.mask)
    regularizer = np.dot(regularizer.T, regularizer)

    #compute M = U.T @ regularizer @ U
    M = U.T @ regularizer @ U

    #get MBB
    M_BB = M[n_C:,n_C:]
    M_BC = M[n_C:,:n_C]

    d = d - (M_BC @ c1)

    c2 = least_squares.lu_factorization(G + M_BB, d)

    #finally build c back from Q@[y,z]

    c = U @ np.hstack([c1,c2])

else:
    print('alpha>0 is inequality constrained regression, solving with lagrange multiplier method...')
    #objective function to minimize with lagrange multiplier lambda
    def h(lambda_):

        A = As.T@As + lambda_*Ah.T@Ah + Omega.T@Omega
        b = As.T@Ys + lambda_*Ah.T@Yh
        #print(f'Solving for lambda={lambda_}')
        c = jnp.linalg.solve(A,b)
        #print('Solved')
        val = jnp.sum((Ah@c - Yh)**2) - alpha
        #print('val:',val)
        return val, c

    def objective(lambda_):
        return h(lambda_)[0]

    val, c = h(0.0)
    print('h(0):', val)
    if val < 0:
        print("Constraint is already satisfied")
    else:
        print('Constraint not satisfied')
        print("Performing grid search for lagrange multiplier...")
        #vmap the objective function
        objective_v = jit(vmap(objective))

        #now plug in a whole set of lambda values
        lambda_values = jnp.linspace(0.0, lambda_max, lambda_points)
        h_vals = objective_v(lambda_values)

        #plot and save the h_vals
        plt.plot(lambda_values, (h_vals))
        plt.xlabel('lambda')
        plt.ylabel('h(lambda)')
        plt.yscale('log')
        plt.savefig(f'{analysis_path}/h_vals.png')


        objective = jit(objective)
        min_idx = jnp.argmin(jnp.abs(h_vals))
        min_lambda = lambda_values[min_idx]
        
        # first, check h for lambda = 0
        print('min_lambda:',min_lambda)
        val, c = h(min_lambda)
        print('min_h:',val)
        #write min lambda to file
        np.savetxt(f'{analysis_path}/min_lambda.txt',np.array([min_lambda]))

c = least_squares.revert_frozen_coefficients(c,
                                        model_uff.n_feats,
                                        model_uff.mask,
                                        model_uff.frozen_c,
                                        model_uff.col_idx)


Ah_unfrozen = least_squares.unfreeze_design_matrix(Ah,
                                        model_uff.n_feats,
                                        model_uff.mask,
                                        model_uff.frozen_c,
                                        model_uff.col_idx)

print("unfrozen num params:",np.shape(Ah_unfrozen)[1])


print('c[0]',c[0])
model_uff.coefficients = c.copy()
#model_uff.coefficients[0] += E0_energy
model_uff.to_json(f"{path}/{model_name_2}.json")


print('testing hard constraint on unconstrained fit...')
# print(list(Ah_unfrozen@c_unconstrained - Yh))

print('testing hard constraint on constrained fit...')
# print(list(Ah_unfrozen@c - Yh))

print("---------------------------")
print('testing unconstrained rmse...')
print('unconstrained fit')
model_uff.coefficients = c_unconstrained.copy()
y_e, p_e, y_f, p_f, rmse_e_unconstrained, rmse_f_unconstrained = model_uff.batched_predict(f'{path}/fit_data.h5', 
                                                           keys=list(df_data.index),
                                                           UFF=True)
print('rmse energy:',rmse_e_unconstrained)
print('rmse force:',rmse_f_unconstrained)

print('---------------------------')
print('testing constrained rmse...')
print('model 1 constrained')

model_uff.coefficients = c.copy()
y_e, p_e, y_f, p_f, rmse_e_constrained, rmse_f_constrained = model_uff.batched_predict(f'{path}/fit_data.h5', 
                                                           keys=list(df_data.index),
                                                           UFF=True)
print('rmse energy constrained:',rmse_e_constrained)
print('rmse force constrained:',rmse_f_constrained)

# save to file
with open(f'{analysis_path}/errors.txt','w') as f:
    f.write(f"alpha: {alpha}\n")
    f.write(f"lambda_root: {min_lambda}\n")
    f.write(f"hval at root: {val}\n")
    f.write(f"unconstrained_energy_err: {rmse_e_unconstrained}\n")
    f.write(f"unconstrained_force_err: {rmse_f_unconstrained}\n")
    f.write(f"constrained_energy_err: {rmse_e_constrained}\n")
    f.write(f"constrained_force_err: {rmse_f_constrained}\n")








