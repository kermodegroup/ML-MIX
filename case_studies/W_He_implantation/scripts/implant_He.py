# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ---------------------------------------------------------------------

from ase.lattice.cubic import BodyCenteredCubic
from matscipy import parameter
import sys
sys.path.append(".")
import params
import ase.io
from ase.io.lammpsdata import write_lammps_data, read_lammps_data
from ase.optimize.precon import PreconLBFGS
import ase
import numpy as np
import os
script_path = os.path.dirname(os.path.abspath(__file__))
from lammps import lammps
import matplotlib.pyplot as plt


parallel = parameter("parallel", True)
try:
    import mpi4py
    rank = mpi4py.MPI.COMM_WORLD.Get_rank()
    nprocs = mpi4py.MPI.COMM_WORLD.Get_size()
except:
    if parallel:
        raise (ModuleNotFoundError, "mpi4py could not be imported!")
    else:
        print("No mpi4py found, proceeding in serial...")
        rank = 0
        nprocs = 1

el = parameter("el")
a0 = parameter("a0")
directions = parameter("directions")
nx = parameter("nx")
ny = parameter("ny")
nz = parameter("nz")
vacuum = parameter("vacuum")
calc = parameter("calc")
He_energy = parameter("He_energy")
lmps_cmds = parameter("lmps_cmds")
nrepeats = parameter("nrepeats")
T = parameter("T")
dump_files = parameter("dump_files")
dump_freq = parameter("dump_freq")
He_place_dist = parameter("He_place_dist")
second_stage_time_per_loop = parameter("second_stage_time_per_loop",0.02) #ps
second_stage_max_time = parameter("second_stage_max_time", 3) #ps
equilibration_time = parameter("equilibration_time", 10000)
equilibration_dump = parameter("equilibration_dump", False)
equilibration_dump_freq = parameter("equilibration_dump_freq", 100)
gpu = parameter("gpu", False)
thermo_freq = parameter("thermo_freq", 100)

if gpu:
    gpu_mode = parameter("gpu_mode")
multi_potential = parameter("multi_potential",False)
if not gpu:
    lmps = lammps()
else:
    args = [
    "-k", "on", "g",f"{nprocs}",  # Kokkos GPU usage
    "-sf", "kk",                 # Use Kokkos style
    "-pk", "kokkos", "newton", "on", "neigh", "half"  # Kokkos package settings
    ]
    lmps = lammps(cmdargs=args)
    if gpu_mode == "mliap":
        import lammps.mliap
        lammps.mliap.activate_mliappy_kokkos(lmps)


if rank == 0:
    supercell = BodyCenteredCubic(size=[1,1,1],symbol=el,latticeconstant=a0,pbc=(1,1,1),directions=directions)
    large_slab = supercell * (nx, ny, nz)

    cell = large_slab.get_cell()
    cell[2,2] += vacuum
    large_slab.set_cell(cell, scale_atoms=False)
    large_slab.set_pbc((1, 1, 0))
    #translate whole slab up by 1/3 of the vacuum
    large_slab.positions[:,2] += vacuum/3

    write_lammps_data(f"W_0K.data", large_slab, atom_style='atomic', velocities=True, specorder=['W', 'He'])

lmps.command("clear")
if multi_potential:
    if not gpu:
        lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/mlmlplugin.so')
        lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/langevinmlmlplugin.so')
        lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/hybridoverlaymlmlplugin.so')
lmps.command("units metal")
lmps.command("dimension 3")
lmps.command("boundary p p f")
lmps.command("atom_style atomic")
lmps.command("atom_modify map yes")
if multi_potential:
    if gpu:
        lmps.command('fix eval_pot all property/atom d_potential_1 d_potential_2 d_eval_1 d_eval_2 ghost yes')
    else:
        lmps.command('fix eval_pot all property/atom i2_potential 2 ghost yes')
        lmps.command('fix eval_arr all property/atom d2_eval 2 ghost yes')

lmps.command("read_data W_0K.data")
lmps.command("mass 1 183.84")
lmps.command("mass 2 4.0026")
lmps.commands_list(lmps_cmds)

lmps.command("fix 1 all nve")
lmps.command("timestep 0.001")
lmps.command(f"thermo {thermo_freq}")
lmps.command(f"fix langevin all langevin {T} {T} 0.01 12345678 zero yes")
if equilibration_dump:
    if not multi_potential:
        lmps.command(f"dump myDump all custom {equilibration_dump_freq} ./dump_eq.lammpstrj id mass type x y z vx vy vz fx fy fz")
    else:
        if not gpu:
            lmps.command(f"dump myDump all custom {equilibration_dump_freq} ./dump_eq.lammpstrj id mass type x y z vx vy vz fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]")
        else:
            lmps.command(f"dump myDump all custom {equilibration_dump_freq} ./dump_eq.lammpstrj id mass type x y z vx vy vz fx fy fz d_potential_1 d_potential_2 d_eval_1 d_eval_2")
lmps.command(f"run {equilibration_time}")

lmps.command("write_data W_He_in.data nofix")

if rank == 0:

    equilibrated_slab = read_lammps_data("W_He_in.data",atom_style="atomic")
    if not os.path.exists("./data"):
        os.makedirs("./data")

    #get the velocity of He atom based on energy
    He_v = np.sqrt(2*He_energy/4.0026)

    # convert He_v to A/ps
    He_v_A_ps = He_v*98.17614
    print(He_v_A_ps)
    #timestep is selected such that atom moves 0.02 A per timestep
    tstep = 0.02/He_v_A_ps
    if tstep > 0.0001:
        tstep = 0.0001
        dist_per_step = 0.0001*He_v_A_ps
        nsteps = int(3*He_place_dist/dist_per_step)
    else:
        #nsteps is selected such that the atom moves 1.5*He_place_dist
        nsteps = int(3*He_place_dist/0.02)

    print(f"IMPLANT TIMESTEP: {tstep}")
    print(f"IMPLANT NSTEP: {nsteps}")

    for n in range(nrepeats):
        # add He at a height of half vacuum above surface
        random_x = np.random.uniform(0, equilibrated_slab.get_cell()[0,0])
        random_y = np.random.uniform(0, equilibrated_slab.get_cell()[1,1])
        velocities_slab = equilibrated_slab.get_velocities()
        positions = equilibrated_slab.get_positions()
        top_W = np.max(positions[equilibrated_slab.symbols == "W"][:,2])
        He_atom = ase.Atom('He', position=(random_x, random_y, top_W + He_place_dist))
        equilibrated_slab_with_He = equilibrated_slab + He_atom
        
        #set the velocity of He atom
        velocities = equilibrated_slab_with_He.get_velocities()
        velocities[:-1] = velocities_slab
        velocities[-1] = [0, 0, -He_v]
        equilibrated_slab_with_He.set_velocities(velocities)
        ase.io.write(f"./data/W_He_in_{n}.xyz", equilibrated_slab_with_He, format='extxyz')
        #write the lammps data file
        write_lammps_data(f"./data/W_He_in_{n}.data", equilibrated_slab_with_He, atom_style='atomic', velocities=True, specorder=['W', 'He'])
        surface_pos = top_W
else:
    tstep = 0.0
    nsteps = 0
    surface_pos = 0.0

tstep = np.array([tstep], dtype=float)
nsteps = np.array([nsteps], dtype=int)
surface_pos = np.array([surface_pos], dtype=float)

if parallel:
    mpi4py.MPI.COMM_WORLD.Bcast(tstep, root=0)
    mpi4py.MPI.COMM_WORLD.Bcast(nsteps, root=0)
    mpi4py.MPI.COMM_WORLD.Bcast(surface_pos, root=0)
tstep = tstep[0]
nsteps = nsteps[0]
surface_pos = surface_pos[0]

reflected_states = np.zeros((nrepeats,), dtype=bool)
implantation_depths = np.zeros((nrepeats,), dtype=float)
final_He_vs = np.zeros((nrepeats,),dtype=float)
for i in range(nrepeats):
    lmps.command("clear")
    if multi_potential:
        if not gpu:
            lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/mlmlplugin.so')
            lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/langevinmlmlplugin.so')
            lmps.command(f'plugin load {script_path}/../../../LAMMPS_plugin/build/hybridoverlaymlmlplugin.so')
    
    lmps.command("units metal")
    lmps.command("dimension 3")
    lmps.command("boundary p p f")
    lmps.command("atom_style atomic")
    lmps.command("atom_modify map yes")
    if multi_potential:
        if gpu:
            lmps.command('fix eval_pot all property/atom d_potential_1 d_potential_2 d_eval_1 d_eval_2 ghost yes')
        else:
            lmps.command('fix eval_pot all property/atom i2_potential 2 ghost yes')
            lmps.command('fix eval_arr all property/atom d2_eval 2 ghost yes')
    lmps.command(f"read_data ./data/W_He_in_{i}.data")
    lmps.command("mass 1 183.84")
    lmps.command("mass 2 4.0026")
    lmps.commands_list(lmps_cmds)

    lmps.command("fix 1 all nve")
    lmps.command(f"timestep {tstep}")

    lmps.command(f"thermo {thermo_freq}")

    if dump_files:
        if not multi_potential:
            lmps.command(f"dump myDump all custom {dump_freq} ./data/dump_{i}.lammpstrj id mass type x y z vx vy vz fx fy fz")
        else:
            if not gpu:
                lmps.command(f"dump myDump all custom {dump_freq} ./data/dump_{i}.lammpstrj id mass type x y z vx vy vz fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]")
            else:
                print(f"setting up dump {i}")
                lmps.command(f"dump myDump all custom {dump_freq} ./data/dump_{i}.lammpstrj id mass type x y z vx vy vz fx fy fz d_potential_1 d_potential_2 d_eval_1 d_eval_2")
    
    lmps.command("group group_he type 2")
    lmps.command("compute vels group_he property/atom vx vy vz")
    lmps.command("compute vel1 group_he reduce ave c_vels[1] c_vels[2] c_vels[3]")
    lmps.command("variable vmag equal sqrt(c_vel1[1]^2 + c_vel1[2]^2 + c_vel1[3]^2)")
    lmps.command("thermo_style custom step temp v_vmag")

    lmps.command(f"dump He_v_dump group_he custom {50} ./data/HeV_dump.lammpstrj id mass type x y z vx vy vz")
    lmps.command(f"run {nsteps}")
    lmps.command("write_data ./data/W_He_out.data nofix")
    loop_count = 0
    v_comp = np.sqrt((3*(1.381e-23)*T)/(4*(1.66e-27)))/100 # ideal gas RMS of He at 1000K in A/ps
    max_loops = int(second_stage_max_time/second_stage_time_per_loop)
    He_vs = []
    while True:
        loop_count += 1
        break_loop = False
        thermalised = False
        max_time = False
        if rank == 0:
            #read the lammps data file
            final_state = read_lammps_data(f"./data/W_He_out.data",atom_style="atomic")
            ase.io.write(f"./data/W_He_out_first_stage_{i}.xyz", final_state, format='extxyz')
            velocities = final_state.get_velocities()
            positions = final_state.get_positions()

            #get the position of the top most W atom
            top_W = np.max(positions[final_state.symbols == "W"][:,2])
            He_pos = positions[-1,2]
            He_v = velocities[-1,2]

            #set new tstep such that He moves 0.02 A per tstep
            He_traj = ase.io.read("./data/HeV_dump.lammpstrj",index=":")
            v_vals = []
            for t in He_traj:
                v_vals.append(np.linalg.norm(t.get_velocities())) #ASE units
            
            mean_v = np.mean(v_vals)
            rms_v = np.sqrt(np.mean(np.array(v_vals)**2))
            rms_v_aps = rms_v*98.17614
            print(f"root mean square He V is {rms_v_aps} A/ps")
            print(f"Comparison He velocity at 1000K is {v_comp}")
            He_vs.append(rms_v_aps)
            #if rms_v_aps<1.5*v_comp:
            #thermalised=True
            new_tstep = 0.02/(mean_v*98.17614) #ps (LAMMPS units)
            if new_tstep > 0.0001:
                new_tstep = 0.0001 # at most 0.1 fs
            print(f"New timestep is {new_tstep}")
            bottom_W = np.min(positions[final_state.symbols == "W"][:,2])

            if ((He_pos - surface_pos) > 2) and (He_v > 0):
                print("He atom has reflected!")
                reflected_states[i] = True
                implantation_depths[i] = 0.0
                break_loop = True
            else:
                reflected_states[i] = False
            
            if ((He_pos - bottom_W) < -2) and (He_v < 0):
                print("He atom has passed through!")
                reflected_states[i] = False
                implantation_depths[i] = -10000.0
                break_loop = True

        
        break_loop = np.array([break_loop], dtype=bool)
        #thermalised = np.array([thermalised], dtype=bool)
        if parallel:
            mpi4py.MPI.COMM_WORLD.Bcast(break_loop, root=0)
        #    mpi4py.MPI.COMM_WORLD.Bcast(thermalised, root=0)
        break_loop = break_loop[0]
        if break_loop:
            break
        #if thermalised:
        #    break

        lmps.command(f"timestep {new_tstep}")
        lmps.command(f"undump He_v_dump")
        lmps.command(f"dump He_v_dump group_he custom {10} ./data/HeV_dump.lammpstrj id mass type x y z vx vy vz")
        lmps.command(f"run {int(second_stage_time_per_loop/new_tstep)}")
        lmps.command(f"write_data ./data/W_He_out.data nofix")

        if loop_count>=max_loops:
            print("Max time reached!")
            max_time = True
            break
    

    if rank == 0:
        if not break_loop:
            #read the lammps data file
            final_state = read_lammps_data(f"./data/W_He_out.data",atom_style="atomic")
            ase.io.write(f"./data/W_He_out_second_stage_{i}.xyz", final_state, format='extxyz')
            velocities = final_state.get_velocities()
            positions = final_state.get_positions()

            #get the position of the top most W atoms
            top_Ws = np.argsort(positions[final_state.symbols == "W"][:,2])
            # now cycle through the top 10 W atoms - if the difference between the top most
            # and the second top is more than 5 A then print that the top W has been ejected
            if (positions[final_state.symbols == "W"][top_Ws[-1],2] - positions[final_state.symbols == "W"][top_Ws[-2],2]) > 5:
                print("WARNING: W ejected!")
            He_pos = positions[-1,2]
            He_v = velocities[-1,2]

            implantation_depths[i] = He_pos - surface_pos
            final_He_vs[i] = He_vs[-1]
            print(f"Implantation depth for {i} is {implantation_depths[i]}")
            np.savetxt(f"./data/He_vs_{i}.txt", np.array(He_vs))
            plt.figure()
            plt.plot(He_vs)
            plt.savefig(f"./data/He_vs_{i}.png")
    
        # write numpy file for arrays on the go
        big_arr = np.zeros((nrepeats, 3))
        big_arr[:,0] = reflected_states
        big_arr[:,1] = implantation_depths
        big_arr[:,2] = final_He_vs
        np.savetxt(f"./data/implantation_depths.txt", big_arr, fmt="%s")


