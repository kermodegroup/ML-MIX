from ase.calculators.lammpslib import LAMMPSlib
import ase.io
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from ase.calculators.singlepoint import SinglePointCalculator

pot = '0ce218cb3e'

cmds = ['pair_style hybrid/overlay pace table spline 5401',
                        f'pair_coeff * * pace {pot}.yace Si', 
                        f'pair_coeff 1 1 table {pot}_pairpot.table Si_Si']

calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)


paths = ["all_expensive","mixed_9"]
image_num = [9, 9]
for j,path in enumerate(paths):
    files = [f"{path}/{i}.lammpstrj" for i in range(image_num[j])]

    images = []
    for i, file in enumerate(files):
        final_image = ase.io.read(file, index=-1)
        final_image.calc = calc
        # now form cell vectors ax, ay, az
        if i == 0:
            ix_iy_iz = np.zeros((len(final_image), 3))
            ref_positions = final_image.positions.copy()
            pos = final_image.positions.copy()
        else:
            cell = final_image.get_cell()
            # now go through each atom
            # and check if it's moved to a different image
            ix_iy_iz = np.zeros((len(final_image), 3))
            pos = final_image.positions.copy()
            ax = cell[0]
            ay = cell[1]
            az = cell[2]
            for j in range(np.shape(pos)[0]):
                # calculate the vector between the reference and current position
                r = pos[j,:] - ref_positions[j,:]
                for k in range(len(r)):
                    if abs(r[k]) > np.linalg.norm(cell[k]/2):
                        ix_iy_iz[j,k] = np.sign(r[k])

                pos[j,:] -= ix_iy_iz[j,0]*ax + ix_iy_iz[j,1]*ay + ix_iy_iz[j,2]*az
        
        final_image.arrays["ix"] = ix_iy_iz[:,0]
        final_image.arrays["iy"] = ix_iy_iz[:,1]
        final_image.arrays["iz"] = ix_iy_iz[:,2]
        final_image.positions = pos
        
        images.append(final_image)     

    energies = []
    for image in tqdm(images):
        energy = image.get_potential_energy()
        forces = image.get_forces()
        image.calc = SinglePointCalculator(image, energy=energy, forces=forces)

    for image in images:
        print(image.get_potential_energy())
    # save to file
    ase.io.write(f"{path}_final_images.xyz", images)


