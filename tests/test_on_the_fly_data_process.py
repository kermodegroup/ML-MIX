import subprocess
import os
my_path = os.path.abspath(__file__)
import ase.io
import numpy as np


def process_data(config_path):
    script_path = os.path.join(os.path.dirname(my_path), '../mix-on-the-fly/modify_dump_otf.py')
    try:
        result = subprocess.run(['python', script_path, config_path], check=True, capture_output=True, text=True)
        print("Script output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error running script:", e.stderr)


def check_data(original_data_path, config_path, fit_elements, mask=False, pae=False, stored_force=False):
    process_data(config_path)
    processed_data = ase.io.read("./otf-fit-data.xyz", index=':')
    original_data = ase.io.read(original_data_path, index=':')
    del_index = []

    # check that mask key is present in all configs
    if mask:
        for t in processed_data:
            assert 'mask' in t.arrays.keys()
    
    # check all fully masked configs have been pruned
    if mask:
        for i in range(len(original_data)):
            mask_arr = original_data[i].arrays['f_mlml_fix']>0.99
            if np.all(~mask_arr):
                del_index.append(i)
        
        for d in del_index:
            del original_data[d]
    
    assert len(processed_data) == len(original_data)

    # check that only the fit elements are present
    for t in processed_data:
        assert np.all([symbol in fit_elements for symbol in t.get_chemical_symbols()])
    

    # check that pae key is present in all configs
    if pae:
        for t in processed_data:
            assert 'pae' in t.arrays.keys()

    # check that forces are present in all configs
    for t in processed_data:
        assert 'fit_forces' in t.arrays.keys()

    # check that the per-atom energies match
    if pae:
        for i in range(len(processed_data)):
            assert np.allclose(processed_data[i].arrays['pae'].flatten(), original_data[i].arrays['c_pae'].flatten())

    # check that the forces match
    for i in range(len(processed_data)):
        if stored_force:
            forces = np.zeros((len(original_data[0]),3))
            forces[:,0] = original_data[i].arrays['f_sf[1]'].flatten()
            forces[:,1] = original_data[i].arrays['f_sf[2]'].flatten()
            forces[:,2] = original_data[i].arrays['f_sf[3]'].flatten()
        else:
            forces = original_data[i].get_forces()
        assert np.allclose(processed_data[i].arrays['fit_forces'], forces)

    # check that mask matches
    if mask:
        for i in range(len(processed_data)):
            mask_arr = original_data[i].arrays['f_mlml_fix']>0.99
            assert np.all(processed_data[i].arrays['mask'].flatten() == mask_arr.flatten())
        
        

def test_process_data():
    check_data('on-fly/dump_data/Fe_dump0.lammpstrj', 
            'on-fly/config_files/Fe_otf_config_0.yaml',
            fit_elements=['Fe'],
            mask=True,
            stored_force=True)

    check_data('on-fly/dump_data/Fe_dump1.lammpstrj', 
            'on-fly/config_files/Fe_otf_config_1.yaml',
            fit_elements=['Fe'],
            pae=True)