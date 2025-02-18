import numpy as np
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import global_plot_settings
global_plot_settings.normal_text()

investigations = ['Fe_attempt','Si_attempt','W_He_attempt']
trials = ['1','2','3']
n_vals = [[16,50],[10,32],[16,50]]
atom_nos = ['8K','250K']

labels = ["Fe ACE_3_20/ACE_2_10", "Si ACE_4_20/ACE_2_10", "W-He ACE_3_20/UF3"]
all_speedups = {}

if not os.path.exists("all_data"):
    os.makedirs("all_data")

for j, investigation in enumerate(investigations):
    speedups = np.zeros((len(trials), 48))

    for k, n in enumerate(n_vals[j]):
        for i, trial in enumerate(trials):
            cheap_times = np.loadtxt(f'./{investigation}/fixed_atoms_increasing_prc/n={n}/trial_{trial}/cheap/results.txt')
            expensive_times = np.loadtxt(f'./{investigation}/fixed_atoms_increasing_prc/n={n}/trial_{trial}/expensive/results.txt')

            print(investigation, n_vals[j], trial)
            print("expensive_times", expensive_times)
            print("cheap_times", cheap_times)
            expensive_times = expensive_times[expensive_times[:, 0].argsort()]
            cheap_times = cheap_times[cheap_times[:, 0].argsort()]

            n_procs = [i for i in range(1, len(expensive_times) + 1)]
            n_steps = np.array(n_procs) * 20

            time_per_step_cheap = cheap_times[:, 1] / n_steps
            time_per_step_expensive = expensive_times[:, 1] / n_steps

            speedup_cheap = time_per_step_cheap[0] / time_per_step_cheap
            speedup_expensive = time_per_step_expensive[0] / time_per_step_expensive

            cheap_vs_expensive_speedup = time_per_step_expensive / time_per_step_cheap

            speedups[i, :len(cheap_vs_expensive_speedup)] = cheap_vs_expensive_speedup[:]

        mean_speedups = np.mean(speedups, axis=0).reshape(-1, 1)
        std_speedups = np.std(speedups, axis=0).reshape(-1, 1)
        print(np.shape(mean_speedups))
        data_vec = np.concatenate((np.array(n_procs).reshape(-1, 1), mean_speedups[:len(cheap_vs_expensive_speedup)], std_speedups[:len(cheap_vs_expensive_speedup)]), axis=1)

        np.savetxt(f"all_data/{investigation}_{atom_nos[k]}_speedups.txt", data_vec)