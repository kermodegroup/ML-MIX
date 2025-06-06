# Simple Examples
These folders contain simple examples which aim to showcase and explain as many features of ML-MIX as possible. For CPU *all-expensive* and *ML/ML* examples are provided in each case.

For the ML/ML examples, a well-commented input script is provided alongside a minimal working one. For users who want to quickly gain a good understanding of how ML-MIX works, it is recommended they read through the commented scripts.

All examples use pure LAMMPS (i.e., no LAMMPS/Python interface), unlike the more complex workflows shown in `case_studies/`.

Two of the examples — `Fe_dumbbell_MD` and `Si_vacancy_relax` — are fleshed out and designed for **CPU-only** runs. A third example, `Si_kokkos`, is just a very simple showcase of a **GPU/Kokkos** demonstration of ML-MIX, using pair_style symmetrix/mace.

Run the CPU examples with:
```bash
mpirun -np 40 lmp -in <example_input_script>.in
```
Run the GPU/Kokkos example with:
```bash
lmp -in test_in_kk.lammps -k on g 1 -sf kk -pk kokkos newton on neigh half
```