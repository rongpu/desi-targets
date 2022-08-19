srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py LRG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py LRG north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py ELG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py ELG north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py QSO south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py QSO north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py BGS_ANY south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py BGS_ANY north
