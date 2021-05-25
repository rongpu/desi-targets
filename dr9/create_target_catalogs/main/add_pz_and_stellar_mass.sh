# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 LRG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 LRG north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 ELG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 ELG north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 QSO south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 QSO north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 BGS_ANY south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass.py sv3 BGS_ANY north
