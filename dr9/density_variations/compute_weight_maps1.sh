# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py BGS_ANY north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py BGS_ANY south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py BGS_BRIGHT north
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py BGS_BRIGHT south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py LRG north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py LRG south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py ELG north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py QSO north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py QSO south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_weight_maps.py ELG south
