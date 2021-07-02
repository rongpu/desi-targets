srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_systematics_maps.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python compute_systematics_maps.py north

srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python test.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python test.py north

srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python test1.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python test1.py north
