srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/compute_systematics_maps.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/compute_systematics_maps.py north

srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/test.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/test.py north

srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/test1.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -L cfs -q interactive python density_variations/test1.py north
