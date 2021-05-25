srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density/compute_systematics_maps.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density/compute_systematics_maps.py north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density1/compute_systematics_maps.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density1/compute_systematics_maps.py north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density2/compute_systematics_maps.py south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density2/compute_systematics_maps.py north
