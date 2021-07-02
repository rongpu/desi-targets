srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python count_randoms.py north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python count_randoms.py south
