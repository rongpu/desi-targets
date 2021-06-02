srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density/count_randoms.py north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density/count_randoms.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density1/count_randoms.py north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density1/count_randoms.py south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density2/count_randoms.py north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python density2/count_randoms.py south

# python density/count_randoms.py north
# python density/count_randoms.py south
# python density1/count_randoms.py north
# python density1/count_randoms.py south
# python density2/count_randoms.py north
# python density2/count_randoms.py south
