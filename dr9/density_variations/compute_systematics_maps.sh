# salloc -N 1 -C cpu -q interactive -t 4:00:00

# python compute_systematics_maps.py south
# python compute_systematics_maps.py north

python test.py south
python test.py north

python test1.py south
python test1.py north
