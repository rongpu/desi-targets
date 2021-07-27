srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 LRG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 LRG north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 BGS_ANY south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 BGS_ANY north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 ELG south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 ELG north
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 QSO south
srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py sv1 QSO north

# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main ELG south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main ELG north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main LRG south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main LRG north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main QSO south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main QSO north
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main BGS_ANY south
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_sweep_columns_to_target_catalogs.py main BGS_ANY north
