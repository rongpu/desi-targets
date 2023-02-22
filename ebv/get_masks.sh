cd git/desi-examples/bright_star_mask/

srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer lrg --input /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer elg --input /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_elgmask_v1.fits.gz

srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python read_pixel_nexp.py --input /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_nexp.fits
srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python read_pixel_nexp.py --input /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_nexp.fits
