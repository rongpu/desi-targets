cd /global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_bitmask.py --tracer lrg --input main_lrg_magnification_north.fits --output main_lrg_magnification_north_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_bitmask.py --tracer lrg --input main_lrg_magnification_south.fits --output main_lrg_magnification_south_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_bitmask.py --tracer lrg --input extended_lrg_magnification_north.fits --output extended_lrg_magnification_north_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_bitmask.py --tracer lrg --input extended_lrg_magnification_south.fits --output extended_lrg_magnification_south_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_nexp.py --input extended_lrg_magnification_north.fits --output extended_lrg_magnification_north_pixel.fits
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_nexp.py --input extended_lrg_magnification_south.fits --output extended_lrg_magnification_south_pixel.fits

cd /global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_bitmask.py --tracer lrg --input dr9_extended_lrg_0.49.0_basic.fits --output dr9_extended_lrg_0.49.0_lrgmask_v1.1.fits.gz
srun -N 1 -C cpu -c 128 -t 04:00:00 -q interactive python ~/git/desi-examples/bright_star_mask/read_pixel_nexp.py --input dr9_extended_lrg_0.49.0_basic.fits --output dr9_extended_lrg_0.49.0_pixel.fits
