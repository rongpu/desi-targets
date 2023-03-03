cd git/desi-examples/bright_star_mask/

python read_pixel_bitmask.py --tracer lrg --input /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_lrgmask_v1.1.fits.gz
python read_pixel_bitmask.py --tracer elg --input /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_elgmask_v1.fits.gz

python read_pixel_nexp.py --input /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_nexp.fits
python read_pixel_nexp.py --input /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_nexp.fits

python read_pixel_bitmask.py --tracer elg --input /pscratch/sd/r/rongpu/ebv/alternative_elg_targets.fits --output /pscratch/sd/r/rongpu/ebv/alternative_elg_targets_elgmask_v1.fits.gz

python read_pixel_bitmask.py --tracer elg --input /pscratch/sd/r/rongpu/ebv/alternative_elg_targets/alternative_elg_targets_desi_ebv.fits --output /pscratch/sd/r/rongpu/ebv/alternative_elg_targets/alternative_elg_targets_desi_ebv_elgmask_v1.fits.gz
