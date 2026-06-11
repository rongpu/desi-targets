# salloc -N 1 -C cpu -t 04:00:00 -q interactive
# cd /global/cfs/cdirs/desicollab/users/rongpu/git/desi-examples/bright_star_mask

python read_pixel_bitmask.py --tracer lrg --input /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_lrg_5.1.0_basic.fits --output /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_lrg_5.1.0_lrgmask_v1.1.fits.gz
python read_pixel_bitmask.py --tracer lrg --input /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_lge_5.1.0_basic.fits --output /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_lge_5.1.0_lrgmask_v1.1.fits.gz
python read_pixel_bitmask.py --tracer elg --input /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_elg_5.1.0_basic.fits --output /global/cfs/cdirs/desicollab/users/rongpu/targets/dr11.0/5.1.0/resolve/dr11_elg_5.1.0_elgmask_v1.fits.gz

python read_pixel_bitmask.py --tracer lrg --input /pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_0.fits --output /pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_0-lrgmask_v1.1.fits.gz
python read_pixel_bitmask.py --tracer elg --input /pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_0.fits --output /pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_0-elgmask_v1.fits.gz
