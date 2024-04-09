# mask_north = zbest['PHOTSYS']=='N'
# mask_south = zbest['PHOTSYS']=='S'

# zfiber & sliding cut extension; 800 targets/sq.deg.

# South
mask_lrg = mask_south.copy()
mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
mask_lrg &= (zfibermag < 21.7)                   # faint limit

lrg_mask_sliding = rmag - w1mag > (w1mag - 17.26) * 1.8  # sliding IR cut
lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.36) * 1.  # low-z sliding IR cut
lrg_mask_sliding |= rmag - w1mag > 3.29
mask_lrg &= lrg_mask_sliding

mask_lowz = (gmag - rmag > 1.3) & ( (gmag - rmag) > -1.55 * (rmag - w1mag) + 3.13)
mask_lowz |= (rmag - w1mag > 1.8)
mask_lrg &= mask_lowz

mask_bright = (gaia_g!=0) & (gaia_g < 18)
mask_lrg &= (~mask_bright)

mask_lrg &= (zfibertotmag>17.5)

lrg_sv3_south = mask_lrg.copy()

# North
mask_lrg = mask_north.copy()
mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
mask_lrg &= (zfibermag < 21.72)                   # faint limit

lrg_mask_sliding = rmag - w1mag > (w1mag - 17.24) * 1.83  # sliding IR cut
lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.33) * 1.  # low-z sliding IR cut
lrg_mask_sliding |= rmag - w1mag > 3.39
mask_lrg &= lrg_mask_sliding

mask_lowz = (gmag - rmag > 1.34) & ( (gmag - rmag) > -1.55 * (rmag - w1mag) + 3.23)
mask_lowz |= (rmag - w1mag > 1.8)
mask_lrg &= mask_lowz

mask_bright = (gaia_g!=0) & (gaia_g < 18)
mask_lrg &= (~mask_bright)

mask_lrg &= (zfibertotmag>17.5)

lrg_sv3_north = mask_lrg.copy()

lrg_sv3 = lrg_sv3_south | lrg_sv3_north
print(np.sum(lrg_sv3))