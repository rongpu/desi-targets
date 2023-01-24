# salloc -N 1 -C cpu -q interactive -t 04:00:00
python rf_photoz_with_magnification.py north 0.99 &> /global/u2/r/rongpu/temp/tmp/north_0.99.log
python rf_photoz_with_magnification.py north 1. &> /global/u2/r/rongpu/temp/tmp/north_1.log
python rf_photoz_with_magnification.py north 1.01 &> /global/u2/r/rongpu/temp/tmp/north_1.01.log
