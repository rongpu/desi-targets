from __future__ import division, print_function
import numpy as np
from astropy.table import Table


def parse_custom_mask(custom_mask_fn='desi_custom_mask_v1.txt'):
    '''
    Parse the custom mask file and return astropy tables.

    Args:
        custom_mask_fn: path to the custom mask file
    Returns:
        circ_mask and rect_mask: circular and rectangular masks in astropy table format
    '''

    with open(custom_mask_fn, 'r') as f:
        lines = list(map(str.strip, f.readlines()))

    circ_mask_arr = []
    rect_mask_arr = []

    for line in lines:
        if line!='' and line[0]!='#':
            line = line[:line.find('#')]
            line = list(map(float, line.split(',')))
            if len(line)==3:
                circ_mask_arr.append(line)
            elif len(line)==4:
                rect_mask_arr.append(line)
            else:
                raise ValueError

    circ_mask_arr = np.array(circ_mask_arr)
    rect_mask_arr = np.array(rect_mask_arr)

    circ_mask = Table(circ_mask_arr, names=['ra', 'dec', 'radius'])
    rect_mask = Table(rect_mask_arr, names=['ramin', 'ramax', 'decmin', 'decmax'])

    return circ_mask, rect_mask
