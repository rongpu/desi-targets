import os

import numpy as np
from scipy.ndimage import map_coordinates

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.utils import isiterable
# from astropy.config import ConfigurationItem

# SFD_MAP_DIR = ConfigurationItem('sfd_map_dir', '.',
#                                 'Directory containing SFD (1998) dust maps, '
#                                 'with names: SFD_dust_4096_[ngp,sgp].fits')

# mapdir = '/global/cfs/cdirs/desi/users/rongpu/useful/sfddata'
# mapdir = '/Users/rongpu/Documents/Data/useful/sfddata'
mapdir = '/global/cfs/cdirs/cosmo/data/dust/v0_1/maps'

def get_values_from_sfd_maps(coordinates, equatorial=True, interpolate=True, order=1, mapdir=mapdir, value='temp'):
    """Get E(B-V) value(s) from Schlegel, Finkbeiner, and Davis 1998 extinction
    maps at the given coordinates.

    Parameters
    ----------
    coordinates : astropy `~astropy.coordinates.coordsystems.SphericalCoordinatesBase` or tuple/list.
        If tuple/list, treated as (RA, Dec) in degrees the ICRS (e.g., "J2000")
        system. RA and Dec can each be float or list or numpy array.
    mapdir : str, optional
        Directory in which to find dust map FITS images, which must be named
        ``SFD_dust_4096_[ngp,sgp].fits``. If `None` (default), the value of
        the SFD_MAP_DIR configuration item is used. By default, this is ``'.'``.
        The value of SFD_MAP_DIR can be set in the configuration file,
        typically located in ``$HOME/.astropy/config/sncosmo.cfg``.
    interpolate : bool
        Interpolate between the map values using
        `scipy.ndimage.map_coordinates`.
    order : int
        Interpolation order, if interpolate=True. Default is 1.

    Returns
    -------
    v : float or `~numpy.ndarray`
        Specific extinction E(B-V) at the given locations.

    """

    # Get mapdir
    if value=='temp':
        fname = os.path.join(mapdir, 'SFD_temp_{0}.fits')
    elif value=='i100':
        fname = os.path.join(mapdir, 'SFD_i100_1024_{0}.fits')
    elif value=='d100':
        fname = os.path.join(mapdir, 'SFD_d100_1024_{0}.fits')
    elif value=='i60':
        fname = os.path.join(mapdir, 'SFD_i60_4096_{0}.fits')
    elif value=='xmap':
        fname = os.path.join(mapdir, 'SFD_xmap_{0}.fits')
    elif value=='Synch_Beta':
        fname = os.path.join(mapdir, 'Synch_Beta_{0}.fits')

    if equatorial:
        # Parse input
        # if not isinstance(coordinates, coord.SphericalCoordinatesBase):
        ra, dec = coordinates
        coordinates = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        # Convert to galactic coordinates.
        coordinates = coordinates.galactic
        l = coordinates.l.radian
        b = coordinates.b.radian
    else:
        l, b = coordinates
        l = np.radians(l)
        b = np.radians(b)

    # Check if l, b are scalar
    return_scalar = False
    if not isiterable(l):
        return_scalar = True
        l, b = np.array([l]), np.array([b])

    # Initialize return array
    v = np.empty_like(l)

    # Treat north (b>0) separately from south (b<0).
    for n, idx, ext in [(1, b >= 0, 'ngp'), (-1, b < 0, 'sgp')]:

        if not np.any(idx): continue
        hdulist = fits.open(fname.format(ext))
        mapd = hdulist[0].data

        # Project from galactic longitude/latitude to lambert pixels.
        # (See SFD98).
        npix = mapd.shape[0]        
        x = (npix / 2 * np.cos(l[idx]) * np.sqrt(1. - n*np.sin(b[idx])) +
             npix / 2 - 0.5)
        y = (-npix / 2 * n * np.sin(l[idx]) * np.sqrt(1. - n*np.sin(b[idx])) +
             npix / 2 - 0.5)
        
        # Get map values at these pixel coordinates.
        if interpolate:
            v[idx] = map_coordinates(mapd, [y, x], order=order)
        else:
            x=np.round(x).astype(np.int)
            y=np.round(y).astype(np.int)
            v[idx] = mapd[y, x]
            
        hdulist.close()
    
    if return_scalar:
        return v[0]
    return v
