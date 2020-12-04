from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import sys, os, glob, time, warnings


def plot_comoving_density(z_phot, area, zmin=0.25, zmax=1.25, dz=0.01, 
    axis=[0.25, 1.25, 0, 0.85], label=None, color=None, lw=None, ls=None, show=True):

    sys.path.append(os.path.expanduser("~")+'/git/desi-lrg-selection')
    from modules import comoving_density

    nbins = int(round((zmax-zmin)/dz))+1
    bins = np.linspace(zmin, zmax, nbins)
    bin_center = (bins[1:]+bins[:-1])/2

    densities = comoving_density(z_phot, bins, area)

    plt.hist(bin_center, bins=bins, weights=densities*1000, histtype='step', label=label, 
        color=color, lw=lw, ls=ls)
    plt.xlabel('Photo-z')
    plt.ylabel('Comoving density * 1000 $(h^3\mathrm{Mpc}^{-3})$')
    plt.grid(alpha=0.2)
    if axis!='auto':
        plt.axis(axis)
    # plt.tight_layout()
    if label is not None:
        plt.legend(loc='upper right')
    if show:
        plt.show()


def density_plot(cat, selection_params=None, mask=None, cmin=None, cmax=None, title=None, lognorm=True, normed=False):

    from matplotlib.colors import LogNorm

    if mask is None:
        mask = np.ones(len(cat), dtype=bool)

    kwargs = {}
    if lognorm:
        kwargs['norm'] = LogNorm()

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    with plt.style.context(('dark_background')):

        # non-stellar cut
        bins = [np.linspace(0.2, 2.8, 100), np.linspace(-1, 3, 100)]
        plt.figure(figsize=(8, 6))
        plt.hist2d((cat['rmag']-cat['zmag'])[mask], (cat['zmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        plt.xlabel('r - z')
        plt.ylabel('z - W1')
        plt.colorbar()
        # plt.axis([0.2, 2.8, -1, 3])
        if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
            x = np.linspace(0, 3)
            y = selection_params['ns_a'] * x + selection_params['ns_b']
            plt.plot(x, y, 'r--', lw=1)
        if title is not None:
            plt.title(title)
        plt.show()

        # r-z and g-r cut
        plt.figure(figsize=(8, 6))
        bins = [np.linspace(0.5, 2.55, 100), np.linspace(-0.5, 3.5, 100)]
        plt.hist2d((cat['rmag']-cat['zmag'])[mask], (cat['gmag']-cat['rmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        plt.colorbar()
        # plt.axis([0.5, 2.55, -0.5, 3.5])
        plt.xlabel('r - z')
        plt.ylabel('g - r')
        if title is not None:
            plt.title(title)
        plt.show()

        # optical sliding cut
        plt.figure(figsize=(8, 6))
        bins = [np.linspace(17, 22.5, 100), np.linspace(0.5, 2.6, 100)]
        plt.hist2d((cat['zmag'])[mask], (cat['rmag']-cat['zmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        plt.colorbar()
        # plt.axis([17, 22.5, 0.5, 2.6])
        plt.xlabel('z')
        plt.ylabel('r - z')
        x = np.linspace(17, 22)
        if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
            y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
            plt.plot(x, y1, 'r--', lw=1)
        if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
            y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
            plt.plot(x, y2, 'r--', lw=1)
        if title is not None:
            plt.title(title)
        plt.show()

        # IR sliding cut
        plt.figure(figsize=(8, 6))
        bins = [np.linspace(16.3, 20.2, 100), np.linspace(0.3, 5., 100)]
        plt.hist2d((cat['w1mag'])[mask], (cat['rmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        plt.colorbar()
        # plt.axis([16.3, 20.2, 0.3, 5.])
        plt.xlabel('W1')
        plt.ylabel('r - W1')
        if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
            x = np.linspace(16, 21)
            y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
            plt.plot(x, y1, 'r--', lw=1)
        if title is not None:
            plt.title(title)
        plt.show()


def density_plot_2x2(cat, selection_params=None, mask=None, cmin=None, cmax=None, title=None, show=True, lognorm=True, normed=False):

    from matplotlib.colors import LogNorm

    if mask is None:
        mask = np.ones(len(cat), dtype=bool)

    kwargs = {}
    if lognorm:
        kwargs['norm'] = LogNorm()

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    with plt.style.context(('dark_background')):

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # non-stellar cut
        bins = [np.linspace(0.2, 2.8, 100), np.linspace(-1, 3, 100)]
        im = axes[0, 0].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['zmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        axes[0, 0].set_xlabel('r - z')
        axes[0, 0].set_ylabel('z - W1')
        # axes[0, 0].set_facecolor('k')
        fig.colorbar(im[3], ax=axes[0, 0])
        # axes[0, 0].axis([0.2, 2.8, -1, 3])
        if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
            x = np.linspace(0, 3)
            y = selection_params['ns_a'] * x + selection_params['ns_b']
            axes[0, 0].plot(x, y, 'r--', lw=1)

        # r-z and g-r cut
        bins = [np.linspace(0.5, 2.55, 100), np.linspace(-0.5, 3.5, 100)]
        im = axes[0, 1].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['gmag']-cat['rmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[0, 1])
        # axes[0, 1].axis([0.5, 2.55, -0.5, 3.5])
        axes[0, 1].set_xlabel('r - z')
        axes[0, 1].set_ylabel('g - r')
        # axes[0, 1].set_facecolor('k')

        # optical sliding cut
        bins = [np.linspace(17, 22.5, 100), np.linspace(0.5, 2.6, 100)]
        im = axes[1, 0].hist2d((cat['zmag'])[mask], (cat['rmag']-cat['zmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 0])
        # axes[1, 0].axis([17, 22.5, 0.5, 2.6])
        axes[1, 0].set_xlabel('z')
        axes[1, 0].set_ylabel('r - z')
        # axes[1, 0].set_facecolor('k')
        x = np.linspace(17, 22)
        if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
            y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y1, 'r--', lw=1)
        if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
            y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y2, 'r--', lw=1)

        # IR sliding cut
        bins = [np.linspace(16.3, 20.2, 100), np.linspace(0.3, 5., 100)]
        im = axes[1, 1].hist2d((cat['w1mag'])[mask], (cat['rmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 1])
        # axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
        axes[1, 1].set_xlabel('W1')
        axes[1, 1].set_ylabel('r - W1')
        # axes[1, 1].set_facecolor('k')
        if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
            x = np.linspace(16, 21)
            y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
            axes[1, 1].plot(x, y1, 'r--', lw=1)

        if title is not None:
            fig.suptitle(title)
        fig.subplots_adjust(top=0.945)

    if show:
        plt.show()
    else:
        return fig, axes



def color_plot(cat, selection_params=None, idx=None, vmin=0.2, vmax=1.0, cmap='Dark2_r', title=None):

    if idx is None:
        idx = np.arange(len(cat))

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    # non-stellar cut
    plt.figure(figsize=(8, 6))
    plt.scatter((cat['rmag']-cat['zmag'])[idx], (cat['zmag']-cat['w1mag'])[idx], 
                c=cat['z_phot'][idx], edgecolor='', s=1.0, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.xlabel('r - z')
    plt.ylabel('z - W1')
    plt.axis([0.2, 2.8, -1, 3])
    if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
        x = np.linspace(0, 3)
        y = selection_params['ns_a'] * x + selection_params['ns_b']
        plt.plot(x, y, 'r--', lw=1)
    if title is not None:
        plt.title(title)
    plt.show()

    # r-z and g-r cut
    plt.figure(figsize=(8, 6))
    plt.scatter((cat['rmag']-cat['zmag'])[idx], (cat['gmag']-cat['rmag'])[idx], 
             c=cat['z_phot'][idx], edgecolor='', s=1.0, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.axis([0.5, 2.55, -0.5, 3.5])
    plt.xlabel('r - z')
    plt.ylabel('g - r')
    if title is not None:
        plt.title(title)
    plt.show()

    # optical sliding cut
    plt.figure(figsize=(8, 6))
    plt.scatter((cat['zmag'])[idx], (cat['rmag']-cat['zmag'])[idx], 
             c=cat['z_phot'][idx], edgecolor='', s=1.0, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.axis([17, 22.5, 0.5, 2.6])
    plt.xlabel('z')
    plt.ylabel('r - z')
    x = np.linspace(17, 22)
    if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
        y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
        plt.plot(x, y1, 'r--', lw=1)
    if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
        y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
        plt.plot(x, y2, 'r--', lw=1)
    if title is not None:
        plt.title(title)
    plt.show()

    # IR sliding cut
    plt.figure(figsize=(8, 6))
    plt.scatter((cat['w1mag'])[idx], (cat['rmag']-cat['w1mag'])[idx], 
             c=cat['z_phot'][idx], edgecolor='', s=1.0, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.axis([16.3, 20.2, 0.3, 5.])
    plt.xlabel('W1')
    plt.ylabel('r - W1')
    if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
        x = np.linspace(16, 21)
        y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
        plt.plot(x, y1, 'r--', lw=1)
    if title is not None:
        plt.title(title)
    plt.show()


def color_plot_2x2(cat, selection_params=None, idx=None, title=None, vmin=0.2, vmax=1.0, cmap='Dark2_r', show=True):

    if idx is None:
        idx = np.arange(len(cat))

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # non-stellar cut
    im = axes[0, 0].scatter((cat['rmag']-cat['zmag'])[idx], (cat['zmag']-cat['w1mag'])[idx], 
                c=cat['z_phot'][idx], s=0.1, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_xlabel('r - z')
    axes[0, 0].set_ylabel('z - W1')
    axes[0, 0].axis([0.2, 2.8, -1, 3])
    fig.colorbar(im, ax=axes[0, 0])
    if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
        x = np.linspace(0, 3)
        y = selection_params['ns_a'] * x + selection_params['ns_b']
        plt.plot(x, y, 'r--', lw=1)

    # r-z and g-r cut
    im = axes[0, 1].scatter((cat['rmag']-cat['zmag'])[idx], (cat['gmag']-cat['rmag'])[idx], 
             c=cat['z_phot'][idx], s=0.1, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_xlabel('r - z')
    axes[0, 1].set_ylabel('g - r')
    axes[0, 1].axis([0.5, 2.55, -0.5, 3.5])
    fig.colorbar(im, ax=axes[0, 1])

    # optical sliding cut
    im = axes[1, 0].scatter((cat['zmag'])[idx], (cat['rmag']-cat['zmag'])[idx], 
             c=cat['z_phot'][idx], s=0.1, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_xlabel('z')
    axes[1, 0].set_ylabel('r - z')
    axes[1, 0].axis([17, 22.5, 0.5, 2.6])
    fig.colorbar(im, ax=axes[1, 0])
    x = np.linspace(17, 22)
    if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
        y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
        plt.plot(x, y1, 'r--', lw=1)
    if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
        y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
        plt.plot(x, y2, 'r--', lw=1)

    # IR sliding cut
    im = axes[1, 1].scatter((cat['w1mag'])[idx], (cat['rmag']-cat['w1mag'])[idx], 
             c=cat['z_phot'][idx], s=0.1, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 1].set_xlabel('W1')
    axes[1, 1].set_ylabel('r - W1')
    axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
    fig.colorbar(im, ax=axes[1, 1])
    if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
        x = np.linspace(16, 21)
        y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
        plt.plot(x, y1, 'r--', lw=1)
    
    if title is not None:
        fig.suptitle(title)
    fig.subplots_adjust(top=0.945)
    
    if show:
        plt.show()
    else:
        return fig, axes


def density_plot_3x2(cat, selection_params=None, mask=None, cmin=None, cmax=None, title=None, show=True, lognorm=True, normed=False):

    from matplotlib.colors import LogNorm

    if mask is None:
        mask = np.ones(len(cat), dtype=bool)

    kwargs = {}
    if lognorm:
        kwargs['norm'] = LogNorm()

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    with plt.style.context(('dark_background')):

        fig, axes = plt.subplots(3, 2, figsize=(18, 20))

        # non-stellar cut
        bins = [np.linspace(0.2, 2.8, 100), np.linspace(-1, 3, 100)]
        im = axes[0, 0].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['zmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        axes[0, 0].set_xlabel('r - z')
        axes[0, 0].set_ylabel('z - W1')
        # axes[0, 0].set_facecolor('k')
        fig.colorbar(im[3], ax=axes[0, 0])
        # axes[0, 0].axis([0.2, 2.8, -1, 3])
        if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
            x = np.linspace(0, 3)
            y = selection_params['ns_a'] * x + selection_params['ns_b']
            axes[0, 0].plot(x, y, 'r--', lw=1)

        # g-r vs r-z
        bins = [np.linspace(0.2, 2.4, 100), np.linspace(0, 2.7, 100)]
        im = axes[0, 1].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['gmag']-cat['rmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[0, 1])
        # axes[0, 1].axis([0.5, 2.55, -0.5, 3.5])
        axes[0, 1].set_xlabel('r - z')
        axes[0, 1].set_ylabel('g - r')
        # axes[0, 1].set_facecolor('k')

        # optical sliding cut
        bins = [np.linspace(17, 22.5, 100), np.linspace(0.5, 2.6, 100)]
        im = axes[1, 0].hist2d((cat['zmag'])[mask], (cat['rmag']-cat['zmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 0])
        # axes[1, 0].axis([17, 22.5, 0.5, 2.6])
        axes[1, 0].set_xlabel('z')
        axes[1, 0].set_ylabel('r - z')
        # axes[1, 0].set_facecolor('k')
        x = np.linspace(17, 22)
        if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
            y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y1, 'r--', lw=1)
        if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
            y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y2, 'r--', lw=1)

        # IR sliding cut
        bins = [np.linspace(16.3, 20.2, 100), np.linspace(0.3, 5., 100)]
        im = axes[1, 1].hist2d((cat['w1mag'])[mask], (cat['rmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 1])
        # axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
        axes[1, 1].set_xlabel('W1')
        axes[1, 1].set_ylabel('r - W1')
        # axes[1, 1].set_facecolor('k')
        if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
            x = np.linspace(16, 21)
            y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
            axes[1, 1].plot(x, y1, 'r--', lw=1)

        # zfiber vs z
        bins = [np.linspace(17, 22.5, 100), np.linspace(18, 23.5, 100)]
        im = axes[2, 0].hist2d((cat['zmag'])[mask], (cat['zfibermag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[2, 0])
        # axes[2, 0].axis([17, 22.5, 18, 23.5])
        axes[2, 0].set_xlabel('z')
        axes[2, 0].set_ylabel('zfiber')
        # axes[2, 0].set_facecolor('k')

        axes[2, 1].axis('off')

        if title is not None:
            fig.suptitle(title)
        fig.subplots_adjust(top=0.96)

    if show:
        plt.show()
    else:
        return fig, axes


def density_plot_3x2_new(cat, selection_params=None, mask=None, cmin=None, cmax=None, title=None, show=True, lognorm=False, normed=False):
    '''
    Replace the g-r vs r-z panel with g-r vs r-W1
    '''

    from matplotlib.colors import LogNorm

    if mask is None:
        mask = np.ones(len(cat), dtype=bool)

    kwargs = {}
    if lognorm:
        kwargs['norm'] = LogNorm()

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    with plt.style.context(('dark_background')):

        fig, axes = plt.subplots(3, 2, figsize=(18, 20))

        # non-stellar cut
        bins = [np.linspace(-0.5, 1, 201), np.linspace(-1, 2.5, 201)]
        im = axes[0, 0].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['zmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        axes[0, 0].set_xlabel('r - z')
        axes[0, 0].set_ylabel('z - W1')
        # axes[0, 0].set_facecolor('k')
        fig.colorbar(im[3], ax=axes[0, 0])
        # axes[0, 0].axis([0.2, 2.8, -1, 3])
        if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
            x = np.linspace(0, 3)
            y = selection_params['ns_a'] * x + selection_params['ns_b']
            axes[0, 0].plot(x, y, 'r--', lw=1)

        # g-r vs r-z
        bins = [np.linspace(-0.5, 1, 201), np.linspace(-0.5, 1, 201)]
        im = axes[0, 1].hist2d((cat['rmag']-cat['zmag'])[mask], (cat['gmag']-cat['rmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[0, 1])
        # axes[0, 1].axis([0.5, 2.55, -0.2, 4.])
        axes[0, 1].set_xlabel('r - z')
        axes[0, 1].set_ylabel('g - r')
        # axes[0, 1].set_facecolor('k')

        # optical sliding cut
        bins = [np.linspace(17.5, 23.5, 201), np.linspace(-1, 1.5, 201)]
        im = axes[1, 0].hist2d((cat['zmag'])[mask], (cat['rmag']-cat['zmag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 0])
        # axes[1, 0].axis([17, 22.5, 0.5, 2.6])
        axes[1, 0].set_xlabel('z')
        axes[1, 0].set_ylabel('r - z')
        # axes[1, 0].set_facecolor('k')
        x = np.linspace(17, 22)
        if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
            y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y1, 'r--', lw=1)
        if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
            y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
            axes[1, 0].plot(x, y2, 'r--', lw=1)

        # IR sliding cut
        bins = [np.linspace(16, 23.5, 201), np.linspace(-1, 4., 201)]
        im = axes[1, 1].hist2d((cat['w1mag'])[mask], (cat['rmag']-cat['w1mag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[1, 1])
        # axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
        axes[1, 1].set_xlabel('W1')
        axes[1, 1].set_ylabel('r - W1')
        # axes[1, 1].set_facecolor('k')
        if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
            x = np.linspace(16, 21)
            y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
            axes[1, 1].plot(x, y1, 'r--', lw=1)

        # zfiber vs z
        bins = [np.linspace(17, 22., 201), np.linspace(18, 23., 201)]
        im = axes[2, 0].hist2d((cat['zmag'])[mask], (cat['zfibermag'])[mask], 
            bins=bins, cmin=cmin, cmax=cmax, density=normed, **kwargs)
        fig.colorbar(im[3], ax=axes[2, 0])
        # axes[2, 0].axis([17, 22.5, 18, 23.5])
        axes[2, 0].set_xlabel('z')
        axes[2, 0].set_ylabel('zfiber')
        # axes[2, 0].set_facecolor('k')

        axes[2, 1].axis('off')

        if title is not None:
            fig.suptitle(title)
        fig.subplots_adjust(top=0.96)

    if show:
        plt.show()
    else:
        return fig, axes
    

def color_plot_3x2(cat, selection_params=None, idx=None, title=None, vmin=0.2, vmax=1.0, cmap='Dark2_r', ms=0.3, show=True):

    if idx is None:
        idx = np.arange(len(cat))

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    fig, axes = plt.subplots(3, 2, figsize=(18, 20))

    # non-stellar cut
    im = axes[0, 0].scatter((cat['rmag']-cat['zmag'])[idx], (cat['zmag']-cat['w1mag'])[idx], 
                c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_xlabel('r - z')
    axes[0, 0].set_ylabel('z - W1')
    axes[0, 0].axis([0.2, 2.8, -1, 3])
    fig.colorbar(im, ax=axes[0, 0])
    if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
        x = np.linspace(0, 3)
        y = selection_params['ns_a'] * x + selection_params['ns_b']
        plt.plot(x, y, 'r--', lw=1)

    # g-r vs r-z
    im = axes[0, 1].scatter((cat['rmag']-cat['zmag'])[idx], (cat['gmag']-cat['rmag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_xlabel('r - z')
    axes[0, 1].set_ylabel('g - r')
    axes[0, 1].axis([0.2, 2.4, 0, 2.7])
    fig.colorbar(im, ax=axes[0, 1])

    # optical sliding cut
    im = axes[1, 0].scatter((cat['zmag'])[idx], (cat['rmag']-cat['zmag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_xlabel('z')
    axes[1, 0].set_ylabel('r - z')
    axes[1, 0].axis([17, 22.5, 0.5, 2.6])
    fig.colorbar(im, ax=axes[1, 0])
    x = np.linspace(17, 22)
    if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
        y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
        plt.plot(x, y1, 'r--', lw=1)
    if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
        y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
        plt.plot(x, y2, 'r--', lw=1)

    # IR sliding cut
    im = axes[1, 1].scatter((cat['w1mag'])[idx], (cat['rmag']-cat['w1mag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 1].set_xlabel('W1')
    axes[1, 1].set_ylabel('r - W1')
    axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
    fig.colorbar(im, ax=axes[1, 1])
    if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
        x = np.linspace(16, 21)
        y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
        plt.plot(x, y1, 'r--', lw=1)
    
    # zfiber vs z
    im = axes[2, 0].scatter((cat['zmag'])[idx], (cat['zfibermag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[2, 0].set_xlabel('z')
    axes[2, 0].set_ylabel('zfiber')
    axes[2, 0].axis([17, 22.5, 18, 23.5])
    fig.colorbar(im, ax=axes[2, 0])

    axes[2, 1].axis('off')

    if title is not None:
        fig.suptitle(title)
    fig.subplots_adjust(top=0.96)
    
    if show:
        plt.show()
    else:
        return fig, axes
    

def color_plot_3x2_new(cat, selection_params=None, idx=None, title=None, vmin=0.2, vmax=1.0, cmap='Dark2_r', ms=0.3, show=True):
    '''
    Replace the g-r vs r-z panel with g-r vs r-W1
    '''
    if idx is None:
        idx = np.arange(len(cat))

    if (selection_params is not None) and ('opt_rz_slope' not in selection_params):
        selection_params['opt_rz_slope'] = 2.

    fig, axes = plt.subplots(3, 2, figsize=(18, 20))

    # non-stellar cut
    im = axes[0, 0].scatter((cat['rmag']-cat['zmag'])[idx], (cat['zmag']-cat['w1mag'])[idx], 
                c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_xlabel('r - z')
    axes[0, 0].set_ylabel('z - W1')
    axes[0, 0].axis([0.2, 2.8, -1, 3])
    fig.colorbar(im, ax=axes[0, 0])
    if (selection_params is not None) and (selection_params['ns_a'] is not None) and (selection_params['ns_b'] is not None):
        x = np.linspace(0, 3)
        y = selection_params['ns_a'] * x + selection_params['ns_b']
        plt.plot(x, y, 'r--', lw=1)

    # g-r vs r-W1
    im = axes[0, 1].scatter((cat['rmag']-cat['w1mag'])[idx], (cat['gmag']-cat['rmag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_xlabel('r - w1')
    axes[0, 1].set_ylabel('g - r')
    axes[0, 1].axis([-0.2, 4.5, 0, 2.7])
    fig.colorbar(im, ax=axes[0, 1])

    # optical sliding cut
    im = axes[1, 0].scatter((cat['zmag'])[idx], (cat['rmag']-cat['zmag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_xlabel('z')
    axes[1, 0].set_ylabel('r - z')
    axes[1, 0].axis([17, 22.5, 0.5, 2.6])
    fig.colorbar(im, ax=axes[1, 0])
    x = np.linspace(17, 22)
    if (selection_params is not None) and (selection_params['opt_rz1'] is not None):
        y1 = (x + selection_params['opt_rz1'])/selection_params['opt_rz_slope']
        plt.plot(x, y1, 'r--', lw=1)
    if (selection_params is not None) and (selection_params['opt_rz2'] is not None):
        y2 = (x + selection_params['opt_rz2'])/selection_params['opt_rz_slope']
        plt.plot(x, y2, 'r--', lw=1)

    # IR sliding cut
    im = axes[1, 1].scatter((cat['w1mag'])[idx], (cat['rmag']-cat['w1mag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 1].set_xlabel('W1')
    axes[1, 1].set_ylabel('r - W1')
    axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
    fig.colorbar(im, ax=axes[1, 1])
    if (selection_params is not None) and (selection_params['ir_a'] is not None) and (selection_params['ir_b'] is not None):
        x = np.linspace(16, 21)
        y1 = (x + selection_params['ir_a'])/selection_params['ir_b']
        plt.plot(x, y1, 'r--', lw=1)
    
    # zfiber vs z
    im = axes[2, 0].scatter((cat['zmag'])[idx], (cat['zfibermag'])[idx], 
             c=cat['z_phot'][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[2, 0].set_xlabel('z')
    axes[2, 0].set_ylabel('zfiber')
    axes[2, 0].axis([17, 22., 18, 23.])
    fig.colorbar(im, ax=axes[2, 0])

    axes[2, 1].axis('off')

    if title is not None:
        fig.suptitle(title)
    fig.subplots_adjust(top=0.96)
    
    if show:
        plt.show()
    else:
        return fig, axes

def color_plot_3x2_simple(cat, color_col, idx=None, title=None, vmin=0.2, vmax=1.0, cmap='Dark2_r', ms=0.3, show=True, figaxis=None, colorbar=True, figsize=(18, 20), grid=False):
    '''
    Replace the g-r vs r-z panel with g-r vs r-W1
    '''
    if idx is None:
        idx = np.arange(len(cat))

    if figaxis is None:
        fig, axes = plt.subplots(3, 2, figsize=figsize)
    else:
        fig, axes = figaxis

    # non-stellar cut
    im = axes[0, 0].scatter((cat['rmag']-cat['zmag'])[idx], (cat['zmag']-cat['w1mag'])[idx], 
                c=cat[color_col][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_xlabel('r - z')
    axes[0, 0].set_ylabel('z - W1')
    axes[0, 0].axis([0.2, 2.8, -1, 3])
    if colorbar:
        fig.colorbar(im, ax=axes[0, 0])
    if grid:
        axes[0, 0].grid(alpha=0.5)

    # g-r vs r-W1
    im = axes[0, 1].scatter((cat['rmag']-cat['w1mag'])[idx], (cat['gmag']-cat['rmag'])[idx], 
             c=cat[color_col][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_xlabel('r - w1')
    axes[0, 1].set_ylabel('g - r')
    axes[0, 1].axis([-0.2, 4.5, 0, 2.7])
    if colorbar:
        fig.colorbar(im, ax=axes[0, 1])
    if grid:
        axes[0, 1].grid(alpha=0.5)

    # optical sliding cut
    im = axes[1, 0].scatter((cat['zmag'])[idx], (cat['rmag']-cat['zmag'])[idx], 
             c=cat[color_col][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_xlabel('z')
    axes[1, 0].set_ylabel('r - z')
    axes[1, 0].axis([17, 22., 0.5, 2.6])
    if colorbar:
        fig.colorbar(im, ax=axes[1, 0])
    if grid:
        axes[1, 0].grid(alpha=0.5)
    x = np.linspace(17, 22)

    # IR sliding cut
    im = axes[1, 1].scatter((cat['w1mag'])[idx], (cat['rmag']-cat['w1mag'])[idx], 
             c=cat[color_col][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 1].set_xlabel('W1')
    axes[1, 1].set_ylabel('r - W1')
    axes[1, 1].axis([16.3, 20.2, 0.3, 5.])
    if colorbar:
        fig.colorbar(im, ax=axes[1, 1])
    if grid:
        axes[1, 1].grid(alpha=0.5)
    
    # zfiber vs z
    im = axes[2, 0].scatter((cat['zmag'])[idx], (cat['zfibermag'])[idx], 
             c=cat[color_col][idx], s=ms, cmap=cmap, vmin=vmin, vmax=vmax)
    axes[2, 0].set_xlabel('z')
    axes[2, 0].set_ylabel('zfiber')
    axes[2, 0].axis([17, 22., 18, 23.])
    if colorbar:
        fig.colorbar(im, ax=axes[2, 0])
    if grid:
        axes[2, 0].grid(alpha=0.5)

    axes[2, 1].axis('off')

    if title is not None:
        fig.suptitle(title)
    fig.subplots_adjust(top=0.96)
    
    if show:
        plt.show()
    else:
        return fig, axes