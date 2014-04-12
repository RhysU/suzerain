#!/usr/bin/env python
"""Usage: spectra.py [OPTIONS...] H5RESTART...
Display spectra averaged across /twopoint_{kx,kz} in all named H5RESTART files.

Options:
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output files *.OUTSUFFIX instead of displaying

Each H5RESTART should have been made by Suzerain perfect_advance (or similar),
meaning that all of the following are well-defined datasets
    /Nx, /Ny, /Nz, /Lx, /Ly, /Lz, /kx, /kz, /twopoint_kx, /twopoint_kz
adhering to a host of ill-documented restrictions.  All shapes must match!
"""
# TODO Accept normalization constants to display results in wall units
# TODO Accept "-y YINDEX     Display data from only /collocation_points_y[YINDEX]"

import collections
import getopt
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import numpy.fft as fft
import sys

# Helper types providing easy result usage like plot(Ex.k, Ex.uu)
# These are the scalar field pairings to be stored in /twopoint_{kx,kz}
PAIRS = ['TT', 'Tu', 'Tv', 'Tw', 'Tr',
               'uu', 'uv', 'uw', 'ur',
                     'vv', 'vw', 'vr',
                           'ww', 'wr',
                                 'rr']
SpectralData = collections.namedtuple('SpectralData', ['y', 'k'] + PAIRS)
PhysicalData = collections.namedtuple('PhysicalData', ['y', 'x'] + PAIRS)


def process(kx, kz, Lx, Lz, Nx, Nz, Rkx, Rkz, y, **kwargs):
    """Distill loaded Rkx, etc. data into easy-to-use form."""
    Rx = fft.irfft(Nx * Rkx, axis=1)  # Non-normalized inverse FFT
    Rz = fft.ifft(Nz * Rkz, axis=1)  # Non-normalized inverse FFT

    # Compute spectra from Rkx using conjugate-symmetry of Rkx
    Ekx = Rkx.copy()
    Ekx[:, 1:, :] += np.conj(Rkx[:, 1:,:])
    assert np.max(np.abs(np.imag(Ekx))) == 0
    Ekx = np.real(Ekx)

    # Compute spectra from Rkz by reflecting and adding the negative
    # wavenumbers
    Ekz          = Rkz[:, 0:(Nz/2), :].copy()
    Ekz[:, 1:, :] += Rkz[:, -1:-(Nz/2):-1,:]
    assert np.max(np.abs(np.imag(Ekz))) < np.finfo(Ekz.dtype).eps
    Ekz = np.real(Ekz)

    # Helper to shorten the following few statements
    def scalarpairs(d):
        return map(lambda i: np.squeeze(d[i, :,:]),
                   np.arange(d.shape[0]))

    # Prepare named view of each spectral scalar pair like 'Ex.uv' tracking y,k
    Rkx = SpectralData(y, kx,          *scalarpairs(Rkx))
    Rkz = SpectralData(y, kz,          *scalarpairs(Rkz))
    Ekx = SpectralData(y, kx,          *scalarpairs(Ekx))
    Ekz = SpectralData(y, kz[:Nz/2+1], *scalarpairs(Ekz))

    # Only half of two-point is interesting as they are periodic and even
    Rx = Rx[:, :(Rx.shape[1]/2+1), :]
    Rz = Rz[:, :(Rz.shape[1]/2+1), :]
    Rx = PhysicalData(y, np.mgrid[:Nx/2+1]*(Lx/Nx), *scalarpairs(Rx))
    Rz = PhysicalData(y, np.mgrid[:Nz/2+1]*(Lz/Nz), *scalarpairs(Rz))

    return (Ekx, Ekz, Rx, Rz, Rkx, Rkz)


def load(h5filenames):
    """Load the data required by process() into a dict from named files.
    Averages of /twopoint_kx and /twopoint_kz are taken across all inputs.
    Other results reflect only the metadata from the last file loaded.
    """
    Rkx = None
    Rkz = None
    d   = {}
    for h5filename in h5filenames:
        h5file = h5py.File(h5filename, 'r')
        d.update(dict(
            kx = h5file['kx'][()],
            kz = h5file['kz'][()],
            Lx = h5file['Lx'][0],
            Ly = h5file['Ly'][0],
            Lz = h5file['Lz'][0],
            Nx = h5file['Nx'][0],
            Ny = h5file['Ny'][0],
            Nz = h5file['Nz'][0],
            y  = h5file['collocation_points_y'][()],
        ))
        if Rkx is None:
            Rkx  = np.squeeze(h5file['twopoint_kx'][()].view(np.complex128))
        else:
            Rkx += np.squeeze(h5file['twopoint_kx'][()].view(np.complex128))
        if Rkz is None:
            Rkz  = np.squeeze(h5file['twopoint_kz'][()].view(np.complex128))
        else:
            Rkz += np.squeeze(h5file['twopoint_kz'][()].view(np.complex128))
        h5file.close()
    Rkx /= len(h5filenames)
    Rkz /= len(h5filenames)
    d.update(dict(
        Rkx = Rkx,
        Rkz = Rkz
    ))
    return d


def contour(k, y, Ek, maxk=None,
            nbins1=7, nbins2=50,
            cmap1='gist_rainbow', cmap2='binary'):
    """Prepare a two-dimensional contour plot a la MK's spectra presentation"""
    if not maxk: maxk = Ek.max()
    lv1  = matplotlib.ticker.MaxNLocator(nbins=nbins1).tick_values(0, maxk)
    lv2  = matplotlib.ticker.MaxNLocator(nbins=nbins2).tick_values(0, maxk)
    fig  = plt.figure()
    ax   = fig.add_subplot(111)
    c1   = ax.contour (k, y, Ek, levels=lv1, cmap=plt.get_cmap(cmap1))
    c2   = ax.contourf(k, y, Ek, levels=lv2, cmap=plt.get_cmap(cmap2))
    cbar = fig.colorbar(c1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    return (fig, ax, c1, c2, cbar)


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    outsuffix = None
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            elif o == "-o":
                outsuffix = a
        if len(args) < 1:
            print >>sys.stderr, "Too few arguments.  See --help."
            return 1
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # Process and distill all incoming data
    (Ekx, Ekz, Rx, Rz, Rkx, Rkz) = process(**load(args))

    # Prepare contour plots for kx
    for index, name in enumerate(PAIRS):
        (fig, ax, c1, c2, cbar) = contour(Ekx.k, Ekx.y, Ekx[index])
        ax.set_title("Spectra: " + name)
        ax.set_xlabel('Streamwise wavenumber')
        ax.set_ylabel('Wall-normal distance')
        if outsuffix:
            fig.savefig(name+'.kx.'+outsuffix)
            plt.close(fig)

    # Prepare contour plots for kx
    for index, name in enumerate(PAIRS):
        (fig, ax, c1, c2, cbar) = contour(Ekz.k, Ekz.y, Ekz[index])
        ax.set_title("Spectra: " + name)
        ax.set_xlabel('Streamwise wavenumber')
        ax.set_ylabel('Wall-normal distance')
        if outsuffix:
            fig.savefig(name+'.kz.'+outsuffix)
            plt.close(fig)

    if not outsuffix:
        plt.show()

    # Pop interactive mode
    plt.interactive(was_interactive)

if __name__ == "__main__":
    sys.exit(main())
