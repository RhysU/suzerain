#!/usr/bin/env python
"""Usage: spectra.py [OPTIONS...] H5RESTART...
Display spectra averaged across /twopoint_{kx,kz} in all named H5RESTART files.

Options:
    -c            Plot contours showing {kx,kz} versus y location
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output files *.OUTSUFFIX instead of displaying
    -p OUTPICKLE  Pickle the two point and spectra into file OUTPICKLE
    -y YINDEX     Plot spectra at wall-normal location YINDEX

Each H5RESTART should have been made by Suzerain perfect_advance (or similar),
meaning that all of the following are well-defined datasets
    /Nx, /Ny, /Nz, /Lx, /Ly, /Lz, /kx, /kz, /twopoint_kx, /twopoint_kz
adhering to a host of ill-documented restrictions.  All shapes must match!
"""
# TODO Accept normalization constants to display results in wall units

import collections
import getopt
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import numpy.fft as fft
import pickle
import sys

# Helper types providing easy result usage like plot(Ex.k, Ex.uu)
# These are the scalar field pairings to be stored in /twopoint_{kx,kz}
PAIRS = ['TT', 'Tu', 'Tv', 'Tw', 'Tr',
               'uu', 'uv', 'uw', 'ur',
                     'vv', 'vw', 'vr',
                           'ww', 'wr',
                                 'rr']
SpectralData  = collections.namedtuple('SpectralData', ['y', 'k'] + PAIRS)
PhysicalData  = collections.namedtuple('PhysicalData', ['y', 'x'] + PAIRS)
ProcessResult = collections.namedtuple('ProcessResult',
                                       ['Ekx', 'Ekz', 'Rx', 'Rz', 'Rkx', 'Rkz'])

def process(kx, kz, Lx, Lz, Nx, Nz, Rkx, Rkz, y, **kwargs):
    """Distill loaded Rkx, etc. data into easy-to-use form."""
    Rx = fft.irfft(Nx * Rkx, axis=1)  # Non-normalized inverse FFT
    Rz = fft.ifft (Nz * Rkz, axis=1)  # Non-normalized inverse FFT

    # Compute spectra from Rkx using conjugate-symmetry of Rkx
    Ekx = Rkx.copy()
    Ekx[:, 1:, :] += np.conj(Rkx[:, 1:,:])
    assert np.max(np.abs(np.imag(Ekx))) == 0
    Ekx = np.real(Ekx)

    # Compute spectra from Rkz by adding reflected negative wavenumbers
    Ekz            = Rkz[:, 0:(Nz/2+1), :].copy()
    Ekz[:, 1:, :] += Rkz[:, -1:-(Nz/2+1):-1,:]
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

    return ProcessResult(Ekx, Ekz, Rx, Rz, Rkx, Rkz)


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
    if Rkx is not None:
        Rkx /= len(h5filenames)
    if Rkz is not None:
        Rkz /= len(h5filenames)
    d.update(dict(
        Rkx = Rkx,
        Rkz = Rkz
    ))
    return d


def contour(k, y, Ek, maxk=None,
            nbins1=7, nbins2=50,
            cm1='gist_rainbow', cm2='binary'):
    """Prepare a two-dimensional contour plot a la MK's spectra presentation"""
    if not maxk: maxk = Ek.max()
    cm1  = plt.get_cmap(cm1)
    cm2  = plt.get_cmap(cm2)
    lv1  = matplotlib.ticker.MaxNLocator(nbins=nbins1).tick_values(0, maxk)
    lv2  = matplotlib.ticker.MaxNLocator(nbins=nbins2).tick_values(0, maxk)
    fig  = plt.figure()
    ax   = fig.add_subplot(111)
    c1   = ax.contour (k, y, np.transpose(Ek), levels=lv1, cmap=cm1)
    c2   = ax.contourf(k, y, np.transpose(Ek), levels=lv2, cmap=cm2)
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
    do_contour = False
    outpickle  = None
    outsuffix  = None
    yindex     = []
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "cho:p:y:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o == "-c":
                do_contour = True
            elif o in ("-h", "--help"):
                print __doc__
                return 0
            elif o == "-o":
                outsuffix = a
            elif o == "-p":
                outpickle = a
            elif o == "-y":
                yindex.append(int(a))
        if len(args) < 1:
            print >>sys.stderr, "Too few arguments.  See --help."
            return 1
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # Load and process all incoming data into easy-to-use form
    data = load(args)
    res  = process(**data)
    (Ekx, Ekz, Rx, Rz, Rkx, Rkz) = res

    # If requested, first pickle results to ease post mortem on crash
    # (it currently requires jumping through hoops to load these pickles).
    if outpickle:
        pickle.dump(res, open(outpickle, "wb"), -1)

    # Prepare spectra at each requested y index
    for j in yindex:
        # Spectra in x direction
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        ax.loglog(Ekx.k, Ekx.TT[:,j], label="T")
        ax.loglog(Ekx.k, Ekx.uu[:,j], label="u")
        ax.loglog(Ekx.k, Ekx.vv[:,j], label="v")
        ax.loglog(Ekx.k, Ekx.ww[:,j], label="w")
        ax.loglog(Ekx.k, Ekx.rr[:,j], label="rho")
        ax.set_xlabel("Streamwise wavenumber")
        ax.set_title("Spectra at y = %g" % Ekx.y[j])
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.spectra.kx.'+outsuffix)
            plt.close(fig)

        # Spectra in z direction
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        ax.loglog(Ekz.k, Ekz.TT[:,j], label="T")
        ax.loglog(Ekz.k, Ekz.uu[:,j], label="u")
        ax.loglog(Ekz.k, Ekz.vv[:,j], label="v")
        ax.loglog(Ekz.k, Ekz.ww[:,j], label="w")
        ax.loglog(Ekz.k, Ekz.rr[:,j], label="r")
        ax.set_xlabel("Spanwise wavenumber")
        ax.set_title("Spectra at y = %g" % Ekx.y[j])
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.spectra.kx.'+outsuffix)
            plt.close(fig)

        # Cospectra in x direction
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        ax.loglog(Ekx.k, Ekx.Tu[:,j], label='Tu')
        ax.loglog(Ekx.k, Ekx.Tv[:,j], label='Tv')
        ax.loglog(Ekx.k, Ekx.Tw[:,j], label='Tw')
        ax.loglog(Ekx.k, Ekx.Tr[:,j], label='Tr')
        ax.loglog(Ekx.k, Ekx.uv[:,j], label='uv')
        ax.loglog(Ekx.k, Ekx.uw[:,j], label='uw')
        ax.loglog(Ekx.k, Ekx.ur[:,j], label='ur')
        ax.loglog(Ekx.k, Ekx.vw[:,j], label='vw')
        ax.loglog(Ekx.k, Ekx.vr[:,j], label='vr')
        ax.loglog(Ekx.k, Ekx.wr[:,j], label='wr')
        ax.set_xlabel("Streamwise wavenumber")
        ax.set_title("Cospectra at y = %g" % Ekx.y[j])
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.cospectra.kx.'+outsuffix)
            plt.close(fig)

        # Cospectra in z direction
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        ax.loglog(Ekz.k, Ekz.Tu[:,j], label='Tu')
        ax.loglog(Ekz.k, Ekz.Tv[:,j], label='Tv')
        ax.loglog(Ekz.k, Ekz.Tw[:,j], label='Tw')
        ax.loglog(Ekz.k, Ekz.Tr[:,j], label='Tr')
        ax.loglog(Ekz.k, Ekz.uv[:,j], label='uv')
        ax.loglog(Ekz.k, Ekz.uw[:,j], label='uw')
        ax.loglog(Ekz.k, Ekz.ur[:,j], label='ur')
        ax.loglog(Ekz.k, Ekz.vw[:,j], label='vw')
        ax.loglog(Ekz.k, Ekz.vr[:,j], label='vr')
        ax.loglog(Ekz.k, Ekz.wr[:,j], label='wr')
        ax.set_xlabel("Streamwise wavenumber")
        ax.set_title("Cospectra at y = %g" % Ekz.y[j])
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.cospectra.kz.'+outsuffix)
            plt.close(fig)

    # Prepare contour plots for kx and kz if requested
    if do_contour:
        for index, name in enumerate(PAIRS):
            index += 2 # Hack based on knowledge of SpectralData
            (fig, ax, c1, c2, cbar) = contour(Ekx.k, Ekx.y, Ekx[index])
            ax.set_title("Spectra: " + name)
            ax.set_xlabel('Streamwise wavenumber')
            ax.set_ylabel('Wall-normal distance')
            if outsuffix:
                fig.savefig(name+'.kx.'+outsuffix)
                plt.close(fig)

        for index, name in enumerate(PAIRS):
            index += 2 # Hack based on knowledge of SpectralData
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
