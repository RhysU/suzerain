#!/usr/bin/env python
"""Usage: spectra.py [OPTIONS...] H5RESTART...
Display spectra averaged across /twopoint_{kx,kz} in all named H5RESTART files.

Options:
    -c            Plot contours showing {kx,kz} versus y location
    -h            Display this help message and exit
    -o OUTSUFFIX  Save the output files *.OUTSUFFIX instead of displaying
    -p OUTPICKLE  Pickle the two point and spectra into file OUTPICKLE
    -P INPICKLE   Load spectra from INPICKLE and not H5RESTART arguments
    -v            Increase verbosity, including status of H5RESTART loads
    -y YINDEX     Plot spectra at wall-normal location YINDEX

Each H5RESTART should have been made by Suzerain perfect_advance (or similar),
meaning that all of the following are well-defined datasets
    /Nx, /Ny, /Nz, /Lx, /Ly, /Lz, /kx, /kz, /twopoint_kx, /twopoint_kz,
and if present, /antioch_constitutive_data,
adhering to a host of ill-documented restrictions.  All shapes must match!
"""
# TODO Accept normalization constants to display results in wall units

from collections import namedtuple
import gb
import getopt
import h5py
import logging
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import numpy.fft as fft
import os
import pickle
import sys


# Prepare a logger to produce messages
logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
log = logging.getLogger(os.path.splitext(os.path.basename(__file__))[0])


# Primitive variables expected to be found in every H5RESTART file
primitive, primitive_pairs = ['T', 'u', 'v', 'w', 'r'], []
for i in xrange(0, len(primitive)):
    for j in xrange(i, len(primitive)):
        primitive_pairs.append(primitive[i] + primitive[j])


# Define, at the module level, SpectralData and friends for primitive-only case
# These will be redefined in process(...), if necessary, for reacting data but
# this definition permits pickling and unpickling to work on non-reacting data.
SpectralData  = namedtuple('SpectralData', ['y', 'k'] + primitive_pairs)
PhysicalData  = namedtuple('PhysicalData', ['y', 'x'] + primitive_pairs)


# The return type from the process(...) method just below.
ProcessResult = namedtuple('ProcessResult',
                           ['Ekx', 'Ekz', 'Rx', 'Rz', 'Rkx', 'Rkz', 'bar'])


def process(kx, kz, Lx, Lz, Nx, Nz, Rkx, Rkz, bar, y, Ns, sn, **kwargs):
    """Distill loaded Rkx, etc. data into easy-to-use form."""

    # Build up the list of variables and then pairwise products of interest
    vars = primitive[:]
    if sn:
        for s in sn:
            vars.append('c'+s)
    pairs = []
    bar_pairs = np.empty(( (len(vars)*(len(vars)+1))/2, len(y)) )
    for i in xrange(0, len(vars)):
        for j in xrange(i, len(vars)):
            bar_pairs[len(pairs),:] = np.multiply(bar[i,:],bar[j,:])
            pairs.append(vars[i] + vars[j])

    # When any species are present, we must redefine namedtuple classes as the
    # default versions provided above are insufficient due to species data.
    # An atrocious hack to permit (un)pickling simple, inert data.  Sorry.
    global SpectralData, PhysicalData
    if sn:
        SpectralData = namedtuple('SpectralData', ['y', 'k'] + pairs)
        PhysicalData = namedtuple('PhysicalData', ['y', 'x'] + pairs)
    else:
        assert(pairs == primitive_pairs)

    Rx = fft.irfft(Nx * Rkx, axis=1)  # Non-normalized inverse FFT
    Rz = fft.ifft (Nz * Rkz, axis=1)  # Non-normalized inverse FFT

    # Compute spectra from Rkx using conjugate-symmetry of Rkx
    Ekx = Rkx.copy()
    Ekx[:, 1:, :] += np.conj(Rkx[:, 1:,:])
    # Subtract mean contribution
    Ekx[:,0,:] -= bar_pairs
    assert np.max(np.abs(np.imag(Ekx))) == 0
    # Convert to real and multiply by Fourier Transform to Fourier Series factor
    Ekx = np.real(Ekx) * Lx / (2*np.pi)

    # Compute spectra from Rkz by adding reflected negative wavenumbers
    Ekz            = Rkz[:, 0:(Nz/2+1), :].copy()
    Ekz[:, 1:, :] += Rkz[:, -1:-(Nz/2+1):-1,:]
    # Subtract mean contribution
    Ekz[:,0,:] -= bar_pairs
    # FIXME: find better criterion for tolerance of assert line below
    assert np.max(np.abs(np.imag(Ekz))) < 1E-10 #np.finfo(Ekz.dtype).eps
    # Convert to real and multiply by Fourier Transform to Fourier Series factor
    Ekz = np.real(Ekz) * Lz / (2*np.pi)

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
    Rz = Rz[:, :(Rz.shape[1]/2+1), :].real # Truncate small imag() part
    # Subtract mean contribution
    for i in xrange(Rx.shape[1]):
        Rx[:,i,:] -= bar_pairs
    for k in xrange(Rz.shape[1]):
        Rz[:,k,:] -= bar_pairs
    Rx = PhysicalData(y, np.mgrid[:Nx/2+1]*(Lx/Nx), *scalarpairs(Rx))
    Rz = PhysicalData(y, np.mgrid[:Nz/2+1]*(Lz/Nz), *scalarpairs(Rz))

    return ProcessResult(Ekx, Ekz, Rx, Rz, Rkx, Rkz, bar)


def load(h5filenames, verbose=False):
    """Load the data required by process() into a dict from named files.
    Averages of /twopoint_kx and /twopoint_kz are taken across all inputs.
    Other results reflect only the metadata from the last file loaded.
    """
    Rkx     = None
    Rkz     = None
    bar     = None
    bar_T   = None
    bar_u   = None
    bar_rho = None
    bar_cs  = None
    d       = {}
    D0      = None
    for ndx, h5filename in enumerate(h5filenames):
        if verbose:
            log.info("Loading %d/%d: %s"
                     % (ndx+1, len(h5filenames), h5filename))

        h5file = h5py.File(h5filename, 'r')
        Ns = 0
        sname = []

        # Grab number of collocation points and B-spline order
        Ny = h5file['Ny'][0]
        k  = h5file['k'][0]

        # Get "mass" matrix and convert to dense format
        D0T_gb = h5file['Dy0T'].value
        D0T    = gb.gb2ge(D0T_gb, Ny, k-2)
        D0     = D0T.transpose()

        if "antioch_constitutive_data" in h5file:
            Ns=h5file['antioch_constitutive_data'].attrs['Ns'][0]
            sname= np.chararray(Ns, itemsize=5)
            for s in xrange(0,Ns):
                sname[s]=h5file['antioch_constitutive_data'].attrs['Species_'+str(s)]
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
            Ns = Ns,
            sn = sname
        ))
        if Rkx is None:
            Rkx  = np.squeeze(h5file['twopoint_kx'][()].view(np.complex128))
        else:
            Rkx += np.squeeze(h5file['twopoint_kx'][()].view(np.complex128))
        if Rkz is None:
            Rkz  = np.squeeze(h5file['twopoint_kz'][()].view(np.complex128))
        else:
            Rkz += np.squeeze(h5file['twopoint_kz'][()].view(np.complex128))
        if bar_T is None:
            bar_T  = np.squeeze(h5file['bar_T'][()].view(np.float64))
        else:
            bar_T += np.squeeze(h5file['bar_T'][()].view(np.float64))
        if bar_u is None:
            bar_u  = np.squeeze(h5file['bar_u'][()].view(np.float64))
        else:
            bar_u += np.squeeze(h5file['bar_u'][()].view(np.float64))
        if bar_rho is None:
            bar_rho  = np.squeeze(h5file['bar_rho'][()].view(np.float64))
        else:
            bar_rho += np.squeeze(h5file['bar_rho'][()].view(np.float64))
        if "antioch_constitutive_data" in h5file:
            if Ns > 1:
                if bar_cs is None:
                    bar_cs  = np.squeeze(h5file['bar_cs_s'][()].view(np.float64))
                else:
                    bar_cs += np.squeeze(h5file['bar_cs_s'][()].view(np.float64))
            else:
                if bar_cs is None:
                    bar_cs  = 0 * np.array(bar_rho).reshape(1,Ny) + 1
                else:
                    bar_cs += 0 * np.array(bar_rho).reshape(1,Ny) + 1
        h5file.close()
    if Rkx is not None:
        Rkx /= len(h5filenames)
    if Rkz is not None:
        Rkz /= len(h5filenames)
    if bar_T is not None:
        bar_T /= len(h5filenames)
        bar_T  = np.dot(D0,bar_T)
    if bar_u is not None:
        bar_u /= len(h5filenames)
        bar_u  = np.transpose(np.dot(D0,np.transpose(bar_u)))
    if bar_rho is not None:
        bar_rho /= len(h5filenames)
        bar_rho = np.dot(D0,bar_rho)
    if bar_cs is not None:
        bar_cs /= len(h5filenames)
        bar_cs  = np.transpose(np.dot(D0,np.transpose(bar_cs)))

    # Pack mean fields, assume no variable remained as None,
    # except for possibly bar_cs
    if bar_cs is None:
        bar = np.concatenate((bar_T, bar_u, bar_rho        )).reshape(5+Ns, Ny)
    else:
        bar = np.concatenate((bar_T, bar_u, bar_rho, bar_cs)).reshape(5+Ns, Ny)

    d.update(dict(
        Rkx = Rkx,
        Rkz = Rkz,
        bar = bar
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
    inpickle   = None
    outpickle  = None
    outsuffix  = None
    verbose    = False
    yindex     = []
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "cho:p:vy:P:", ["help"])
        except getopt.error as msg:
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
            elif o == "-P":
                inpickle = a
            elif o == "-v":
                verbose = True
            elif o == "-y":
                yindex.append(int(a))
        if len(args) < 1 and not inpickle:
            print >>sys.stderr, "Too few arguments.  See --help."
            return 1
    except Usage as err:
        print >>sys.stderr, err.msg
        return 2

    # Push interactive mode off (in case we get used from IPython)
    was_interactive = plt.isinteractive()
    plt.interactive(False)

    # Load and process incoming data
    if inpickle:
        log.info("Loading processed results from pickle file")
        res = pickle.load(open(inpickle, "rb"))
    else:
        log.info("Loading data from restart files")
        data = load(args, verbose)
        log.info("Processing loaded data")
        res  = process(**data)

    # Unpack processed information into conveniently named variables
    (Ekx, Ekz, Rx, Rz, Rkx, Rkz, bar) = res

    # If requested, first pickle results to ease post mortem on crash
    if outpickle:
        log.info("Pickling data to %s" % (outpickle,))
        pickle.dump(res, open(outpickle, "wb"), -1)

    log.info("Preparing spectra at each requested y index")
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
        ax.set_xlim([min(Ekx.k[Ekx.k != 0].min(), Ekz.k[Ekz.k != 0].min()),
                     max(Ekx.k[Ekx.k != 0].max(), Ekz.k[Ekz.k != 0].max())]);
        ax.set_ylim([np.finfo(Ekx.k.dtype).eps, 1]);
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
        ax.loglog(Ekz.k, Ekz.rr[:,j], label="rho")
        ax.set_xlabel("Spanwise wavenumber")
        ax.set_title("Spectra at y = %g" % Ekx.y[j])
        ax.set_xlim([min(Ekx.k[Ekx.k != 0].min(), Ekz.k[Ekz.k != 0].min()),
                     max(Ekx.k[Ekx.k != 0].max(), Ekz.k[Ekz.k != 0].max())]);
        ax.set_ylim([np.finfo(Ekz.k.dtype).eps, 1]);
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.spectra.kz.'+outsuffix)
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
        ax.set_xlim([min(Ekx.k[Ekx.k != 0].min(), Ekz.k[Ekz.k != 0].min()),
                     max(Ekx.k[Ekx.k != 0].max(), Ekz.k[Ekz.k != 0].max())]);
        ax.set_ylim([np.finfo(Ekx.k.dtype).eps, 1]);
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
        ax.set_xlim([min(Ekx.k[Ekx.k != 0].min(), Ekz.k[Ekz.k != 0].min()),
                     max(Ekx.k[Ekx.k != 0].max(), Ekz.k[Ekz.k != 0].max())]);
        ax.set_ylim([np.finfo(Ekz.k.dtype).eps, 1]);
        ax.legend()
        if outsuffix:
            fig.savefig(str(j)+'.cospectra.kz.'+outsuffix)
            plt.close(fig)

    # If requested....
    if do_contour:
        log.info("Preparing contour plots for kx and kz")

        for index, name in enumerate(Ekx._fields):
            index += 2 # Hack based on knowledge of SpectralData
            (fig, ax, c1, c2, cbar) = contour(Ekx.k, Ekx.y, Ekx[index])
            ax.set_title("Spectra: " + name)
            ax.set_xlabel('Streamwise wavenumber')
            ax.set_ylabel('Wall-normal distance')
            if outsuffix:
                fig.savefig(name+'.kx.'+outsuffix)
                plt.close(fig)

        for index, name in enumerate(Ekz._fields):
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
