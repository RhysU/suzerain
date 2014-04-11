#!/usr/bin/env python
import collections
import h5py
import numpy as np
import numpy.fft as fft
import sys

# Load scenario parameters from restart file
s = h5py.File('restart00000.h5', 'r')
y = np.squeeze(s['collocation_points_y'][()])
Nx, Ny, Nz = s['Nx'][0], s['Ny'][0], s['Nz'][0]
Lx, Ly, Lz = s['Lx'][0], s['Ly'][0], s['Lz'][0]
kx, kz     = s['kx'], s['kz']

# Load two point versus (f, kx, yj) and (f, kx, yj) and transform to physical
Rkx = np.squeeze(s['twopoint_kx'][()].view(np.complex128))
Rkz = np.squeeze(s['twopoint_kz'][()].view(np.complex128))
Rx  = fft.irfft(Nx * Rkx, axis=1)  # Nonnormalized iFFT
Rz  = fft.ifft (Nz * Rkz, axis=1)  # Nonnormalized iFFT

# Compute spectra from Rkx using conjugate-symmetry of Rkx
Ekx          = Rkx.copy()
Ekx[:,1:,:] += np.conj(Rkx[:,1:,:])
assert np.max(np.abs(np.imag(Ekx))) == 0
Ekx = np.real(Ekx)

# Compute spectra from Rkz by reflecting and adding the negative wavenumbers
Ekz          = Rkz[:,0:(Nz/2),:].copy()
Ekz[:,1:,:] += Rkz[:,-1:-(Nz/2):-1,:]
assert np.max(np.abs(np.imag(Ekz))) < np.finfo(Ekz.dtype).eps
Ekz = np.real(Ekz)

# Prepare named view of each spectral scalar pair like 'Ex.uv' tracking y,k
SpectralData = collections.namedtuple('SpectralData',
                                      ['y', 'k', 'shape',
                                       'TT', 'Tu', 'Tv', 'Tw', 'Tr',
                                             'uu', 'uv', 'uw', 'ur',
                                                   'vv', 'vw', 'vr',
                                                         'ww', 'wr',
                                                               'rr' ])
def scalarpairs(d):
    return map(lambda i: np.squeeze(d[i,:,:]),
               np.arange(d.shape[0]))
Rkx = SpectralData(y, kx,          Rkx.shape[1:], *scalarpairs(Rkx))
Rkz = SpectralData(y, kz,          Rkz.shape[1:], *scalarpairs(Rkz))
Ekx = SpectralData(y, kx,          Ekx.shape[1:], *scalarpairs(Ekx))
Ekz = SpectralData(y, kz[:Nz/2+1], Ekx.shape[1:], *scalarpairs(Ekz))

# Prepare named view of each physical scalar pair like 'Rx.uv' tracking y,x
PhysicalData = collections.namedtuple('PhysicalData',
                                      ['y', 'x', 'shape',
                                       'TT', 'Tu', 'Tv', 'Tw', 'Tr',
                                             'uu', 'uv', 'uw', 'ur',
                                                   'vv', 'vw', 'vr',
                                                         'ww', 'wr',
                                                               'rr' ])
# Only half the two-point is interesting as they are periodic and even
Rx = Rx[:,:(Rx.shape[1]/2+1),:]
Rz = Rz[:,:(Rz.shape[1]/2+1),:]
Rx = PhysicalData(y, np.mgrid[:Nx/2+1]*(Lx/Nx), Rx.shape[1:], *scalarpairs(Rx))
Rz = PhysicalData(y, np.mgrid[:Nz/2+1]*(Lz/Nz), Rz.shape[1:], *scalarpairs(Rz))
