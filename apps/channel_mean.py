#!/usr/bin/env python
import cvxopt
import cvxopt.blas
import matplotlib.pyplot as plt
import numpy
import sys
import tables

def process(filename):

    # Open restart file
    h = tables.openFile(filename, "r")

    # Obtain real-valued coefficients on (0,0)-mode pencil in Y
    c_rho  = cvxopt.matrix(numpy.array(h.root.rho.read() [0,0,:,0]))
    c_rhou = cvxopt.matrix(numpy.array(h.root.rhou.read()[0,0,:,0]))
    c_rhov = cvxopt.matrix(numpy.array(h.root.rhov.read()[0,0,:,0]))
    c_rhow = cvxopt.matrix(numpy.array(h.root.rhow.read()[0,0,:,0]))
    c_rhoe = cvxopt.matrix(numpy.array(h.root.rhoe.read()[0,0,:,0]))

    # Obtain banded operator to compute point values at collocation points
    x   = h.root.collocation_points.read()
    m   = h.root.Dy0.attrs.m[0]
    n   = h.root.Dy0.attrs.n[0]
    kl  = h.root.Dy0.attrs.kl[0]
    ku  = h.root.Dy0.attrs.ku[0]
    Dy0 = cvxopt.matrix(h.root.Dy0.read().transpose())

    # Compute values at collocation points
    rho = cvxopt.matrix(0.0, c_rho.size)
    cvxopt.blas.gbmv(Dy0, m, kl, c_rho, rho)

    rhou = cvxopt.matrix(0.0, c_rhou.size)
    cvxopt.blas.gbmv(Dy0, m, kl, c_rhou, rhou)

    rhov = cvxopt.matrix(0.0, c_rhov.size)
    cvxopt.blas.gbmv(Dy0, m, kl, c_rhov, rhov)

    rhow = cvxopt.matrix(0.0, c_rhow.size)
    cvxopt.blas.gbmv(Dy0, m, kl, c_rhow, rhow)

    rhoe = cvxopt.matrix(0.0, c_rhoe.size)
    cvxopt.blas.gbmv(Dy0, m, kl, c_rhoe, rhoe)

    # Plot several mean quantities
    plt.plot(x, rho, label='rho')
    plt.plot(x, rhou, label='rhou')
    plt.plot(x, rhoe, label='rhoe')
    plt.legend()
    plt.show()

    # Close restart file
    h.close()

if __name__ == '__main__':
    for filename in sys.argv[1:]:
        print 'Processing ' + filename
        process(filename)
