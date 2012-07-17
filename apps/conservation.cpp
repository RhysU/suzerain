//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// conservation.cpp: compute discrete conservation error for given basis
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <esio/esio.h>
#include <suzerain/bspline.hpp>
#include <suzerain/math.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/program_options.hpp>
#include <suzerain/utility.hpp>

#include "logging.hpp"
#include "precision.hpp"
#include "support.hpp"

#pragma warning(disable:383 1572)

// Global B-spline related-details initialized in main()
static boost::shared_ptr<suzerain::bspline>      b;
static boost::shared_ptr<suzerain::bsplineop>    bop;    // Collocation
static boost::shared_ptr<suzerain::bsplineop>    gop;    // Galerkin L2
static boost::shared_ptr<suzerain::bsplineop_lu> boplu;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    logging::initialize(MPI_COMM_WORLD,             // Initialize logging
                        channel::log4cxx_config);
    WARN0_ENABLE();                                 // Disable chattiness

    // Process incoming program arguments from command line, input files
    double L       = 2;
    int    N       = 16;
    int    k       = 5;
    {
        suzerain::ProgramOptions options(
                "Check conservative nature of discrete B-spline operators");
        using boost::program_options::value;
        options.add_options()
            ("L", value(&L)->default_value(L),
             "Domain size")
            ("N", value(&N)->default_value(N),
             "Number of degrees of freedom in the basis")
            ("k", value(&k)->default_value(k),
             "Wall-normal B-spline order (4 indicates piecewise cubics)")
        ;
        std::vector<std::string> positional = options.process(argc, argv);

        if (positional.size() != 0) {
            FATAL("No non-option arguments permitted");
            return EXIT_FAILURE;
        }
        if (L <= 0) {
            FATAL("L <= 0 not permitted");
            return EXIT_FAILURE;
        }
        if (N <= 0) {
            FATAL("N <= 0 not permitted");
            return EXIT_FAILURE;
        }
        if (k < 3 /* quadratics */) {
            FATAL("k >= 3 required for one non-trivial spatial derivatives");
            return EXIT_FAILURE;
        }
    }

    // Modify IEEE settings after startup complete as startup relies on NaNs
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    Eigen::VectorXr  exp(N);
    exp.setZero(N, 1);
    exp[0] = -1;
    exp[N-1] = 1;

    Eigen::VectorXr  vec(N);
    Eigen::MatrixXXr mat(N,N);
    for (real_t htdelta = 0; htdelta < 7; htdelta += 0.1) {
        channel::create(N, k, 0.0, L, htdelta, b, bop);
        b->integration_coefficients(0, vec.data());
        boplu = boost::make_shared<suzerain::bsplineop_lu>(*bop);
        boplu->factor_mass(*bop);

        mat = Eigen::MatrixXXr::Identity(N,N);
        boplu->solve(N, mat.data(), 1, N);
        bop->apply(1, N, 1.0, mat.data(), 1, N);
        boplu->solve(N, mat.data(), 1, N);
        vec = vec.transpose() * mat;

        std::cout << htdelta << "\t"
                  << (vec - exp).norm() << "\t"
                  << vec.transpose() << std::endl;
    }

    return EXIT_SUCCESS;
}
