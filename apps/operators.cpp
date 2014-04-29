//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

/** @file
 * Compute discrete conservation error for given basis
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/math.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/utility.hpp>

#pragma warning(disable:383 1572)

using namespace suzerain;

// Global B-spline related-details initialized in main()
static shared_ptr<bspline>      b;
static shared_ptr<bsplineop>    cop;    // Collocation
static shared_ptr<bsplineop>    gop;    // Galerkin L2
static shared_ptr<bsplineop_lu> boplu;

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    support::logging::initialize(MPI_COMM_WORLD,    // Initialize logging
                               support::log4cxx_config);
    WARN0_ENABLE();                                 // Disable chattiness

    // Process incoming program arguments from command line, input files
    double L       = 2;
    int    N       = 16;
    int    k       = 5;
    {
        support::program_options options(
                "Check discrete B-spline operator quality against sin(y+1)");
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
        if (k < 4 /* cubics */) {
            FATAL("k >= 4 required for one non-trivial spatial derivatives");
            return EXIT_FAILURE;
        }
    }

    // Modify IEEE settings after startup complete as startup relies on NaNs
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(mpi::comm_rank(MPI_COMM_WORLD));

    // Over the stretching parameters of interest...
    for (real_t htdelta = -3.5; htdelta < 3.5; htdelta += 0.1) {

        // ...generate the desired bases
        support::create_bsplines(N, k, 0.0, L, htdelta, b, cop);
        boplu = make_shared<bsplineop_lu>(*cop);
        boplu->factor_mass(*cop);

        // ...obtain coefficient representation of sin(y+1)
        VectorXr coeff(N);
        for (int j = 0; j < N; ++j) {
            coeff[j] = std::sin(b->collocation_point(j) + 1);
        }
        boplu->solve(1, coeff.data(), 1, N);

        // ...compute the second derivative in place
        cop->apply(2, 1, 1.0, coeff.data(), 1, N);

        // ...find the absolute error relative to analytical -sin(y+1)
        for (int j = 0; j < N; ++j) {
            coeff[j] += std::sin(b->collocation_point(j) + 1);
            coeff[j] = std::abs(coeff[j]);
        }

        // ...and output the maximum and sum of absolute errors
        std::cout << htdelta << "\t"
                  << coeff.maxCoeff() << "\t"
                  << coeff.sum() << std::endl;
    }

    return EXIT_SUCCESS;
}
