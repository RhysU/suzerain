/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#include <suzerain/diffwave.h>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/fftw.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/utility.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);       // Tear down BLAS

struct MPIFixture {
    MPIFixture()  { MPI_Init(NULL, NULL); }
    ~MPIFixture() { MPI_Finalize(); }
};
BOOST_GLOBAL_FIXTURE(MPIFixture);               // Setup and tear down MPI

// Performs an MPI_Barrier call at construction and destruction
struct MPIBarrierRAII {
    const MPI_Comm comm_;
    MPIBarrierRAII(MPI_Comm comm) : comm_(comm) { MPI_Barrier(comm_); }
    ~MPIBarrierRAII()                           { MPI_Barrier(comm_); }
};


// Test accumulate and apply in one shot to speed test execution,
// otherwise we have to setup and teardown P3DFFT twice
BOOST_AUTO_TEST_SUITE(accumulateAndApply)

using suzerain::pencil_grid;

// Succeed on small relative error small absolute error
static bool smallerr(const double expected,
                     const double actual,
                     const double relerr,
                     const double maxrelerr,
                     const double maxabserr)
{
    SUZERAIN_UNUSED(expected); // Included for predicate logging
    SUZERAIN_UNUSED(actual);   // Included for predicate logging

    return (relerr <= maxrelerr) || (std::abs(expected - actual) < maxabserr);
}

// Return value is worst observed relative error
static double test_accumulateAndApply_helper(const pencil_grid &pg,
                                             const int dxcnt,
                                             const int dzcnt,
                                             const double Lx,
                                             const double Ly,
                                             const double Lz,
                                             const int Nx,
                                             const int Nz,
                                             const double maxrelerr,
                                             const double maxabserr)
{
    // Note: Incoming pencil_grid describes dealiased extents
    const int nproc = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Create a uniform grid on [0, Lx) x [0, Ly) x [0, Lz)
    std::valarray<double> gridx(pg.local_physical_extent[0]);
    std::valarray<double> gridy(pg.local_physical_extent[1]);
    std::valarray<double> gridz(pg.local_physical_extent[2]);
    for (int i = 0; i < pg.local_physical_extent[0]; ++i) {
        gridx[i] = double(i+pg.local_physical_start[0]);
    }
    for (int j = 0; j < pg.local_physical_extent[1]; ++j) {
        gridy[j] = double(j+pg.local_physical_start[1]);
    }
    for (int k = 0; k < pg.local_physical_extent[2]; ++k) {
        gridz[k] = double(k+pg.local_physical_start[2]);
    }
    gridx *= Lx/pg.global_physical_extent[0];
    gridy *= Ly/pg.global_physical_extent[1];
    gridz *= Lz/pg.global_physical_extent[2];

    // Retrieve storage size, strides, and allocate working arrays
    const int nelem = pg.local_physical_storage();
    boost::shared_array<double>
        A(suzerain::blas::calloc_as<double>(nelem), suzerain::blas::free),
        B(suzerain::blas::calloc_as<double>(nelem), suzerain::blas::free);

    // Create composable test functions in the X and Z directions
    // Note that maximum wavenumber (inclusive) is (N-1)/2.
    periodic_function<> fx1(
            pg.global_physical_extent[0], Nx/2, M_PI/3, Lx, 11);
    periodic_function<> fx2(
            pg.global_physical_extent[0], Nx/2, M_PI/7, Lx, 13);
    periodic_function<> fz1(
            pg.global_physical_extent[2], Nz/2, M_PI/4, Lz, 17);
    periodic_function<> fz2(
            pg.global_physical_extent[2], Nz/2, M_PI/9, Lz, 19);

    // Track worst point-wise error magnitude for either test Gives an idea of
    // any lost precision and how loose any test tolerances may be in practice.
    double worstrelerr = 0;

    // P3DFFT does not normalize after transformations;
    // we need the rescaling factor handy.
    const double scale = pg.global_physical_extent[0]
                       * pg.global_physical_extent[2];

    // -----------------------------------------------------
    BOOST_TEST_MESSAGE("Test suzerain_diffwave_accumulate");
    // -----------------------------------------------------

    // Populate the synthetic fields
    std::fill(A.get(), A.get() + pg.local_physical_extent.prod(),
              std::numeric_limits<double>::quiet_NaN());
    std::fill(B.get(), B.get() + pg.local_physical_extent.prod(),
              std::numeric_limits<double>::quiet_NaN());
    {
        double *pA = A.get(), *pB = B.get();
        for (int j = 0; j < (int) gridy.size(); ++j) {
            for (int k = 0; k < (int) gridz.size(); ++k) {
                for (int i = 0; i < (int) gridx.size(); ++i) {
                    *pA++ = fx1.physical_evaluate(gridx[i])
                          * fz1.physical_evaluate(gridz[k])
                          * (j+1);
                    *pB++ = fx2.physical_evaluate(gridx[i])
                          * fz2.physical_evaluate(gridz[k])
                          * (j+1);
                }
            }
        }
    }

    // Transform to wave space
    pg.transform_physical_to_wave(A.get());
    pg.transform_physical_to_wave(B.get());

    // Differentiate-and-accumulate
    // Build FFT forward-and-inverse normalization factor into accumulation
    {
        const std::complex<double> alpha = 2.0/scale;
        const std::complex<double> *x    = (const std::complex<double> *) A.get();
        const std::complex<double> beta  = 3.0/scale;
              std::complex<double> *y    = (std::complex<double> *) B.get();
        const int Ny          = pg.global_physical_extent[1];
        const int dNx         = pg.global_physical_extent[0];
        const int dkbx        = pg.local_wave_start[0];
        const int dkex        = pg.local_wave_end[0];
        const int dNz         = pg.global_physical_extent[2];
        const int dkbz        = pg.local_wave_start[2];
        const int dkez        = pg.local_wave_end[2];

        suzerain_diffwave_accumulate(
            dxcnt, dzcnt, alpha, x, beta, y, Lx, Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
    }

    // Transform to physical space
    pg.transform_wave_to_physical(A.get());
    pg.transform_wave_to_physical(B.get());

    // Ensure the synthetic accumulated field came back cleanly
    {
        const double *pB = B.get();
        for (int j = 0; j < (int) gridy.size(); ++j) {
            for (int k = 0; k < (int) gridz.size(); ++k) {
                for (int i = 0; i < (int) gridx.size(); ++i) {
                    {
                        const double cfx1
                            = fx1.physical_evaluate(gridx[i], dxcnt);
                        const double cfz1
                            = fz1.physical_evaluate(gridz[k], dzcnt);
                        const double cfx2
                            = fx2.physical_evaluate(gridx[i]);
                        const double cfz2
                            = fz2.physical_evaluate(gridz[k]);
                        const double expected_B
                            = 2*cfx1*cfz1*(j+1) + 3*cfx2*cfz2*(j+1);

                        double relerr = relative_error(*pB, expected_B);
                        BOOST_CHECK_PREDICATE(smallerr,
                                (expected_B)(*pB)(relerr)(maxrelerr)(maxabserr));
                        worstrelerr = std::max(worstrelerr, relerr);

                        ++pB;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------
    BOOST_TEST_MESSAGE("Test suzerain_diffwave_apply");
    // -----------------------------------------------------

    // Populate the synthetic fields
    std::fill(A.get(), A.get() + pg.local_physical_extent.prod(),
              std::numeric_limits<double>::quiet_NaN());
    {
        double *pA = A.get();
        for (int j = 0; j < (int) gridy.size(); ++j) {
            for (int k = 0; k < (int) gridz.size(); ++k) {
                for (int i = 0; i < (int) gridx.size(); ++i) {
                    *pA++ = fx1.physical_evaluate(gridx[i])
                          * fz1.physical_evaluate(gridz[k])
                          * (j+1);
                }
            }
        }
    }

    // Dump physical space contents on one processor for dxcnt == dzcnt == 0
    if (nproc == 1 && dxcnt == 0 && dzcnt == 0) {
        for (int i = 0; i < pg.local_physical_extent.prod(); ++i) {
            BOOST_TEST_MESSAGE("DEBUG: Physical space data linear index " << i
                               << ": " << A[i]);
        }
    }

    // Transform to wave space
    pg.transform_physical_to_wave(A.get());

    // Dump wave space contents on one processor for dxcnt == dzcnt == 0
    if (nproc == 1 && dxcnt == 0 && dzcnt == 0) {
        for (int i = 0; i < pg.local_wave_extent.prod(); ++i) {
            BOOST_TEST_MESSAGE("DEBUG: Wave space data at linear index " << i
                               << ": " << A[2*i] << "+" << A[2*i+1] << "i");
        }
    }

    // Differentiate-and-accumulate
    // Build FFT forward-and-inverse normalization factor into accumulation
    {
        const std::complex<double> alpha = 2.0/scale;
              std::complex<double> *x    = (std::complex<double> *) A.get();
        const int Ny   = pg.global_physical_extent[1];
        const int dNx  = pg.global_physical_extent[0];
        const int dkbx = pg.local_wave_start[0];
        const int dkex = pg.local_wave_end[0];
        const int dNz  = pg.global_physical_extent[2];
        const int dkbz = pg.local_wave_start[2];
        const int dkez = pg.local_wave_end[2];

        suzerain_diffwave_apply(
            dxcnt, dzcnt, alpha, x, Lx, Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
    }

    // Transform to physical space
    pg.transform_wave_to_physical(A.get());

    // Ensure the synthetic accumulated field came back cleanly
    {
        const double *pA = A.get();
        for (int j = 0; j < (int) gridy.size(); ++j) {
            for (int k = 0; k < (int) gridz.size(); ++k) {
                for (int i = 0; i < (int) gridx.size(); ++i) {
                    {
                        const double cfx1
                            = fx1.physical_evaluate(gridx[i], dxcnt);
                        const double cfz1
                            = fz1.physical_evaluate(gridz[k], dzcnt);
                        const double expected_A
                            = 2*cfx1*cfz1*(j+1);

                        double relerr = relative_error(*pA, expected_A);
                        BOOST_CHECK_PREDICATE(smallerr,
                                (expected_A)(*pA)(relerr)(maxrelerr)(maxabserr));
                        worstrelerr = std::max(worstrelerr, relerr);

                        ++pA;
                    }
                }
            }
        }
    }

    return worstrelerr;
}

static void test_accumulateAndApply(const int Ny,
                                    const int Nx,
                                    const int dNx,
                                    const int Nz,
                                    const int dNz)
{
    const int MAX_DXCNT_INCLUSIVE = 4;
    const int MAX_DZCNT_INCLUSIVE = 4;

    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    const int np     = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    suzerain::array<int,3> global_physical_extents = {{ dNx, Ny, dNz }};
    suzerain::array<int,2> processor_grid = {{ 0, 0 }};

    if (!procid) {
        BOOST_TEST_MESSAGE("Establishing pencil_grid for extents "
                           << global_physical_extents);
    }
#if defined SUZERAIN_HAVE_P3DFFT && defined SUZERAIN_HAVE_UNDERLING
    {
        // Perform an implementation-to-implementation shootout
        suzerain::pencil_grid_p3dfft p3(
                global_physical_extents, processor_grid,
                suzerain::fftw::estimate, suzerain::fftw::estimate);
        suzerain::pencil_grid_underling un(
                global_physical_extents, processor_grid,
                suzerain::fftw::estimate, suzerain::fftw::estimate);

#define EIGEN_CHECK_EQUAL(a,b) \
        BOOST_CHECK_EQUAL_COLLECTIONS(a.data(), a.data() + a.size(), \
                                      b.data(), b.data() + b.size())

        // Global extents should always match
        EIGEN_CHECK_EQUAL(p3.global_physical_extent,
                          un.global_physical_extent);
        EIGEN_CHECK_EQUAL(p3.global_wave_extent,
                          un.global_wave_extent);

        if (np == 1) {
            // On a single processor all details should match
            EIGEN_CHECK_EQUAL(p3.local_physical_start,
                              un.local_physical_start);
            EIGEN_CHECK_EQUAL(p3.local_physical_end,
                              un.local_physical_end);
            EIGEN_CHECK_EQUAL(p3.local_physical_extent,
                              un.local_physical_extent);
            EIGEN_CHECK_EQUAL(p3.local_wave_start,
                              un.local_wave_start);
            EIGEN_CHECK_EQUAL(p3.local_wave_end,
                              un.local_wave_end);
            EIGEN_CHECK_EQUAL(p3.local_wave_extent,
                              un.local_wave_extent);

            // Underling should require slightly bigger buffers
            BOOST_CHECK_LE(p3.local_physical_storage(),
                           un.local_physical_storage());
            BOOST_CHECK_LE(p3.local_wave_storage(),
                           un.local_wave_storage());
        }

#undef EIGEN_CHECK_EQUAL

    }
#endif

    suzerain::pencil_grid_default pg(
            global_physical_extents, processor_grid,
            suzerain::fftw::estimate, suzerain::fftw::estimate);

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            const double maxrelerr = 3.0*std::pow(10, -11 + (dxcnt+dzcnt)/3.0);
            const double maxabserr = 1e-12;
            if (!procid) {
                BOOST_TEST_MESSAGE("Testing"
                                    << " Ny=" << Ny
                                    << ", Nx=" << Nx
                                    << ", dNx=" << dNx
                                    << ", Nz=" << Nz
                                    << ", dNz=" << dNz
                                    << " for dx=" << dxcnt
                                    << ", dz=" << dzcnt
                                    << ", np=" << np
                                    << ", maxrelerr=" << maxrelerr);
            }
            MPIBarrierRAII barrier(MPI_COMM_WORLD);
            const double worstrelerr = test_accumulateAndApply_helper(
                pg, dxcnt, dzcnt, 4*M_PI, 2, 4*M_PI/3, Nx, Nz,
                maxrelerr, maxabserr);
            BOOST_TEST_MESSAGE("Rank " << procid
                               << " observed local maximum relative error = "
                               << std::scientific
                               << std::setw(16)
                               << worstrelerr);
        }
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_1D_not_dealiased )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*                       Ny,  Nx, dNx,  Nz, dNz */
        test_accumulateAndApply(  1,   1,   1,   8,   8);
        test_accumulateAndApply(  1,   1,   1,   7,   7);
        test_accumulateAndApply(  1,   8,   8,   1,   1);
        test_accumulateAndApply(  1,   7,  12,   1,   1);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_1D_dealiased_z )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*                       Ny,  Nx, dNx,  Nz, dNz */
        test_accumulateAndApply(  1,   1,   1,   6,   8);
        test_accumulateAndApply(  1,   1,   1,   6,   9);
        test_accumulateAndApply(  1,   1,   1,   5,   8);
        test_accumulateAndApply(  1,   1,   1,   5,   9);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_1D_dealiased_x )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*                       Ny,  Nx, dNx,  Nz, dNz */
        test_accumulateAndApply(  1,   6,   8,   1,   1);
        test_accumulateAndApply(  1,   6,  10,   1,   1);
        test_accumulateAndApply(  1,   5,   8,   1,   1);
        test_accumulateAndApply(  1,   5,  10,   1,   1);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_2D_not_dealiased )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*                       Ny,  Nx, dNx,  Nz, dNz */
        test_accumulateAndApply(  1,   4,   4,   4,   4);
        test_accumulateAndApply(  1,   8,   8,   8,   8);
        test_accumulateAndApply(  1,   8,   8,   7,   7);
        test_accumulateAndApply(  1,   8,   8,   8,   8);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_not_dealiased )
{
    /*                       Ny,  Nx, dNx,  Nz, dNz */
    test_accumulateAndApply(  3,  24,  24,  40,  40 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_z )
{
    /*                       Ny,  Nx, dNx,  Nz, dNz */
    test_accumulateAndApply(  2,  24,  24,  40,  60 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_x )
{
    /*                       Ny,  Nx, dNx,  Nz, dNz */
    test_accumulateAndApply(  2,  24,  36,  40,  40 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_xz )
{
    /*                       Ny,  Nx, dNx,  Nz, dNz */
    test_accumulateAndApply(  2,  24,  36,  40,  40 );
}

BOOST_AUTO_TEST_SUITE_END()
