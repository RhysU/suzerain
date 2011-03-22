#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/diffwave.h>
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

static bool smallrelerror(const double expected,
                          const double actual,
                          const double relerror,
                          const double maxrelerror)
{
    SUZERAIN_UNUSED(expected); // Included for predicate logging
    SUZERAIN_UNUSED(actual);   // Included for predicate logging
    return relerror <= maxrelerror;
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
                                             const double maxrelerror)
{
    // Note: Incoming pencil_grid describes dealiased extents

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
    boost::shared_array<double> A(new double[nelem]), B(new double[nelem]);

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
    double worstrelerror = 0;

    // P3DFFT does not normalize after transformations;
    // we need the rescaling factor handy.
    const double scale = pg.global_physical_extent[0]
                      * pg.global_physical_extent[2];

    // ---------------------------------
    // Test suzerain_diffwave_accumulate
    // ---------------------------------

    // Populate the synthetic fields
    std::fill(A.get(),A.get()+nelem,std::numeric_limits<double>::quiet_NaN());
    std::fill(B.get(),B.get()+nelem,std::numeric_limits<double>::quiet_NaN());
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
        const double alpha[2] = { 2.0/scale, 0.0/scale };
        const double (*x)[2]  = reinterpret_cast<const double (*)[2]>(A.get());
        const double beta[2]  = { 3.0/scale, 0.0/scale };
        double (*y)[2]        = reinterpret_cast<double (*)[2]>(B.get());
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

                        double relerror = relative_error(*pB, expected_B);
                        BOOST_CHECK_PREDICATE(smallrelerror,
                                (expected_B)(*pB)(relerror)(maxrelerror));
                        worstrelerror = std::max(worstrelerror, relerror);

                        ++pB;
                    }
                }
            }
        }
    }

    // ---------------------------------
    // Test suzerain_diffwave_apply
    // ---------------------------------

    // Populate the synthetic fields
    std::fill(A.get(),A.get()+nelem,std::numeric_limits<double>::quiet_NaN());
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

    // Transform to wave space
    pg.transform_physical_to_wave(A.get());

    // Differentiate-and-accumulate
    // Build FFT forward-and-inverse normalization factor into accumulation
    {
        const double alpha[2] = { 2.0/scale, 0.0/scale };
        double (*x)[2]        = reinterpret_cast<double (*)[2]>(A.get());
        const int Ny          = pg.global_physical_extent[1];
        const int dNx         = pg.global_physical_extent[0];
        const int dkbx        = pg.local_wave_start[0];
        const int dkex        = pg.local_wave_end[0];
        const int dNz         = pg.global_physical_extent[2];
        const int dkbz        = pg.local_wave_start[2];
        const int dkez        = pg.local_wave_end[2];

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

                        double relerror = relative_error(*pA, expected_A);
                        BOOST_CHECK_PREDICATE(smallrelerror,
                                (expected_A)(*pA)(relerror)(maxrelerror));
                        worstrelerror = std::max(worstrelerror, relerror);

                        ++pA;
                    }
                }
            }
        }
    }

    return worstrelerror;
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

    boost::array<int,3> global_physical_extents = {{ dNx, Ny, dNz }};
    boost::array<int,2> processor_grid = {{ 0, 0 }};
    pencil_grid pg(global_physical_extents, processor_grid);

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            const double maxrelerror
                = 3.0*std::pow(10, -11 + (dxcnt+dzcnt)/3.0);
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
                                    << ", maxrelerror=" << maxrelerror);
            }
            MPIBarrierRAII barrier(MPI_COMM_WORLD);
            const double worstrelerror = test_accumulateAndApply_helper(
                pg, dxcnt, dzcnt, 4*M_PI, 2, 4*M_PI/3, Nx, Nz, maxrelerror);
            BOOST_TEST_MESSAGE("Rank " << procid
                               << " observed local maximum relative error = "
                               << std::scientific
                               << std::setw(16)
                               << worstrelerror);
        }
    }
}

BOOST_AUTO_TEST_CASE( p3dfft_sanity_check )
{
    BOOST_REQUIRE(p3dfft_using_stride1());
    BOOST_REQUIRE_EQUAL(p3dfft_get_precision(), 2);
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
        test_accumulateAndApply(  1,   7,   7,   1,   1);
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
        test_accumulateAndApply(  1,   6,   9,   1,   1);
        test_accumulateAndApply(  1,   5,   8,   1,   1);
        test_accumulateAndApply(  1,   5,   9,   1,   1);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_2D_not_dealiased )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*                       Ny,  Nx, dNx,  Nz, dNz */
        test_accumulateAndApply(  1,   8,   8,   8,   8);
        test_accumulateAndApply(  1,   7,   7,   7,   7);
        test_accumulateAndApply(  1,   7,   7,   8,   8);
        test_accumulateAndApply(  1,   8,   8,   7,   7);
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
