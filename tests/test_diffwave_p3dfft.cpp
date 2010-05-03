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


BOOST_AUTO_TEST_SUITE(accumulate)

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
static double test_accumulate_helper(const pencil_grid &pg,
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
    typedef pencil_grid::index index;
    typedef pencil_grid::size_type size_type;

    // Create a uniform grid on [0, Lx) x [0, Ly) x [0, Lz)
    std::valarray<double> gridx(pg.local_physical_extent()[0]);
    std::valarray<double> gridy(pg.local_physical_extent()[1]);
    std::valarray<double> gridz(pg.local_physical_extent()[2]);
    for (index i = 0; i < pg.local_physical_extent()[0]; ++i) {
        gridx[i] = double(i+pg.local_physical_start()[0]);
    }
    for (index j = 0; j < pg.local_physical_extent()[1]; ++j) {
        gridy[j] = double(j+pg.local_physical_start()[1]);
    }
    for (index k = 0; k < pg.local_physical_extent()[2]; ++k) {
        gridz[k] = double(k+pg.local_physical_start()[2]);
    }
    gridx *= Lx/pg.global_extents()[0];
    gridy *= Ly/pg.global_extents()[1];
    gridz *= Lz/pg.global_extents()[2];

    // Retrieve storage size, strides, and allocate working arrays
    const int nelem = pg.local_physical_storage();
    boost::shared_array<double> A(new double[nelem]), B(new double[nelem]);
    std::fill(A.get(),A.get()+nelem,std::numeric_limits<double>::quiet_NaN());
    std::fill(B.get(),B.get()+nelem,std::numeric_limits<double>::quiet_NaN());

    // Create composable test functions in the X and Z directions
    // Note that maximum wavenumber (inclusive) is (N-1)/2.
    periodic_function<> fx1(pg.global_extents()[0], Nx/2, M_PI/3, Lx, 11);
    periodic_function<> fx2(pg.global_extents()[0], Nx/2, M_PI/7, Lx, 13);
    periodic_function<> fz1(pg.global_extents()[2], Nz/2, M_PI/4, Lz, 17);
    periodic_function<> fz2(pg.global_extents()[2], Nz/2, M_PI/9, Lz, 19);

    // Populate the synthetic fields
    {
        double *pA = A.get(), *pB = B.get();
        for (index j = 0; j < gridy.size(); ++j) {
            for (index k = 0; k < gridz.size(); ++k) {
                for (index i = 0; i < gridx.size(); ++i) {
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
    const size_type scale = pg.global_extents()[0] *pg.global_extents()[2];
    {
        const double alpha[2] = { 2.0/scale, 0.0/scale };
        const double (*x)[2]  = reinterpret_cast<const double (*)[2]>(A.get());
        const double beta[2]  = { 3.0/scale, 0.0/scale };
        double (*y)[2]        = reinterpret_cast<double (*)[2]>(B.get());
        const int Ny          = pg.global_extents()[1];
        const int dNx         = pg.global_extents()[0];
        const int dkbx        = pg.local_wave_start()[0];
        const int dkex        = pg.local_wave_end()[0];
        const int dNz         = pg.global_extents()[2];
        const int dkbz        = pg.local_wave_start()[2];
        const int dkez        = pg.local_wave_end()[2];

        suzerain_diffwave_accumulate(
            dxcnt, dzcnt, alpha, x, beta, y, Lx, Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
    }

    // Transform to physical space
    pg.transform_wave_to_physical(A.get());
    pg.transform_wave_to_physical(B.get());

    // Track worst point-wise error magnitude
    double worstrelerror = 0;

    // Ensure the synthetic accumulated field came back cleanly
    {
        const double *pB = B.get();
        for (index j = 0; j < gridy.size(); ++j) {
            for (index k = 0; k < gridz.size(); ++k) {
                for (index i = 0; i < gridx.size(); ++i) {
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

    return worstrelerror;
}

static void test_accumulate(const int Ny,
                            const int Nx,
                            const int dNx,
                            const int Nz,
                            const int dNz)
{
    const int MAX_DXCNT_INCLUSIVE = 4;
    const int MAX_DZCNT_INCLUSIVE = 4;

    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    const int np     = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    pencil_grid::size_type_3d global_extents = { dNx, Ny, dNz };
    pencil_grid::size_type_2d processor_grid = { 0, 0 };
    pencil_grid pg(global_extents, processor_grid);

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            const double maxrelerror = 5*std::pow(10, -11 + (dxcnt+dzcnt)/4.0);
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
            const double worstrelerror = test_accumulate_helper(
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
        /*               Ny,  Nx, dNx,  Nz, dNz */
        test_accumulate(  1,   1,   1,   8,   8);
        test_accumulate(  1,   1,   1,   7,   7);
        test_accumulate(  1,   8,   8,   1,   1);
        test_accumulate(  1,   7,   7,   1,   1);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_1D_dealiased_z )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*               Ny,  Nx, dNx,  Nz, dNz */
        test_accumulate(  1,   1,   1,   6,   8);
        test_accumulate(  1,   1,   1,   6,   9);
        test_accumulate(  1,   1,   1,   5,   8);
        test_accumulate(  1,   1,   1,   5,   9);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_1D_dealiased_x )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*               Ny,  Nx, dNx,  Nz, dNz */
        test_accumulate(  1,   6,   8,   1,   1);
        test_accumulate(  1,   6,   9,   1,   1);
        test_accumulate(  1,   5,   8,   1,   1);
        test_accumulate(  1,   5,   9,   1,   1);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_quasi_2D_not_dealiased )
{
    // Check small grids with odd sizes only in the single rank case.
    // P3DFFT cannot distribute them on non-1x1 Cartesian processor grids.
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) == 1) {
        /*               Ny,  Nx, dNx,  Nz, dNz */
        test_accumulate(  1,   8,   8,   8,   8);
        test_accumulate(  1,   7,   7,   7,   7);
        test_accumulate(  1,   7,   7,   8,   8);
        test_accumulate(  1,   8,   8,   7,   7);
    }
}

BOOST_AUTO_TEST_CASE( accumulate_not_dealiased )
{
    /*               Ny,  Nx, dNx,  Nz, dNz */
    test_accumulate(  3,  24,  24,  40,  40 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_z )
{
    /*               Ny,  Nx, dNx,  Nz, dNz */
    test_accumulate(  2,  24,  24,  40,  60 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_x )
{
    /*               Ny,  Nx, dNx,  Nz, dNz */
    test_accumulate(  2,  24,  36,  40,  40 );
}

BOOST_AUTO_TEST_CASE( accumulate_dealiased_xz )
{
    /*               Ny,  Nx, dNx,  Nz, dNz */
    test_accumulate(  2,  24,  36,  40,  40 );
}

BOOST_AUTO_TEST_SUITE_END()
