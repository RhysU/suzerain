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

static void test_accumulate_helper(const pencil_grid &pg,
                                   const int dxcnt,
                                   const int dzcnt,
                                   const double Lx,
                                   const double Ly,
                                   const double Lz,
                                   const int Nx,
                                   const int Nz)
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
                          + fz1.physical_evaluate(gridz[k]);
                    *pB++ = fx2.physical_evaluate(gridx[i])
                          + fz2.physical_evaluate(gridz[k]);
                }
            }
        }
    }

    // Transform to wave space
    pg.transform_physical_to_wave(A.get());
    pg.transform_physical_to_wave(B.get());

    // Differentiate-and-accumulate
    {
        const double alpha[2] = { 2.0, 0.0 };
        const double (*x)[2]  = reinterpret_cast<const double (*)[2]>(A.get());
        const double beta[2]  = { 3.0, 0.0 };
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

    // Ensure the synthetic fields came back cleanly
    const size_type scale = pg.global_extents()[0] *pg.global_extents()[2];
    const double    close = std::pow(10, -9 + (dxcnt+dzcnt)/2.5)
                          * std::sqrt((double) scale)
                          * 5;
    {
        double *pA = A.get(), *pB = B.get();
        for (index j = 0; j < gridy.size(); ++j) {
            for (index k = 0; k < gridz.size(); ++k) {
                for (index i = 0; i < gridx.size(); ++i) {
                    {
                        const double cfx1 = fx1.physical_evaluate(gridx[i]);
                        const double cfz1 = fz1.physical_evaluate(gridz[k]);
                        const double expected_A = scale*(cfx1 + cfz1);
                        BOOST_CHECK_CLOSE(*pA++, expected_A, close);
                    }

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
                            = scale*(2*cfx1 + 2*cfz1 + 3*cfx2 + 3*cfz2);
                        BOOST_CHECK_CLOSE(*pB++, expected_B, close);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( accumulate )
{
    const int MAX_DXCNT_INCLUSIVE = 4;
    const int MAX_DZCNT_INCLUSIVE = 4;

    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    boost::array<int,5> c[] = {
        /* Ny,  Nx, dNx,  Nz, dNz */
        {   1,   1,   1,   6,   6 }

//      {   4,  24,  24,  40,  40 }
//     ,{   4,  24,  24,  40,  60 } // Dealiased Z
//     ,{   4,  24,  36,  40,  40 } // Dealiased X
//     ,{   4,  24,  36,  40,  60 } // Dealiased X,Z
    };

    for (int l = 0; l < sizeof(c)/sizeof(c[0]); ++l) {
        // global_extents use dealiased dimensions
        pencil_grid::size_type_3d global_extents = { c[l][2],c[l][0],c[l][4] };
        pencil_grid::size_type_2d processor_grid = { 0, 0 };
        pencil_grid pg(global_extents, processor_grid);
        assert(pg.local_wave_extent()[1] == c[l][0]); // P3DFFT using STRIDE1?

        for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
            for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
                if (!procid) {
                    BOOST_TEST_MESSAGE("Testing " << c[l]
                                       << " for dxcnt = " << dxcnt
                                       << ", dzcnt = " << dzcnt
                                       << " using " << nproc << " processor"
                                       << (nproc > 1 ? "s" : "") );
                }
                MPIBarrierRAII barrier(MPI_COMM_WORLD);
                test_accumulate_helper(
                    pg, dxcnt, dzcnt, 4*M_PI, 2, 4*M_PI/3, c[l][1], c[l][3]);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
