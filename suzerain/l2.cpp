//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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
 * @copydoc l2.hpp
 */

#include <suzerain/l2.hpp>

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/error.h>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

std::vector<field_L2xyz>
compute_field_L2xyz(
        const contiguous_state<4,complex_t> &state,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& gop)
{
    // Ensure state storage meets this routine's assumptions
    // Notice state.shape()[0] may be any value
    using boost::numeric_cast;
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed Galerkin field_L2xyz operator matrices
    SUZERAIN_ENSURE(gop.get()->method==SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only want non-dealiased X-direction modes to contribute to L2
    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Only want non-dealiased Z-direction modes to contribute to L2
    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Temporary storage for inner product computations
    //
    // Though u^H M u is always real for SPD M, floating point errors may
    // produce a slightly non-SPD M so we permit complex intermediate results
    // up until the return structure is packed.
    VectorXc tmp;
    tmp.setZero(grid.N.y());

    // Contiguous temporary storage for accumulating and broadcasting results
    VectorXc buf(2*state.shape()[0]);
    Map<VectorXc> total2(buf.data() + 0,                state.shape()[0]);
    Map<VectorXc> mean2 (buf.data() + state.shape()[0], state.shape()[0]);

    // Compute the local L2 contribution towards each L^2_xyz norm squared
    // Computation uses partial sums at each loop to reduce swamping which is
    // more-or-less recursive summation using large partitioning factors.
    total2.setZero();
    for (size_t k = 0; k < state.shape()[0]; ++k) {
        complex_t jsum = 0;
        for (int j = 0; j < 2; ++j) {
            complex_t nsum = 0;
            for (int n = mzb[j]; n < mze[j]; ++n) {
                complex_t msum = 0;
                for (int m = mxb[0]; m < mxe[0]; ++m) {
                    const complex_t * u_mn
                        = &state[k][0][m - dgrid.local_wave_start.x()]
                                      [n - dgrid.local_wave_start.z()];
                    gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
                    complex_t dot = blas::dot(grid.N.y(), u_mn,       1,
                                                          tmp.data(), 1);
                    if (m > 0 && m < grid.dN.x()/2) {
                        dot *= 2;
                    }
                    msum += dot;
                }
                nsum += msum;
            }
            jsum += nsum;
        }
        total2[k] += jsum;
    }
    total2 *= grid.L.x() * grid.L.z();

    // Compute the mean contributions to each L^2_xyz norm squared.
    if (dgrid.has_zero_zero_modes()) {
        for (size_t k = 0; k < state.shape()[0]; ++k) {
            const complex_t * u_mn = &state[k][0][0][0];
            gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
            mean2[k] = blas::dot(grid.N.y(), u_mn, 1, tmp.data(), 1);
        }
        mean2 *= grid.L.x() * grid.L.z();
    } else {
        mean2.setZero();
    }

    // Allreduce total2 and mean2 simultaneously.
    // (Though it moves more overall data, Allreduce chosen as Redmine #3079
    // suspected a race-related hang from earlier Reduce/Bcast implementation).
    //
    // Reduction operation is complex-valued, but OpenMPI pre-1.7.3 can bomb
    // unless we keep it real: https://svn.open-mpi.org/trac/ompi/ticket/3127.
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, buf.data(),
            (sizeof(complex_t)/sizeof(real_t))*buf.size(),
            mpi::datatype<real_t>(), MPI_SUM, MPI_COMM_WORLD));

    // Obtain fluctuating2 = total2 - mean2 and pack the return structure
    std::vector<field_L2xyz> retval(state.shape()[0]);
    for (size_t k = 0; k < retval.size(); ++k) {
        retval[k].mean        = std::sqrt(std::abs(mean2[k]            ));
        retval[k].fluctuating = std::sqrt(std::abs(total2[k] - mean2[k]));
    }

    return retval;
}

//
// Implementation very much resembles compute_field_L2xyz.
// Begin by understanding that first and then compare the differences.
//
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop)
{
    // Ensure state storage meets this routine's assumptions
    // Notice state.shape()[0] may be any value
    using boost::numeric_cast;
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed collocation-based operator matrices
    SUZERAIN_ENSURE(cop.get()->method==SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    // Only want non-dealiased X-direction modes to contribute to L2
    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Only want non-dealiased Z-direction modes to contribute to L2
    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Temporary complex-valued storage for matrix-vector products
    ArrayXc tmp;
    tmp.setZero(grid.N.y());

    // Contiguous temporary storage for accumulating and broadcasting results
    // Fast columns store wall-normal swaths of data.
    // Real-valued since cwiseAbs2 always produces real-valued results.
    ArrayXXr buf(grid.N.y(), 2*state.shape()[0]);
    ArrayXXr::ColsBlockXpr total2 = buf.leftCols (state.shape()[0]);
    ArrayXXr::ColsBlockXpr mean2  = buf.rightCols(state.shape()[0]);

    // Compute the local L2 contribution towards each L^2_xz norm squared
    // Computation uses partial sums at each loop to reduce swamping which is
    // more-or-less recursive summation using large partitioning factors.
    ArrayXr jsum, nsum, msum;
    total2.setZero();
    for (size_t k = 0; k < state.shape()[0]; ++k) {
        jsum.setZero(grid.N.y());
        for (int j = 0; j < 2; ++j) {
            nsum.setZero(grid.N.y());
            for (int n = mzb[j]; n < mze[j]; ++n) {
                msum.setZero(grid.N.y());
                for (int m = mxb[0]; m < mxe[0]; ++m) {
                    const complex_t * u_mn
                        = &state[k][0][m - dgrid.local_wave_start.x()]
                                      [n - dgrid.local_wave_start.z()];
                    if (m > 0 && m < grid.dN.x()/2) {
                        // Conjugate symmetry uses root-2 so Abs2 below gives 2
                        cop.accumulate(0, M_SQRT2, u_mn, 1, 0., tmp.data(), 1);
                    } else {
                        cop.accumulate(0, 1.0,     u_mn, 1, 0., tmp.data(), 1);
                    }
                    msum += tmp.abs2();
                }
                nsum += msum;
            }
            jsum += nsum;
        }
        total2.col(k) += jsum;
    }
    total2 *= grid.L.x() * grid.L.z();

    // Compute the mean contributions to each L^2_xyz norm squared.
    if (dgrid.has_zero_zero_modes()) {
        for (size_t k = 0; k < state.shape()[0]; ++k) {
            const complex_t * u_mn = &state[k][0][0][0];
            cop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
            mean2.col(k) = tmp.cwiseAbs2();
        }
        mean2 *= grid.L.x() * grid.L.z();
    } else {
        mean2.setZero();
    }

    // Allreduce total2 and mean2 simultaneously.
    // (Though it moves more overall data, Allreduce chosen as Redmine #3079
    // suspected a race-related hang from earlier Reduce/Bcast implementation).
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, buf.data(),
            buf.size(), mpi::datatype<real_t>(), MPI_SUM, MPI_COMM_WORLD));

    // Obtain fluctuating2 = total2 - mean2 and pack the return structure
    std::vector<field_L2xz> retval(state.shape()[0]);
    for (size_t k = 0; k < retval.size(); ++k) {
        retval[k].mean        = mean2.col(k).sqrt();
        retval[k].fluctuating = (total2.col(k) - mean2.col(k)).sqrt();
    }

    return retval;
}

//
// Implementation very much resembles previous compute_field_L2xz
// with this routine assuming B-spline mass matrix already applied
//
std::vector<field_L2xz>
compute_field_L2xz(
        const contiguous_state<4,complex_t> &state,
        const specification_grid& grid,
        const pencil_grid& dgrid)
{
    // Ensure state storage meets this routine's assumptions
    // Notice state.shape()[0] may be any value
    using boost::numeric_cast;
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Only want non-dealiased X-direction modes to contribute to L2
    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Only want non-dealiased Z-direction modes to contribute to L2
    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Contiguous temporary storage for accumulating and broadcasting results
    // Fast columns store wall-normal swaths of data.
    // Real-valued since cwiseAbs2 always produces real-valued results.
    ArrayXXr buf(grid.N.y(), 2*state.shape()[0]);
    ArrayXXr::ColsBlockXpr total2 = buf.leftCols (state.shape()[0]);
    ArrayXXr::ColsBlockXpr mean2  = buf.rightCols(state.shape()[0]);

    // Compute the local L2 contribution towards each L^2_xz norm squared
    // Computation uses partial sums at each loop to reduce swamping which is
    // more-or-less recursive summation using large partitioning factors.
    ArrayXr jsum, nsum, msum;
    total2.setZero();
    for (size_t k = 0; k < state.shape()[0]; ++k) {
        jsum.setZero(grid.N.y());
        for (int j = 0; j < 2; ++j) {
            nsum.setZero(grid.N.y());
            for (int n = mzb[j]; n < mze[j]; ++n) {
                msum.setZero(grid.N.y());
                for (int m = mxb[0]; m < mxe[0]; ++m) {
                    Map<const ArrayXc> u_mn(
                            &state[k][0][m - dgrid.local_wave_start.x()]
                                        [n - dgrid.local_wave_start.z()],
                            grid.N.y());
                    if (m > 0 && m < grid.dN.x()/2) { // Conjugate symmetry
                        msum += 2*u_mn.abs2();
                    } else {
                        msum +=   u_mn.abs2();
                    }
                }
                nsum += msum;
            }
            jsum += nsum;
        }
        total2.col(k) += jsum;
    }
    total2 *= grid.L.x() * grid.L.z();

    // Compute the mean contributions to each L^2_xyz norm squared.
    if (dgrid.has_zero_zero_modes()) {
        for (size_t k = 0; k < state.shape()[0]; ++k) {
            Map<const ArrayXc> u_mn(&state[k][0][0][0], grid.N.y());
            mean2.col(k) = u_mn.cwiseAbs2();
        }
        mean2 *= grid.L.x() * grid.L.z();
    } else {
        mean2.setZero();
    }

    // Allreduce total2 and mean2 simultaneously.
    // (Though it moves more overall data, Allreduce chosen as Redmine #3079
    // suspected a race-related hang from earlier Reduce/Bcast implementation).
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(),
            mpi::datatype<real_t>(), MPI_SUM, MPI_COMM_WORLD));

    // Obtain fluctuating2 = total2 - mean2 and pack the return structure
    std::vector<field_L2xz> retval(state.shape()[0]);
    for (size_t k = 0; k < retval.size(); ++k) {
        retval[k].mean        = mean2.col(k).sqrt();
        retval[k].fluctuating = (total2.col(k) - mean2.col(k)).sqrt();
    }

    return retval;
}

void
compute_twopoint_xlocal(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index si,
        const contiguous_state<4,complex_t>::index sj,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        complex_t * const out)
{
    // Ensure state matches assumptions; state.shape()[0] may be any value
    using boost::numeric_cast;
    SUZERAIN_ENSURE((int) state.shape()[1] == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE((int) state.shape()[2] == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE((int) state.shape()[3] == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE(state.strides()[1] == (int) 1);
    SUZERAIN_ENSURE(state.strides()[2] == (int) state.shape()[1]);
    SUZERAIN_ENSURE(state.strides()[3] == (int) ( state.shape()[1]
                                                 *state.shape()[2]));

    // Grab relevant parallel decomposition information
    const int Ny   = dgrid.global_wave_extent.y();
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();

    // Prepare a view of the output storage and zero it prior to accumulation
    Map<MatrixXXc> o(out, Ny, Nx/2+1);
    o.setZero();

    // Sum rank-local contribution to Fourier-transformed two-point correlation
    // We /do/ sum dealiasing/Nyquist modes so clear a priori if that bugs you
    // (though notice only Nx/2+1 and not dNx/2+1 modes are reported in result)
    for (int n = dkbz; n < dkez; ++n) {
        for (int m = dkbx; m < dkex; ++m) {
            const int wm = inorder::wavenumber(dNx, m);
            if (wm >= o.cols()) continue;                            // Ignored

            Map<const VectorXc> u_mn(&state[si][0][m - dkbx][n - dkbz], Ny);
            Map<const VectorXc> v_mn(&state[sj][0][m - dkbx][n - dkbz], Ny);

            o.col(wm) += u_mn.conjugate().cwiseProduct(v_mn);
        }
    }

    // Ensure real-valued R_uv(x) arises from Fourier-transforming result
    // (occurs trivially for R_uu but dealiasing can break property for R_uv)
    o.col(0).imag().setZero();
    if (Nx%2 == 0) {
        o.col(Nx/2).imag().setZero();
    }
}

void
compute_twopoint_zlocal(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index si,
        const contiguous_state<4,complex_t>::index sj,
        const specification_grid& grid,
        const pencil_grid& dgrid,
        complex_t * const out)
{
    // Ensure state matches assumptions; state.shape()[0] may be any value
    using boost::numeric_cast;
    SUZERAIN_ENSURE((int) state.shape()[1] == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE((int) state.shape()[2] == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE((int) state.shape()[3] == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE(state.strides()[1] == (int) 1);
    SUZERAIN_ENSURE(state.strides()[2] == (int) state.shape()[1]);
    SUZERAIN_ENSURE(state.strides()[3] == (int) ( state.shape()[1]
                                                 *state.shape()[2]));

    // Grab relevant parallel decomposition information
    const int Ny   = dgrid.global_wave_extent.y();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();

    // Prepare a view of the output storage and zero it prior to accumulation
    Map<MatrixXXc> o(out, Ny, Nz);
    o.setZero();

    // Sum rank-local contribution to Fourier-transformed two-point correlation
    // We /do/ sum dealiasing/Nyquist modes so clear a priori if that bugs you
    for (int n = dkbz; n < dkez; ++n) {
        if (!inorder::wavenumber_translatable(Nz, dNz, n)) continue;  // Ignored
        const int wn_idx = inorder::index(Nz, inorder::wavenumber(dNz, n));

        for (int m = dkbx; m < dkex; ++m) {

            Map<const VectorXc> u_mn(&state[si][0][m - dkbx][n - dkbz], Ny);
            Map<const VectorXc> v_mn(&state[sj][0][m - dkbx][n - dkbz], Ny);

            if (m == 0) {
                o.col(wn_idx) +=   u_mn.conjugate().cwiseProduct(v_mn);
            } else if (inorder::wavenumber_nyquist_index(dNx, m)) {
                o.col(wn_idx) +=   u_mn.cwiseProduct(v_mn.conjugate());
            } else {
                // Admittedly, there are faster ways to compute this product
                o.col(wn_idx) +=   u_mn.conjugate().cwiseProduct(v_mn)
                                 + u_mn.cwiseProduct(v_mn.conjugate());
            }
        }
    }

    // Ensure real-valued R_uv(z) arises from Fourier-transforming result
    // (occurs trivially for R_uu but dealiasing can break property for R_uv)
    o.col(0).imag().setZero();
    if (Nz%2 == 0) {
        o.col(Nz/2).imag().setZero();
    }
}

shared_array<complex_t>
compute_twopoint_x(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index nf,
        const specification_grid& grid,
        const pencil_grid& dgrid)
{
    SUZERAIN_ENSURE((int) state.shape()[0] <= (int) nf);

    // Allocate contiguous storage for all the pairwise results
    // As two-point is real-valued, only positive wavenumbers computed
    const int npairs  = (nf*(nf+1))/2;
    const int bufpair = grid.N.y() * (grid.N.x()/2+1);
    shared_array<complex_t> retval(
            (complex_t*)suzerain_blas_malloc(sizeof(complex_t)*bufpair*npairs),
            suzerain_blas_free);

    // Compute result for each unique pair of scalar fields
    for (int si = 0; si < nf; ++si) {
        for (int sj = si; sj < nf; ++sj) {
            const int ndxpair = nf*si + sj - (si*(si+1))/2;
            compute_twopoint_xlocal(state, si, sj, grid, dgrid,
                                    &retval[bufpair*ndxpair]);
        }
    }

    // Perform the global summation across all pairs at once
    //
    // Reduction operation is complex-valued, but OpenMPI pre-1.7.3 can bomb
    // unless we keep it real: https://svn.open-mpi.org/trac/ompi/ticket/3127.
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, retval.get(),
            (sizeof(complex_t)/sizeof(real_t))*bufpair*npairs,
            mpi::datatype<real_t>(), MPI_SUM, MPI_COMM_WORLD));

    return retval;
}

shared_array<complex_t>
compute_twopoint_z(
        const contiguous_state<4,complex_t> &state,
        const contiguous_state<4,complex_t>::index nf,
        const specification_grid& grid,
        const pencil_grid& dgrid)
{
    SUZERAIN_ENSURE((int) state.shape()[0] <= (int) nf);

    // Allocate contiguous storage for all the pairwise results
    // As two-point is real-valued, only positive wavenumbers computed
    const int npairs  = (nf*(nf+1))/2;
    const int bufpair = grid.N.y() * grid.N.z();
    shared_array<complex_t> retval(
            (complex_t*)suzerain_blas_malloc(sizeof(complex_t)*bufpair*npairs),
            suzerain_blas_free);

    // Compute result for each unique pair of scalar fields
    for (int si = 0; si < nf; ++si) {
        for (int sj = si; sj < nf; ++sj) {
            const int ndxpair = nf*si + sj - (si*(si+1))/2;
            compute_twopoint_zlocal(state, si, sj, grid, dgrid,
                                    &retval[bufpair*ndxpair]);
        }
    }

    // Perform the global summation across all pairs at once
    //
    // Reduction operation is complex-valued, but OpenMPI pre-1.7.3 can bomb
    // unless we keep it real: https://svn.open-mpi.org/trac/ompi/ticket/3127.
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, retval.get(),
            (sizeof(complex_t)/sizeof(real_t))*bufpair*npairs,
            mpi::datatype<real_t>(), MPI_SUM, MPI_COMM_WORLD));

    return retval;
}

} // end namespace suzerain
