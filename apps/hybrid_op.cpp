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
// explicit_op.hpp: Operators for channel simulation
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "nonlinear.hpp"
#include "hybrid_op.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#pragma warning(disable:383 1572)

#define COUNTOF(expr) SUZERAIN_COUNTOF(expr)

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
real_t twopiover(const real_t L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

namespace channel {

void HybridIsothermalLinearOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index) const
{
    using suzerain::inorder::wavenumber;
    namespace field = channel::field;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    // const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    // const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const real_t twopioverLx = twopiover(scenario.Lx);  // Weird looking...
    const real_t twopioverLz = twopiover(scenario.Lz);  // ...for FP control

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);

    // Scratch for "in-place" suzerain_rholut_imexop_accumulate usage
    Eigen::VectorXc tmp(Ny*field::count);
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Iterate across local wavenumbers and apply operator "in-place".
    // Does not shortcircuit on only-dealiased state (TODO should it?)
    for (int n = dkbz; n < dkez; ++n) {
        const real_t kn = twopioverLz*wavenumber(dNz, n);
        for (int m = dkbx; m < dkex; ++m) {
            const real_t km = twopioverLx*wavenumber(dNx, m);

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Copy pencil into temporary storage
            suzerain::blas::copy(field::count*Ny, p, 1, tmp.data(), 1);

            // Accumulate result back into state storage
            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    tmp.data() + field::ndx::rho *Ny,
                    tmp.data() + field::ndx::rhou*Ny,
                    tmp.data() + field::ndx::rhov*Ny,
                    tmp.data() + field::ndx::rhow*Ny,
                    tmp.data() + field::ndx::rhoe*Ny,
                    0.0,
                    p + field::ndx::rho *Ny,
                    p + field::ndx::rhou*Ny,
                    p + field::ndx::rhov*Ny,
                    p + field::ndx::rhow*Ny,
                    p + field::ndx::rhoe*Ny);
        }
    }
}

void HybridIsothermalLinearOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output,
        const component delta_t,
        const std::size_t substep_index) const
{
    using suzerain::inorder::wavenumber;
    namespace field = channel::field;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    // const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    // const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const real_t twopioverLx = twopiover(scenario.Lx);  // Weird looking...
    const real_t twopioverLz = twopiover(scenario.Lz);  // ...for FP control

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == input.shape()[1])) return;

    // Input and output state storage has contiguous wall-normal scalars?
    // Furthermore, input has contiguous wall-normal pencils of all state?
    SUZERAIN_ENSURE(output.isIsomorphic(input));
    SUZERAIN_ENSURE(input.shape()   [0] ==  field::count);
    SUZERAIN_ENSURE(input.strides() [0] == (unsigned) Ny);
    SUZERAIN_ENSURE(input.shape()   [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(input.strides() [1] ==             1);
    SUZERAIN_ENSURE(output.strides()[1] ==             1);

    // Scratch for suzerain_rholut_imexop_accumulate usage
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Iterate across local wavenumbers and apply operator "in-place".
    // Does not shortcircuit on only-dealiased state (TODO should it?)
    for (int n = dkbz; n < dkez; ++n) {
        const real_t kn = twopioverLz*wavenumber(dNz, n);
        for (int m = dkbx; m < dkex; ++m) {
            const real_t km = twopioverLx*wavenumber(dNx, m);

            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    &input [field::ndx::rho ][0][m - dkbx][n - dkbz],
                    &input [field::ndx::rhou][0][m - dkbx][n - dkbz],
                    &input [field::ndx::rhov][0][m - dkbx][n - dkbz],
                    &input [field::ndx::rhow][0][m - dkbx][n - dkbz],
                    &input [field::ndx::rhoe][0][m - dkbx][n - dkbz],
                    beta,
                    &output[field::ndx::rho ][0][m - dkbx][n - dkbz],
                    &output[field::ndx::rhou][0][m - dkbx][n - dkbz],
                    &output[field::ndx::rhov][0][m - dkbx][n - dkbz],
                    &output[field::ndx::rhow][0][m - dkbx][n - dkbz],
                    &output[field::ndx::rhoe][0][m - dkbx][n - dkbz]);

        }
    }
}

void HybridIsothermalLinearOperator::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    // Shorthand
    using Eigen::Map;
    using Eigen::ArrayXc;
    using Eigen::ArrayXi;
    using Eigen::ArrayXr;
    using Eigen::ArrayXXc;
    using suzerain::inorder::wavenumber;
    namespace field = channel::field;
    namespace ndx   = field::ndx;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);
    SUZERAIN_UNUSED(iota);

    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    // const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    // const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const std::size_t wall_lower  = 0;
    const std::size_t wall_upper  = Ny - 1;
    const real_t twopioverLx = twopiover(scenario.Lx);  // Weird looking...
    const real_t twopioverLz = twopiover(scenario.Lz);  // ...for FP control

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);

    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions

    // Scratch for suzerain_rholut_imexop-based "inversion" using ZGBSVX
    suzerain_bsmbsm A = suzerain_bsmbsm_construct(
            (int) field::count, Ny, bop.max_kl(), bop.max_ku());
    ArrayXXc papt(A.LD,      A.N);                     // Holds PAP^T
    ArrayXXc lu  (A.LD+A.KL, A.N);                     // Holds LU of PAP^T
    ArrayXc buf(std::max(A.ld*A.n, 4*A.N));            // For packc, b, x, work
    complex_t * const b    = buf.data();               // ...parcel out b
    complex_t * const x    = buf.data() +   A.N;       // ...parcel out x
    complex_t * const work = buf.data() + 2*A.N;       // ...parcel out work
    ArrayXr rcrwork(3*A.N);                            // For r, c, rwork
    real_t * const r     = rcrwork.data();             // ...parcel out r
    real_t * const c     = rcrwork.data() +   A.N;     // ...parcel out c
    real_t * const rwork = rcrwork.data() + 2*A.N;     // ...parcel out rwork
    ArrayXi ipiv(A.N);                                 // For ipiv

    // Pack reference details for suzerain_rholut_imexop routines
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Prepare mapped PAP^T indices for boundary condition implementation
    const int papt_noslip[2][3] = {
        {
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhou)*A.n + wall_lower ),
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhov)*A.n + wall_lower ),
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhow)*A.n + wall_lower )
        }, {
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhou)*A.n + wall_upper ),
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhov)*A.n + wall_upper ),
            suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhow)*A.n + wall_upper )
        }
    };
    const int papt_rho[2] = {
        suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rho )*A.n + wall_lower ),
        suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rho )*A.n + wall_upper ),
    };
    const int papt_rhoe[2] = {
        suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhoe)*A.n + wall_lower ),
        suzerain_bsmbsm_qinv(A.S, A.n, ((int) ndx::rhoe)*A.n + wall_upper ),
    };

    // Iterate across local wavenumbers and "invert" operator "in-place".
    // Does not shortcircuit on only-dealiased state (TODO should it?)
    for (int n = dkbz; n < dkez; ++n) {
        const real_t kn = twopioverLz*wavenumber(dNz, n);
        for (int m = dkbx; m < dkex; ++m) {
            const real_t km = twopioverLx*wavenumber(dNx, m);

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Form wavenumber-dependent PAP^T within papt
            suzerain_rholut_imexop_packc(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    ndx::rho, ndx::rhou, ndx::rhov, ndx::rhow, ndx::rhoe,
                    buf.data(), &A, papt.data());

            // b := P p
            suzerain_bsmbsm_zaPxpby('N', A.S, A.n, 1.0, p, 1, 0.0, b, 1);

            // channel_treatment step (8) sets no-slip conditions
            // on wall collocation points.
            //
            // channel_treatment step (9) sets isothermal conditions at walls
            // using rho_wall = e_wall * gamma * (gamma - 1).
            //
            // Access is one pass per row (w/ others hopefully nearby in cache).
            // Attempt made to not unnecessarily disturb matrix conditioning.
            // With apologies to whomever next needs to decipher this horror.

            // Zero RHS of momentum and energy equations at lower, upper walls
            for (size_t wall = 0; wall < 2; ++wall) {
                for (size_t eqn = 0; eqn < COUNTOF(papt_noslip[wall]); ++eqn) {
                    b[papt_noslip[wall][eqn]] = 0;
                }
                b[papt_rhoe[wall]] = 0;
            }

            // Modify the equations within PAP^T for lower, upper walls
            for (size_t wall = 0; wall < 2; ++wall) {
                int begin, end, inc;

                // Zero all row entries but the mass matrix one on the diagonal
                // Then momentum equations are like const*rho{u,v,w} = 0.
                for (size_t eqn = 0; eqn < COUNTOF(papt_noslip[wall]); ++eqn) {
                    complex_t * const row = (complex_t *) suzerain_gbmatrix_row(
                            A.N, A.N, A.KL, A.KU, (void *) papt.data(), A.LD,
                            sizeof(complex_t), papt_noslip[wall][eqn],
                            &begin, &end, &inc);
                    for (int rowndx = begin; rowndx < end; ++rowndx) {
                        if (rowndx != papt_noslip[wall][eqn]) {
                            row[rowndx*inc] = 0;
                        }
                    }
                }

                // Set constraint const*rho - const*gamma*(gamma-1)*rhoe = 0
                complex_t * const row = (complex_t *) suzerain_gbmatrix_row(
                        A.N, A.N, A.KL, A.KU, (void *) papt.data(), A.LD,
                        sizeof(complex_t), papt_rhoe[wall], &begin, &end, &inc);
                for (int rowndx = begin; rowndx < end; ++rowndx) {
                    if (rowndx == papt_rho[wall]) {
                        // NOP
                    } else if (rowndx == papt_rhoe[wall]) {
                        row[rowndx*inc] = row[papt_rho[wall]*inc]
                                        * s.gamma*(1 - s.gamma);
                    } else {
                        row[rowndx*inc] = 0;
                    }
                }

            }

            // x := (LU)^-1 b using GBSVX
            //
            // "Why GBSVX?" you ask.  "It's expensive!" you rightly observe.
            //
            // Well, we might as well get our money's worth out of the
            // factorization once we've paid for it.  If GBSVX is too
            // expensive, GBSV is a less intensive option.  I have unfounded
            // hunches that (a) the incremental cost from GBSV to GBSVX is low
            // once everything is in L2 cache and (b) always solving the
            // equations well using iterative refinement will save us from the
            // woes of modest resolution and possibly less-than-perfectly
            // conditioned operators stemming from
            // reference-values-during-transients and high wavenumbers.
            //
            // TODO Track and report statistics on rcond, ferr, and berr
            real_t rcond, ferr, berr;
            char equed;
            const int info = suzerain_lapack_zgbsvx(/* fact  */ 'E',
                                                    /* trans */ 'N',
                                                    /* n     */ A.N,
                                                    /* kl    */ A.KL,
                                                    /* ku    */ A.KU,
                                                    /* nrhs  */ 1,
                                                    /* ab    */ papt.data(),
                                                    /* ldab  */ A.LD,
                                                    /* afb   */ lu.data(),
                                                    /* ldafb */ A.LD + A.KL,
                                                    /* ipiv  */ ipiv.data(),
                                                    /* equed */ &equed,
                                                    /* r     */ r,
                                                    /* c     */ c,
                                                    /* b     */ b,
                                                    /* ldb   */ A.N,
                                                    /* x     */ x,
                                                    /* ldx   */ A.N,
                                                    /* rcond */ &rcond,
                                                    /* ferr  */ &ferr,
                                                    /* berr  */ &berr,
                                                    /* work  */ work,
                                                    /* rwork */ rwork);
            if (info) {
                char buffer[80];
                snprintf(buffer, sizeof(buffer),
                        "suzerain_lapack_zgbsvx reported error %d", info);
#ifndef NDEBUG
                std::cerr << buffer << ": Banded PAP^T = " << std::endl;
                std::cerr << papt << std::endl;
#endif
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            }

            // p := P^T x
            suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1.0, x, 1, 0.0, p, 1);
        }
    }

    // State leaves method as coefficients in X, Y, and Z directions
}

std::vector<real_t> HybridNonlinearOperator::applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return channel::applyNonlinearOperator<true,  channel::linearize::rhome>
            (*this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return channel::applyNonlinearOperator<false, channel::linearize::rhome>
            (*this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} // end namespace channel
