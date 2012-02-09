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
// explicit_op.hpp: Operators for channel_explicit
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "nonlinear.hpp"
#include "hybrid_op.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#pragma warning(disable:383 1572)

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

HybridIsothermalLinearOperator::HybridIsothermalLinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
    : suzerain::OperatorBase<real_t>(scenario, grid, dgrid, b, bop),
      common(common)
{
    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;
}

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
    Eigen::Vector3c tmp(Ny*field::count);
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
    SUZERAIN_UNUSED(substep_index);

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
    const bool has_zero_zero_modes = dgrid.has_zero_zero_modes();
    const real_t twopioverLx = twopiover(scenario.Lx);  // Weird looking...
    const real_t twopioverLz = twopiover(scenario.Lz);  // ...for FP control

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);

    // Means of the implicit momentum and energy forcing coefficients(!) are
    // maintained across each individual time step for sampling the statistics
    // /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock via
    // ILowStorageMethod::iota a la mean += iota*(sample - mean).  Note that
    // the sample must be divided by delta_t to account for step sizes.

    // See channel_treatment writeup for information on the steps below.
    // Steps may appear out of order relative to the writeup.
    // For phi real, the implicit operator is purely real for km == kn == 0.
    // Success of the zero-zero mode treatment relies on real-valued phi!
    SUZERAIN_ENSURE(phi.imag() == 0);

    // channel_treatment step (1) done during nonlinear operator application
    // via shared OperatorCommonBlock storage space.

    // Integral constraints enabled only when parameters are non-inf, non-NaN.
    // Allow disabling these to meet manufactured solution verification needs.
    const bool constrain_bulk_rhou
            = (boost::math::isnormal)(scenario.bulk_rhou);
    const bool constrain_bulk_rho
            = (boost::math::isnormal)(scenario.bulk_rho);

    // channel_treatment step (2) loads ones and mean streamwise velocity at
    // collocation points into the imaginary part of the constant (zero zero)
    // mode coefficients.  No forcing occurs at the lower or upper walls.
    if (constrain_bulk_rhou && has_zero_zero_modes) {

        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        mean_rhou.imag()[wall_lower] = 0;
        mean_rhou.imag().segment(1, Ny-2).setOnes();
        mean_rhou.imag()[wall_upper] = 0;

        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.imag()[wall_lower]      = 0;
        mean_rhoe.imag().segment(1, Ny-2) = common.u().segment(1, Ny-2);
        mean_rhoe.imag()[wall_upper]      = 0;
    }

    // B-spline numerics are not conservative when the mean density field has
    // any asymmetry.  Combat the bulk density's very, very small tendency to
    // drift by adding an integral constraint that the bulk density stay fixed
    // at a target value.  Approach follows that of the bulk momentum forcing.
    if (constrain_bulk_rho && has_zero_zero_modes) {
        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        mean_rho.imag().setConstant(1);
    }

    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions
    //
    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.
    //
    // channel_treatment step (9) sets isothermal conditions on wall
    // collocation points using e_wall = rho_wall / (gamma * (gamma - 1)).
    const real_t inv_gamma_gamma1 = 1/(scenario.gamma*(scenario.gamma - 1));

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

            // Apply isothermal, no slip conditions to tuple (papt, b)
            // FIXME Implement boundary conditions

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
            suzerain_lapack_zgbsvx(/* fact  */ 'E',
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

            // p := P^T x
            suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1.0, x, 1, 0.0, p, 1);
        }
    }

    if (constrain_bulk_rhou && has_zero_zero_modes) {

        // channel_treatment steps (4), (5), and (6) determine and
        // apply the appropriate bulk momentum forcing to achieve
        // a target value.
        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou.matrix());
        const real_t varphi = (scenario.bulk_rhou - bulk.real()) / bulk.imag();
        mean_rhou.real() += varphi * mean_rhou.imag();
        common.f() += iota*((varphi/delta_t)*mean_rhou.imag() - common.f());
        mean_rhou.imag().setZero();

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() += (varphi*scenario.Ma*scenario.Ma) * mean_rhoe.imag();
        common.f_dot_u() += iota*(   (varphi/delta_t)*mean_rhoe.imag()
                                   - common.f_dot_u());
        mean_rhoe.imag().setZero();

        // channel_treatment steps (8) and (9) already performed above
    }

    // Complete the bulk density target forcing
    if (constrain_bulk_rho && has_zero_zero_modes) {

        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rho.matrix());
        const real_t varphi = (scenario.bulk_rho - bulk.real()) / bulk.imag();
        mean_rho.real() += varphi * mean_rho.imag();
        mean_rho.imag().setZero();

        // Note that only the continuity equation is forced which neglects
        // contributions to the momentum and energy equations.  Neglecting
        // these terms is acceptable as the forcing is very small and vanishes
        // when the mean field is symmetric in the wall-normal direction.

        // Statistics for this quantity are not tracked.
    }

    // No volumetric energy forcing in performed current substep
    common.qb() += iota*(/* (Zero/delta_t) */ - common.qb());

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
