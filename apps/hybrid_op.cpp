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

// Must appear before Suzerain header includes to obtain timer hooks
#include "timers.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/complex.hpp>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#include "hybrid_op.hpp"
#include "logging.hpp"
#include "nonlinear.hpp"

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

void HybridIsothermalLinearOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index) const
{
    GRVY_TIMER_BEGIN("applyMassPlusScaledOperator");

    using suzerain::inorder::wavenumber;
    using suzerain::inorder::wavenumber_abs;
    using suzerain::inorder::wavenumber_max;
    namespace field = channel::field;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
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
    // Short circuits on wavenumbers present only for dealiasing
    for (int n = dkbz; n < dkez; ++n) {
        if (wavenumber_abs(dNz, n) > wavenumber_max(Nz)) continue;
        const real_t kn = twopioverLz*wavenumber(dNz, n);
        for (int m = dkbx; m < dkex; ++m) {
            if (wavenumber_abs(dNx, m) > wavenumber_max(Nx)) continue;
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

    GRVY_TIMER_END("applyMassPlusScaledOperator");
}

void HybridIsothermalLinearOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output,
        const component delta_t,
        const std::size_t substep_index) const
{
    GRVY_TIMER_BEGIN("accumulateMassPlusScaledOperator");

    using suzerain::inorder::wavenumber;
    using suzerain::inorder::wavenumber_abs;
    using suzerain::inorder::wavenumber_max;
    namespace field = channel::field;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
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
    // Short circuits on wavenumbers present only for dealiasing
    for (int n = dkbz; n < dkez; ++n) {
        if (wavenumber_abs(dNz, n) > wavenumber_max(Nz)) continue;
        const real_t kn = twopioverLz*wavenumber(dNz, n);
        for (int m = dkbx; m < dkex; ++m) {
            if (wavenumber_abs(dNx, m) > wavenumber_max(Nx)) continue;
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

    GRVY_TIMER_END("accumulateMassPlusScaledOperator");
}

/**
 * A functor for applying isothermal conditions for PAP^T-based operators.
 * Encapsulated in a class so that the same logic may be applied to real-
 * and complex-valued operators and right hand sides.
 *
 * Helps accomplish the following steps:
 * <ul>
 * <li>
 *     channel_treatment step (8) sets no-slip conditions
 *     on wall collocation points.
 * </li><li>
 * </li>
 *     channel_treatment step (9) sets isothermal conditions at walls
 *     using rho_wall = e_wall * gamma * (gamma - 1).
 * </ul>
 */
class IsothermalNoSlipPAPTEnforcer
{
    enum { nwalls = 2, nmomentum = 3 };

    // Indices within PAP^T at which to apply boundary conditions
    // Computed once within constructor and then repeatedly used
    int rho[nwalls], noslip[nwalls][nmomentum], rhoe[nwalls];

    // Precomputed coefficient based on the isothermal equation of state
    real_t gamma_times_one_minus_gamma;

public:

    IsothermalNoSlipPAPTEnforcer(const suzerain_bsmbsm &A,
                                 const suzerain_rholut_imexop_scenario &s)
        : gamma_times_one_minus_gamma(s.gamma * (1 - s.gamma))
    {
        // Starting offset to named scalars in InterleavedState pencil
        namespace ndx = channel::field::ndx;
        const int start_rho  = static_cast<int>(ndx::rho )*A.n;
        const int start_rhou = static_cast<int>(ndx::rhou)*A.n;
        const int start_rhov = static_cast<int>(ndx::rhov)*A.n;
        const int start_rhow = static_cast<int>(ndx::rhow)*A.n;
        const int start_rhoe = static_cast<int>(ndx::rhoe)*A.n;

        // Relative to start_foo what is the offset to lower, upper walls
        const int wall[nwalls] = { 0, A.n - 1};

        // Prepare indices within PAP^T corresponding to the walls.
        // Uses that A_{i,j} maps to {PAP^T}_{{q^-1}(i),{q^(-1)}(j)}.
        for (int i = 0; i < nwalls; ++i) {
            rho   [i]    = suzerain_bsmbsm_qinv(A.S, A.n, start_rho +wall[i]);
            noslip[i][0] = suzerain_bsmbsm_qinv(A.S, A.n, start_rhou+wall[i]);
            noslip[i][1] = suzerain_bsmbsm_qinv(A.S, A.n, start_rhov+wall[i]);
            noslip[i][2] = suzerain_bsmbsm_qinv(A.S, A.n, start_rhow+wall[i]);
            rhoe  [i]    = suzerain_bsmbsm_qinv(A.S, A.n, start_rhoe+wall[i]);
        }
    }

    /**
     * Zero RHS of momentum and energy equations at lower, upper walls.
     * Must be used in conjunction with op().
     */
    template<typename T>
    void rhs(T * const b)
    {
        for (int wall = 0; wall < nwalls; ++wall) {
            for (int eqn = 0; eqn < nmomentum; ++eqn) {
                b[noslip[wall][eqn]] = 0;
            }
            b[rhoe[wall]] = 0;
        }
    }

    /**
     * Modify the equations within PAP^T for lower, upper walls.
     * Must be done in conjunction with rhs().
     */
    template<typename T>
    void op(const suzerain_bsmbsm& A, T * const papt, int papt_ld)
    {
        // Access is one pass per row (w/ others hopefully nearby in cache).
        // Attempt made to not unnecessarily disturb matrix conditioning.

        for (int wall = 0; wall < nwalls; ++wall) {
            int begin, end, inc;

            // Zero all row entries but the mass matrix one on the diagonal
            // Then momentum equations are like const*rho{u,v,w} = 0.
            for (size_t eqn = 0; eqn < nmomentum; ++eqn) {
                T * const row = (T *) suzerain_gbmatrix_row(
                        A.N, A.N, A.KL, A.KU, (void *) papt, papt_ld,
                        sizeof(T), noslip[wall][eqn], &begin, &end, &inc);
                for (int rowndx = begin; rowndx < end; ++rowndx) {
                    if (rowndx != noslip[wall][eqn]) {
                        row[rowndx*inc] = 0;
                    }
                }
            }

            // Set constraint const*rho - const*gamma*(gamma-1)*rhoe = 0
            T * const row = (T *) suzerain_gbmatrix_row(
                    A.N, A.N, A.KL, A.KU, (void *) papt, papt_ld,
                    sizeof(T), rhoe[wall], &begin, &end, &inc);
            // Necessary to ensure constraint possible in degenerate case
            complex_t &rhocoeff = row[rho[wall]*inc];
            if (rhocoeff == complex_t(0)) rhocoeff = 1;
            // Scan row and adjust coefficients for constraint
            for (int rowndx = begin; rowndx < end; ++rowndx) {
                if (rowndx == rho[wall]) {
                    // NOP
                } else if (rowndx == rhoe[wall]) {
                    row[rowndx*inc] = rhocoeff*gamma_times_one_minus_gamma;
                } else {
                    row[rowndx*inc] = 0;
                }
            }
        }
    }
};

void HybridIsothermalLinearOperator::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    GRVY_TIMER_BEGIN("invertMassPlusScaledOperator");

    // Shorthand
    using suzerain::inorder::wavenumber;
    using suzerain::inorder::wavenumber_abs;
    using suzerain::inorder::wavenumber_max;
    namespace field = channel::field;
    namespace ndx   = field::ndx;
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);
    SUZERAIN_UNUSED(iota);

    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    const int Ny   = dgrid.global_wave_extent.y();
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
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

    // Details for suzerain_rholut_imexop-based "inversion" using ?GBSVX
    // Macros used to automatically increase paranoia during debug builds
    suzerain_bsmbsm A = suzerain_bsmbsm_construct(
            (int) field::count, Ny, bop.max_kl(), bop.max_ku());
#ifndef NDEBUG
# define SCRATCH_C(type, name, ...) \
         type name = type::Constant(__VA_ARGS__, suzerain::complex::NaN<type::Scalar>())
# define SCRATCH_R(type, name, ...) \
         type name = type::Constant(__VA_ARGS__, std::numeric_limits<type::Scalar>::quiet_NaN())
# define SCRATCH_I(type, name, ...) \
         type name = type::Constant(__VA_ARGS__, -12345)
#else
# define SCRATCH_C(type, name, ...) type name(__VA_ARGS__)
# define SCRATCH_R(type, name, ...) type name(__VA_ARGS__)
# define SCRATCH_I(type, name, ...) type name(__VA_ARGS__)
#endif
    SCRATCH_C(Eigen::ArrayXXc, buf,     A.ld,      A.n);  // For packc calls
    SCRATCH_C(Eigen::ArrayXXc, papt,    A.LD,      A.N);  // Holds PAP^T
    SCRATCH_C(Eigen::ArrayXXc, lu,      A.LD+A.KL, A.N);  // Holds LU of PAP^T
    SCRATCH_I(Eigen::ArrayXi,  ipiv,    A.N);             // Linear solve...
    SCRATCH_R(Eigen::ArrayXr,  r,       A.N);
    SCRATCH_R(Eigen::ArrayXr,  c,       A.N);
    SCRATCH_C(Eigen::ArrayXc,  b,       A.N);
    SCRATCH_C(Eigen::ArrayXc,  x,       A.N);
    SCRATCH_C(Eigen::ArrayXc,  work,  2*A.N);
    SCRATCH_R(Eigen::ArrayXr,  rwork,   A.N);
#undef SCRATCH_C
#undef SCRATCH_R
#undef SCRATCH_I

    // Pack reference details for suzerain_rholut_imexop routines
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Prepare an almost functor mutating RHS and PAP^T to enforce BCs.
    IsothermalNoSlipPAPTEnforcer bc_enforcer(A, s);

    // How will we solve the linear system of equations?
    enum solve_types { gbsvx, gbsv };
    static const solve_types solve_type = gbsvx;

    // Iterate across local wavenumbers and "invert" operator "in-place".
    // Short circuits on wavenumbers present only for dealiasing
    for (int n = dkbz; n < dkez; ++n) {
        if (wavenumber_abs(dNz, n) > wavenumber_max(Nz)) continue;
        const real_t kn = twopioverLz*wavenumber(dNz, n);

        for (int m = dkbx; m < dkex; ++m) {
            if (wavenumber_abs(dNx, m) > wavenumber_max(Nx)) continue;
            const real_t km = twopioverLx*wavenumber(dNx, m);

            // Form complex-valued, wavenumber-dependent PAP^T within papt
            GRVY_TIMER_BEGIN("implicit operator assembly");
            switch (solve_type) {
            case gbsvx:
                suzerain_rholut_imexop_packc(
                        phi, km, kn, &s, &ref, &ld, bop.get(),
                        ndx::rho, ndx::rhou, ndx::rhov, ndx::rhow, ndx::rhoe,
                        buf.data(), &A, papt.data());
                break;
            case gbsv:
                suzerain_rholut_imexop_packf(
                        phi, km, kn, &s, &ref, &ld, bop.get(),
                        ndx::rho, ndx::rhou, ndx::rhov, ndx::rhow, ndx::rhoe,
                        buf.data(), &A, lu.data());
                break;
            }
            GRVY_TIMER_END("implicit operator assembly");

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Given state pencil "p" the rest of the solve loop looks like
            //
            //     b := P p          using suzerain_bsmbsm_?aPxpby
            //     apply BC to RHS   using IsothermalNoSlipPAPTEnforcer
            //     apply BC to PAP^T using IsothermalNoSlipPAPTEnforcer
            //     x := (LU)^-1 b    using ?gbsvx which factorizes PAP^T
            //                       or ?gbsv which factorizes in place
            //     p := P^T x        using suzerain_bsmbsm_?aPxpby
            //
            // "Why gbsvx?" you ask.  "It's expensive!" you rightly observe.
            //
            // Well, we might as well get our money's worth out of the
            // factorization once we've paid for it.  If gbsvx is too
            // expensive, gbsv is a less intensive option.  I have unfounded
            // hunches that (a) the incremental cost from gbsv to gbsvx is low
            // once everything is in L2 cache and (b) always solving the
            // equations well using iterative refinement will save us from the
            // woes of modest resolution and possibly less-than-perfectly
            // conditioned operators stemming (from reference values during
            // transients and high wavenumbers).  So we'll run with gbsvx until
            // it's demonstrably a bad idea.

            const char *method;               // Used for error reporting
            static const char fact    = 'E';  // Common inputs for ?gbsvx
            static const char trans   = 'N';
            int info;                         // Common outputs for ?gbsvx
            char equed;
            real_t rcond, ferr[1], berr[1];

            GRVY_TIMER_BEGIN("implicit operator solve");
            switch (solve_type) {
            case gbsvx:
                method = "zgbsvx";
                break;
            case gbsv:
                method = "zgbsv";
                break;
            }
            GRVY_TIMER_BEGIN("suzerain_bsmbsm_zaPxpby");
            suzerain_bsmbsm_zaPxpby('N', A.S, A.n, 1, p, 1, 0, b.data(), 1);
            GRVY_TIMER_END("suzerain_bsmbsm_zaPxpby");
            GRVY_TIMER_BEGIN("implicit operator BCs");
            bc_enforcer.rhs(b.data());
            switch (solve_type) {
                case gbsvx:
                    bc_enforcer.op(A, papt.data(), papt.colStride());
                    break;
                case gbsv:
                    bc_enforcer.op(A, lu.data() + A.KL, lu.colStride());
                    break;
            }
            GRVY_TIMER_END("implicit operator BCs");
            switch (solve_type) {
            case gbsvx:
                info = suzerain_lapack_zgbsvx(fact, trans, A.N, A.KL, A.KU, 1,
                    papt.data(), papt.colStride(), lu.data(), lu.colStride(),
                    ipiv.data(), &equed, r.data(), c.data(),
                    b.data(), A.N, x.data(), A.N,
                    &rcond, ferr, berr, work.data(), rwork.data());
                GRVY_TIMER_BEGIN("suzerain_bsmbsm_zaPxpby");
                suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1, x.data(), 1, 0, p, 1);
                GRVY_TIMER_END("suzerain_bsmbsm_zaPxpby");
                break;
            case gbsv:
                info = suzerain_lapack_zgbsv(A.N, A.KL, A.KU, 1,
                    lu.data(), lu.colStride(), ipiv.data(), b.data(), A.N);
                GRVY_TIMER_BEGIN("suzerain_bsmbsm_zaPxpby");
                suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1, b.data(), 1, 0, p, 1);
                GRVY_TIMER_END("suzerain_bsmbsm_zaPxpby");
                break;
            }
            GRVY_TIMER_END("implicit operator solve");

            char buffer[128];
            if (info == 0) {
                // Success
            } else if (info < 0) {
                snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_%s reported error in argument %d",
                    method, -info);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            } else if (info <= A.N) {
                snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_%s reported singularity in PAP^T row %d"
                    " corresponding to A row %d for state scalar %d",
                    method, info - 1, suzerain_bsmbsm_q(A.S, A.n, info-1),
                    suzerain_bsmbsm_q(A.S, A.n, info-1) / A.n);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            } else if (info == A.N+1) {
                snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_%s reported condition number like %g for "
                    " m=%d, n=%d with km=%g, kn=%g",
                    method, 1/rcond, m, n, km, kn);
                WARN(buffer); // Warn user but continue...
            } else {
                snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_%s reported unknown error %d",
                    method, info);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            }

            if (solve_type == gbsvx) {
                // TODO Track and statistics on rcond, equed, ferr, and berr
            }
        }
    }

    // State leaves method as coefficients in X, Y, and Z directions

    GRVY_TIMER_END("invertMassPlusScaledOperator");
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
