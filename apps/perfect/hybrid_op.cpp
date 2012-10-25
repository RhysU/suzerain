//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// explicit_op.hpp: Operators for channel simulation
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/complex.hpp>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timers.h>

#include "../logging.hpp"
#include "hybrid_op.hpp"
#include "nonlinear.hpp"

#pragma warning(disable:383 1572)

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
suzerain::real_t twopiover(const suzerain::real_t L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

namespace suzerain { namespace perfect {

void HybridIsothermalLinearOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("applyMassPlusScaledOperator");

    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    namespace field = support::field;
    SUZERAIN_UNUSED(method);
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
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);

    // Scratch for "in-place" suzerain_rholut_imexop_accumulate usage
    VectorXc tmp(Ny*field::count);
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Iterate across local wavenumbers and apply operator "in-place"
    // Short circuit "continues" occur for Nyquist and non-dealiased modes...
    // ...where the former will be zeroed during the later invert call.
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        if (std::abs(wn) > wavenumber_absmin(Nz)) continue;
        const real_t kn = twopioverLz*wn;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            if (std::abs(wm) > wavenumber_absmin(Nx)) continue;
            const real_t km = twopioverLx*wm;

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Copy pencil into temporary storage
            blas::copy(field::count*Ny, p, 1, tmp.data(), 1);

            // Accumulate result back into state storage automatically
            // adjusting for when input imaginary part a priori should be zero
            SUZERAIN_TIMER_BEGIN("suzerain_rholut_imexop_accumulate");
            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    wn == 0 && wm == 0,
                    tmp.data() + field::ndx::e   * Ny,
                    tmp.data() + field::ndx::mx  * Ny,
                    tmp.data() + field::ndx::my  * Ny,
                    tmp.data() + field::ndx::mz  * Ny,
                    tmp.data() + field::ndx::rho * Ny,
                    0.0,
                    p + field::ndx::e   * Ny,
                    p + field::ndx::mx  * Ny,
                    p + field::ndx::my  * Ny,
                    p + field::ndx::mz  * Ny,
                    p + field::ndx::rho * Ny);
            SUZERAIN_TIMER_END("suzerain_rholut_imexop_accumulate");
        }
    }
}

void HybridIsothermalLinearOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        ContiguousState<4,complex_t> &output,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("accumulateMassPlusScaledOperator");

    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    namespace field = support::field;
    SUZERAIN_UNUSED(method);
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
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

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

    // Iterate across local wavenumbers and apply operator "in-place"
    // Short circuit "continues" occur for Nyquist and non-dealiased modes...
    // ...where the former will be zeroed during the later invert call.
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        if (std::abs(wn) > wavenumber_absmin(Nz)) continue;
        const real_t kn = twopioverLz*wn;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            if (std::abs(wm) > wavenumber_absmin(Nx)) continue;
            const real_t km = twopioverLx*wm;

            // Accumulate result automatically adjusting for when input
            // imaginary part a priori should be zero
            SUZERAIN_TIMER_BEGIN("suzerain_rholut_imexop_accumulate");
            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    wn == 0 && wm == 0,
                    &input [field::ndx::e  ][0][m - dkbx][n - dkbz],
                    &input [field::ndx::mx ][0][m - dkbx][n - dkbz],
                    &input [field::ndx::my ][0][m - dkbx][n - dkbz],
                    &input [field::ndx::mz ][0][m - dkbx][n - dkbz],
                    &input [field::ndx::rho][0][m - dkbx][n - dkbz],
                    beta,
                    &output[field::ndx::e   ][0][m - dkbx][n - dkbz],
                    &output[field::ndx::mx  ][0][m - dkbx][n - dkbz],
                    &output[field::ndx::my  ][0][m - dkbx][n - dkbz],
                    &output[field::ndx::mz  ][0][m - dkbx][n - dkbz],
                    &output[field::ndx::rho ][0][m - dkbx][n - dkbz]);
            SUZERAIN_TIMER_END("suzerain_rholut_imexop_accumulate");

        }
    }
}

/**
 * A functor for applying isothermal conditions for PA^TP^T-based operators.
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
class IsothermalNoSlipPATPTEnforcer
{
    enum { nwalls = 2, nmomentum = 3 };

    // Indices within PA^TP^T at which to apply boundary conditions
    // Computed once within constructor and then repeatedly used
    int rho[nwalls], noslip[nwalls][nmomentum], e[nwalls];

    // Precomputed coefficient based on the isothermal equation of state
    real_t gamma_times_one_minus_gamma;

public:

    IsothermalNoSlipPATPTEnforcer(const suzerain_bsmbsm &A,
                                  const suzerain_rholut_imexop_scenario &s)
        : gamma_times_one_minus_gamma(s.gamma * (1 - s.gamma))
    {
        // Starting offset to named scalars in InterleavedState pencil
        namespace ndx = support::field::ndx;
        const int start_e   = static_cast<int>(ndx::e  ) * A.n;
        const int start_mx  = static_cast<int>(ndx::mx ) * A.n;
        const int start_my  = static_cast<int>(ndx::my ) * A.n;
        const int start_mz  = static_cast<int>(ndx::mz ) * A.n;
        const int start_rho = static_cast<int>(ndx::rho) * A.n;

        // Relative to start_foo what is the offset to lower, upper walls
        const int wall[nwalls] = { 0, A.n - 1};

        // Prepare indices within PA^TP^T corresponding to the walls.
        // Uses that A_{i,j} maps to {PA^TP^T}_{{q^-1}(i),{q^(-1)}(j)}.
        for (int i = 0; i < nwalls; ++i) {
            e     [i]    = suzerain_bsmbsm_qinv(A.S, A.n, start_e   + wall[i]);
            noslip[i][0] = suzerain_bsmbsm_qinv(A.S, A.n, start_mx  + wall[i]);
            noslip[i][1] = suzerain_bsmbsm_qinv(A.S, A.n, start_my  + wall[i]);
            noslip[i][2] = suzerain_bsmbsm_qinv(A.S, A.n, start_mz  + wall[i]);
            rho   [i]    = suzerain_bsmbsm_qinv(A.S, A.n, start_rho + wall[i]);
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
            b[e[wall]] = 0;
        }
    }

    /**
     * Modify the equations within PA^TP^T for lower, upper walls where each
     * contiguous column within PA^TP^T contains one equation.  This storage is
     * ideal from a cache locality perspective for this boundary condition.
     * Must be done in conjunction with rhs().
     */
    template<typename T>
    void op(const suzerain_bsmbsm& A_T, T * const patpt, int patpt_ld)
    {
        // Attempt made to not unnecessarily disturb matrix conditioning.

        for (int wall = 0; wall < nwalls; ++wall) {
            int begin, end;

            // Zero off-diagonals so momentum equations are const*rho{u,v,w}=0.
            for (size_t eqn = 0; eqn < nmomentum; ++eqn) {
                T * const col = (T *) suzerain_gbmatrix_col(
                        A_T.N, A_T.N, A_T.KL, A_T.KU, (void *) patpt, patpt_ld,
                        sizeof(T), noslip[wall][eqn], &begin, &end);
                for (int i = begin; i < end; ++i) {
                    if (i != noslip[wall][eqn]) {
                        col[i] = 0;
                    }
                }
            }

            // Set constraint const*rho - const*gamma*(gamma-1)*rhoE = 0
            T * const col = (T *) suzerain_gbmatrix_col(
                    A_T.N, A_T.N, A_T.KL, A_T.KU, (void *) patpt, patpt_ld,
                    sizeof(T), e[wall], &begin, &end);
            // Necessary to ensure constraint possible in degenerate case
            complex_t &rhocoeff = col[rho[wall]];
            if (rhocoeff == complex_t(0)) rhocoeff = 1;
            // Scan row and adjust coefficients for constraint
            for (int i = begin; i < end; ++i) {
                if (i == e[wall]) {
                    col[i] = rhocoeff*gamma_times_one_minus_gamma;
                } else if (i == rho[wall]) {
                    // NOP
                } else {
                    col[i] = 0;
                }
            }
        }
    }
};

void HybridIsothermalLinearOperator::invertMassPlusScaledOperator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> * const ic0) const
{
    SUZERAIN_TIMER_BEGIN("invertMassPlusScaledOperator");

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    namespace field = support::field;
    namespace ndx   = field::ndx;
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

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
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

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
    SCRATCH_C(       ArrayXXc, buf,     A.ld,      A.n); // For packc calls
    SCRATCH_C(       ArrayXXc, patpt,   A.LD,      A.N); // Holds PA^TP^T
    SCRATCH_C(       ArrayXXc, lu,      A.LD+A.KL, A.N); // Holds LU of PA^TP^T
    SCRATCH_I(Eigen::ArrayXi,  ipiv,    A.N);            // Linear solve...
    SCRATCH_R(       ArrayXr,  r,       A.N);
    SCRATCH_R(       ArrayXr,  c,       A.N);
    SCRATCH_C(       ArrayXc,  b,       A.N);
    SCRATCH_C(       ArrayXc,  x,       A.N);
    SCRATCH_C(       ArrayXc,  work,  2*A.N);
    SCRATCH_R(       ArrayXr,  rwork,   A.N);
#undef SCRATCH_C
#undef SCRATCH_R
#undef SCRATCH_I

    // Pack reference details for suzerain_rholut_imexop routines
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Prepare an almost functor mutating RHS and PA^TP^T to enforce BCs.
    IsothermalNoSlipPATPTEnforcer bc_enforcer(A, s);

    // Solver-related operational details
    const char *mname = "UNKNOWN";          // Used for error reporting
    char fact = spec.equil() ? 'E' : 'N';   // Equilibrate?
    switch (spec.method()) {
    case spec_zgbsv::zgbsv:   mname = "suzerain_lapackext_zgbsv";   break;
    case spec_zgbsv::zgbsvx:  mname = "suzerain_lapack_zgbsvx";     break;
    case spec_zgbsv::zcgbsvx: mname = "suzerain_lapackext_zcgbsvx"; break;
    }
    static const char trans = 'T';  // Un-transpose transposed operator
    int info;                       // Common output for all solvers
    char equed;                     // zgbsvx equilibration type
    real_t rcond, ferr, berr;       // zgbsvx outputs for one RHS
    real_t afrob, tolsc, res;       // zcgbsvx outputs for one RHS...
    int apprx, aiter, siter, diter; // ...Ditto

    // Iterate across local wavenumbers and "invert" operator "in-place"
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        const real_t kn = twopioverLz*wn;

        // Factorization reuse will not aid us across large jumps in km
        fact = (spec.reuse() && dkex - dkbx > 1) ? 'N' : fact;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            const real_t km = twopioverLx*wm;

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Short circuiting did NOT occur for Nyquist or dealiasing modes...
            if (   std::abs(wn) > wavenumber_absmin(Nz)
                || std::abs(wm) > wavenumber_absmin(Nx)) {
                memset(p, 0, A.N*sizeof(p[0]));             // so we may zero,
                fact = spec.reuse() ? 'N' : fact;           // mark reuse moot,
                continue;                                   // and then bail.
            }

            SUZERAIN_TIMER_BEGIN("implicit operator miscellaneous");

            // Form complex-valued, wavenumber-dependent PA^TP^T within patpt.
            // This is the transpose of the implicit operator we desire.
            SUZERAIN_TIMER_BEGIN("implicit operator assembly");
            if (spec.method() == spec_zgbsv::zgbsv) {  // In-place LUP
                suzerain_rholut_imexop_packf(
                        phi, km, kn, &s, &ref, &ld, bop.get(),
                        ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                        buf.data(), &A, lu.data());
            } else {                                   // Out-of-place LUP
                suzerain_rholut_imexop_packc(
                        phi, km, kn, &s, &ref, &ld, bop.get(),
                        ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                        buf.data(), &A, patpt.data());
            }
            SUZERAIN_TIMER_END("implicit operator assembly");

            // Given state pencil "p" the rest of the solve loop looks like
            //
            //     b := P p            using suzerain_bsmbsm_?aPxpby
            //     apply BC to RHS     using IsothermalNoSlipPATPTEnforcer
            //     apply BC to PA^TP^T using IsothermalNoSlipPATPTEnforcer
            //     x := (LU)^-T b      using ?gbsvx which factorizes PA^TP^T
            //                         or ?gbsv which factorizes in place
            //     p := P^T x          using suzerain_bsmbsm_?aPxpby

            SUZERAIN_TIMER_BEGIN("suzerain_bsmbsm_zaPxpby");
            suzerain_bsmbsm_zaPxpby('N', A.S, A.n, 1, p, 1, 0, b.data(), 1);
            SUZERAIN_TIMER_END("suzerain_bsmbsm_zaPxpby");

            SUZERAIN_TIMER_BEGIN("implicit operator BCs");
            bc_enforcer.rhs(b.data());
            if (spec.method() == spec_zgbsv::zgbsv) {  // In-place LUP
                bc_enforcer.op(A, lu.data() + A.KL, lu.colStride());
            } else {                                   // Out-of-place LUP
                bc_enforcer.op(A, patpt.data(), patpt.colStride());
            }
            SUZERAIN_TIMER_END("implicit operator BCs");

            switch (spec.method()) {
            default:
                SUZERAIN_ERROR_VOID("unknown solve_type", SUZERAIN_ESANITY);

            case spec_zgbsv::zgbsv:                    // In-place LUP, solve
                SUZERAIN_TIMER_BEGIN(mname);
                info = suzerain_lapackext_zgbsv(trans, A.N, A.KL, A.KU, 1,
                    lu.data(), lu.colStride(), ipiv.data(), b.data(), A.N);
                SUZERAIN_TIMER_END(mname);
                break;

            case spec_zgbsv::zgbsvx:                   // Out-of-place
                SUZERAIN_TIMER_BEGIN(mname);
                info = suzerain_lapack_zgbsvx(fact, trans, A.N, A.KL, A.KU, 1,
                    patpt.data(), patpt.colStride(), lu.data(), lu.colStride(),
                    ipiv.data(), &equed, r.data(), c.data(),
                    b.data(), A.N, x.data(), A.N,
                    &rcond, &ferr, &berr, work.data(), rwork.data());
                SUZERAIN_TIMER_END(mname);
                // TODO Statistics on rcond, equed, ferr, and berr
                break;

            case spec_zgbsv::zcgbsvx:                  // Out-of-place
                SUZERAIN_TIMER_BEGIN(mname);
                fact  = spec.reuse() ? fact : 'N';
                apprx = fact == 'N' ? 0 : 1;
                aiter = spec.aiter();
                afrob = -1;
                siter = spec.siter();
                diter = spec.diter();
                tolsc = spec.tolsc();
                info  = suzerain_lapackext_zcgbsvx(&fact, &apprx, aiter, trans,
                        A.N, A.KL, A.KU, patpt.data(), &afrob, lu.data(),
                        ipiv.data(), b.data(), x.data(), &siter, &diter,
                        &tolsc, work.data(), &res);
                SUZERAIN_TIMER_END(mname);
                // TODO Statistics on fact, apprx, siter, diter, tolsc, res
                break;
            }

            SUZERAIN_TIMER_BEGIN("suzerain_bsmbsm_zaPxpby");
            if (spec.method() == spec_zgbsv::zgbsv) {  // In-place solve
                suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1, b.data(), 1, 0, p, 1);
            } else {                                   // Out-of-place solve
                suzerain_bsmbsm_zaPxpby('T', A.S, A.n, 1, x.data(), 1, 0, p, 1);
            }
            SUZERAIN_TIMER_END("suzerain_bsmbsm_zaPxpby");

            // Report any errors that occurred during the solve
            char buffer[128];
            if (info == 0) {
                // Success
            } else if (info < 0) {
                snprintf(buffer, sizeof(buffer),
                    "%s reported error in argument %d",
                    mname, -info);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            } else if (info <= A.N) {
                snprintf(buffer, sizeof(buffer),
                    "%s reported singularity in PA^TP^T row %d"
                    " corresponding to A row %d for state scalar %d",
                    mname, info - 1, suzerain_bsmbsm_q(A.S, A.n, info-1),
                    suzerain_bsmbsm_q(A.S, A.n, info-1) / A.n);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            } else if (info == A.N+1 && spec.method() == spec_zgbsv::zgbsvx) {
                snprintf(buffer, sizeof(buffer),
                    "%s reported condition number like %g for "
                    " m=%d, n=%d with km=%g, kn=%g",
                    mname, 1/rcond, m, n, km, kn);
                WARN(buffer); // Warn user but continue...
            } else {
                snprintf(buffer, sizeof(buffer),
                    "%s reported unknown error %d", mname, info);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            }

            SUZERAIN_TIMER_END("implicit operator miscellaneous");
        }
    }

    // State leaves method as coefficients in X, Y, and Z directions

    SUZERAIN_TIMER_END("invertMassPlusScaledOperator");
}

std::vector<real_t> HybridNonlinearOperator::applyOperator(
            const real_t time,
            ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return applyNonlinearOperator<true,  linearize::rhome>
            (this->scenario.alpha,
             this->scenario.beta,
             this->scenario.gamma,
             this->scenario.Ma,
             this->scenario.Pr,
             this->scenario.Re,
             *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return applyNonlinearOperator<false, linearize::rhome>
            (this->scenario.alpha,
             this->scenario.beta,
             this->scenario.gamma,
             this->scenario.Ma,
             this->scenario.Pr,
             this->scenario.Re,
             *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} /* namespace perfect */ } /* namespace suzerain */
