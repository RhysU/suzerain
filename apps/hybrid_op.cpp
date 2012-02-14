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
#include <suzerain/complex.hpp>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#include "logging.hpp"

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
            for (int rowndx = begin; rowndx < end; ++rowndx) {
                if (rowndx == rho[wall]) {
                    // NOP
                } else if (rowndx == rhoe[wall]) {
                    row[rowndx*inc] = row[rho[wall]*inc]
                                    * gamma_times_one_minus_gamma;
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
    // Shorthand
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

    // Iterate across local wavenumbers and "invert" operator "in-place".
    // Does not shortcircuit on only-dealiased state (TODO should it?)
    for (int n = dkbz; n < dkez; ++n) {
        const real_t kn = twopioverLz*wavenumber(dNz, n);

        for (int m = dkbx; m < dkex; ++m) {
            const real_t km = twopioverLx*wavenumber(dNx, m);

            // Form complex-valued, wavenumber-dependent PAP^T within papt
            suzerain_rholut_imexop_packc(
                    phi, km, kn, &s, &ref, &ld, bop.get(),
                    ndx::rho, ndx::rhou, ndx::rhov, ndx::rhow, ndx::rhoe,
                    buf.data(), &A, papt.data());

            // Debug: determine what submatrices are contributing NaNs!
#ifndef NDEBUG
            {
                bool papt_contained_no_NaN_entries = true;
                for (int i = 0; i < A.N; ++i) {
                    const int qi = suzerain_bsmbsm_q(A.S, A.n, i);
                    for (int j = 0; j < A.N; ++j) {
                        const int qj = suzerain_bsmbsm_q(A.S, A.n, j);
                        if (suzerain_gbmatrix_in_band(A.LD,A.KL,A.KU,i,j)) {
                            int o = suzerain_gbmatrix_offset(A.LD,A.KL,A.KU,i,j);
                            if (papt(o) != papt(o)) { // isnan
                                WARN("NaN PAP^T_{"<<i<<","<<j<<"} from "
                                     <<"submatrix ("<<qi/A.n<<","<<qj/A.n<<") "
                                     <<"element ("<<qi%A.n<<","<<qj%A.n<<")");
                                papt_contained_no_NaN_entries = false;
                            }
                        }
                    }
                }
                assert(papt_contained_no_NaN_entries);
            }
#endif

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Given state pencil "p" the rest of the solve loop looks like
            //
            //     b := P p          using suzerain_bsmbsm_?aPxpby
            //     apply BC to RHS   using IsothermalNoSlipPAPTEnforcer
            //     apply BC to PAP^T using IsothermalNoSlipPAPTEnforcer
            //     x := (LU)^-1 b    using ?gbsvx which factorizes PAP^T
            //     p := P^T x        using suzerain_bsmbsm_?aPxpby
            //
            // where (km != 0 && kn != 0) uses complex-valued logic with one
            // right hand side.  The zero-zero modes require different handling
            // to prevent coupling between the real-valued state and the
            // imaginary-valued "state" used for integral constraints.  In some
            // cases those two could be solved together, but not here.
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
            // conditioned operators stemming (from reference values during
            // transients and high wavenumbers).  So we'll run with GBSVX until
            // it's demonstrably a bad idea.

            const char *method;               // Used for error reporting
            static const char fact    = 'E';  // Common inputs for {z,d}gbsvx
            static const char trans   = 'N';
            int info;                         // Common outputs for {z,d}gbsvx
            char equed;
            real_t rcond, ferr[2], berr[2];

            if (km || kn) {  // Complex-valued solve with one right hand side

                method = "zgbsvx";
                suzerain_bsmbsm_zaPxpby(
                        'N', A.S, A.n, 1., p, 1, 0., b.data(), 1);
                bc_enforcer.rhs(b.data());
                bc_enforcer.op(A, papt.data(), papt.colStride());
                info = suzerain_lapack_zgbsvx(fact, trans, A.N, A.KL, A.KU, 1,
                    papt.data(), papt.colStride(), lu.data(), lu.colStride(),
                    ipiv.data(), &equed, r.data(), c.data(),
                    b.data(), b.size(), x.data(), x.size(),
                    &rcond, ferr, berr, work.data(), rwork.data());
                suzerain_bsmbsm_zaPxpby(
                        'T', A.S, A.n, 1., x.data(), 1, 0., p, 1);

            } else {         // Two real-valued solves for zero-zero mode

                // Process one right hand side for each of real(p) and imag(p).
                //
                // The "Complex" operator is actually real-valued, we repack
                // PAP^T as a real operator for the solution process.  This
                // allows two real RHS in the same storage as one complex RHS;
                // convince yourself that DGBSVX can use ZGBSVX working storage.
                method = "dgbsvx";
                suzerain_bsmbsm_daPxpby(
                        'N', A.S, A.n, 1., ((real_t*)p),            2,
                                       0., ((real_t*)b.data()),     1);
                suzerain_bsmbsm_daPxpby(
                        'N', A.S, A.n, 1., ((real_t*)p)+1,          2,
                                       0., ((real_t*)b.data())+A.N, 1);
                bc_enforcer.rhs(((real_t*)b.data()));
                bc_enforcer.rhs(((real_t*)b.data())+A.N);
                for (int i = 0; i < papt.size(); ++i) {  // Pack real operator
                    ((real_t*)papt.data())[i] = papt(i).real();
                }
                bc_enforcer.op(A, (real_t*)papt.data(), papt.colStride());
                info = suzerain_lapack_dgbsvx(fact, trans, A.N, A.KL, A.KU, 2,
                            (real_t*)papt.data(), papt.colStride(),
                            (real_t*)lu.data(), lu.colStride(),
                            ipiv.data(), &equed, r.data(), c.data(),
                            (real_t*)b.data(), b.size(),
                            (real_t*)x.data(), x.size(), &rcond, ferr, berr,
                            (real_t*)work.data(), (int*)rwork.data());
                suzerain_bsmbsm_daPxpby(
                        'T', A.S, A.n, 1., ((real_t*)x.data()),     1,
                                       0., ((real_t*)p),            2);
                suzerain_bsmbsm_daPxpby(
                        'T', A.S, A.n, 1., ((real_t*)x.data())+A.N, 1,
                                       0., ((real_t*)p)+1,          2);
                ferr[0] = std::max(ferr[0], ferr[1]);
                berr[0] = std::max(berr[0], berr[1]);
            }

            if (info) {
                char buffer[80];
                snprintf(buffer, sizeof(buffer),
                         "suzerain_lapack_%s reported error %d", method, info);
                SUZERAIN_ERROR_VOID(buffer, SUZERAIN_ESANITY);
            }

            // TODO Track and report statistics on rcond, equed, ferr, and berr
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
