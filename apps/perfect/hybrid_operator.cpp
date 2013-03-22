//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc hybrid_operator.hpp
 */

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
#include <suzerain/ndx.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/timers.h>

#include "hybrid_operator.hpp"
#include "nonlinear_operator.hpp"

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

isothermal_hybrid_linear_operator::isothermal_hybrid_linear_operator(
        const zgbsv_specification& spec,
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : operator_base(grid, dgrid, cop, b)
    , solver(bsmbsm_solver::build(suzerain_bsmbsm_construct(
                5, dgrid.global_wave_extent.y(), cop.max_kl(), cop.max_ku()),
                spec, 1))
    , scenario(scenario)
    , common(common)
    , who("operator.L")
{
    // NOP
}

isothermal_hybrid_linear_operator::~isothermal_hybrid_linear_operator()
{
    // TODO Allreduce bsmbsm_solver statistics from all ranks
    if (solver) {
        std::vector<std::string> summary = solver->summarize_statistics();
        for (std::size_t i = 0; i < summary.size(); ++i) {
            INFO0(who, summary[i]);
        }
    }
}

void isothermal_hybrid_linear_operator::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("apply_mass_plus_scaled_operator");
    SUZERAIN_ENSURE(common.linearization != linearize::none);

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < solver->S);
    assert(static_cast<int>(ndx::mx ) < solver->S);
    assert(static_cast<int>(ndx::my ) < solver->S);
    assert(static_cast<int>(ndx::mz ) < solver->S);
    assert(static_cast<int>(ndx::rho) < solver->S);

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
    SUZERAIN_ENSURE(state.shape()  [0] == (unsigned) solver->S);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned)        Ny);
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned)        Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==                    1);

    // Scratch for "in-place" suzerain_rholut_imexop_accumulate usage
    VectorXc tmp(solver->N);
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
            blas::copy(solver->N, p, 1, tmp.data(), 1);

            // Accumulate result back into state storage automatically
            // adjusting for when input imaginary part a priori should be zero
            SUZERAIN_TIMER_SCOPED("suzerain_rholut_imexop_accumulate");
            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, cop.get(),
                    tmp.data() + ndx::e   * Ny,
                    tmp.data() + ndx::mx  * Ny,
                    tmp.data() + ndx::my  * Ny,
                    tmp.data() + ndx::mz  * Ny,
                    tmp.data() + ndx::rho * Ny,
                    0.0,
                    p + ndx::e   * Ny,
                    p + ndx::mx  * Ny,
                    p + ndx::my  * Ny,
                    p + ndx::mz  * Ny,
                    p + ndx::rho * Ny);
        }
    }
}

void isothermal_hybrid_linear_operator::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("accumulate_mass_plus_scaled_operator");
    SUZERAIN_ENSURE(common.linearization != linearize::none);

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < solver->S);
    assert(static_cast<int>(ndx::mx ) < solver->S);
    assert(static_cast<int>(ndx::my ) < solver->S);
    assert(static_cast<int>(ndx::mz ) < solver->S);
    assert(static_cast<int>(ndx::rho) < solver->S);

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
    SUZERAIN_ENSURE(output.is_isomorphic(input));
    SUZERAIN_ENSURE(input.shape()   [0] == (unsigned) solver->S);
    SUZERAIN_ENSURE(input.strides() [0] == (unsigned)        Ny);
    SUZERAIN_ENSURE(input.shape()   [1] == (unsigned)        Ny);
    SUZERAIN_ENSURE(input.strides() [1] ==                    1);
    SUZERAIN_ENSURE(output.strides()[1] ==                    1);

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
            SUZERAIN_TIMER_SCOPED("suzerain_rholut_imexop_accumulate");
            suzerain_rholut_imexop_accumulate(
                    phi, km, kn, &s, &ref, &ld, cop.get(),
                    &input [ndx::e  ][0][m - dkbx][n - dkbz],
                    &input [ndx::mx ][0][m - dkbx][n - dkbz],
                    &input [ndx::my ][0][m - dkbx][n - dkbz],
                    &input [ndx::mz ][0][m - dkbx][n - dkbz],
                    &input [ndx::rho][0][m - dkbx][n - dkbz],
                    beta,
                    &output[ndx::e   ][0][m - dkbx][n - dkbz],
                    &output[ndx::mx  ][0][m - dkbx][n - dkbz],
                    &output[ndx::my  ][0][m - dkbx][n - dkbz],
                    &output[ndx::mz  ][0][m - dkbx][n - dkbz],
                    &output[ndx::rho ][0][m - dkbx][n - dkbz]);

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

    IsothermalNoSlipPATPTEnforcer(const suzerain_bsmbsm &A_T,
                                  const suzerain_rholut_imexop_scenario &s)
        : gamma_times_one_minus_gamma(s.gamma * (1 - s.gamma))
    {
        // Starting offset to named scalars in interleaved_state pencil
        const int e0   = static_cast<int>(ndx::e  ) * A_T.n;
        const int mx0  = static_cast<int>(ndx::mx ) * A_T.n;
        const int my0  = static_cast<int>(ndx::my ) * A_T.n;
        const int mz0  = static_cast<int>(ndx::mz ) * A_T.n;
        const int rho0 = static_cast<int>(ndx::rho) * A_T.n;

        // Relative to foo0 what is the offset to lower, upper walls
        const int wall[nwalls] = { 0, A_T.n - 1};

        // Prepare indices within PA^TP^T corresponding to the walls.
        // Uses that {A^T}_{i,j} maps to {PA^TP^T}_{{q^-1}(i),{q^(-1)}(j)}.
        for (int i = 0; i < nwalls; ++i) {
            e     [i]    = suzerain_bsmbsm_qinv(A_T.S, A_T.n, e0   + wall[i]);
            noslip[i][0] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, mx0  + wall[i]);
            noslip[i][1] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, my0  + wall[i]);
            noslip[i][2] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, mz0  + wall[i]);
            rho   [i]    = suzerain_bsmbsm_qinv(A_T.S, A_T.n, rho0 + wall[i]);
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

void isothermal_hybrid_linear_operator::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("invert_mass_plus_scaled_operator");
    SUZERAIN_ENSURE(common.linearization != linearize::none);

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < solver->S);
    assert(static_cast<int>(ndx::mx ) < solver->S);
    assert(static_cast<int>(ndx::my ) < solver->S);
    assert(static_cast<int>(ndx::mz ) < solver->S);
    assert(static_cast<int>(ndx::rho) < solver->S);

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
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned)        Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==                    1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned)        Ny);
    SUZERAIN_ENSURE(state.shape()  [0] == (unsigned) solver->S);

    // Compute how many additional mean constraints we must solve
    // Ensure conformant, mean constraints are arriving on the correct rank
    const std::size_t nconstraints = ic0 ? ic0->shape()[2]*ic0->shape()[3] : 0;
    if (nconstraints) {
        SUZERAIN_ENSURE(dgrid.has_zero_zero_modes());
        SUZERAIN_ENSURE(ic0->shape()  [1] == (unsigned)        Ny);
        SUZERAIN_ENSURE(ic0->strides()[1] ==                    1);
        SUZERAIN_ENSURE(ic0->strides()[0] == (unsigned)        Ny);
        SUZERAIN_ENSURE(ic0->shape()  [0] == (unsigned) solver->S);
    }

    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions

    // Pack reference details for suzerain_rholut_imexop routines
    suzerain_rholut_imexop_scenario s(this->imexop_s());
    suzerain_rholut_imexop_ref   ref;
    suzerain_rholut_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Prepare an almost functor mutating RHS and PA^TP^T to enforce BCs.
    IsothermalNoSlipPATPTEnforcer bc_enforcer(*solver, s);

    // Prepare a scratch buffer for packc/packf usage
    ArrayXXc buf(solver->ld, solver->n);

    // Iterate across local wavenumbers and "invert" operator "in-place"
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        const real_t kn = twopioverLz*wn;

        // Factorization reuse will not aid us across large jumps in km
        if (dkex - dkbx > 1) solver->apprx(false);

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            const real_t km = twopioverLx*wm;

            // Get pointer to (.,m,n)-th state pencil
            complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

            // Short circuiting didn't yet occur for Nyquist/dealiasing modes...
            if (   std::abs(wn) > wavenumber_absmin(Nz)
                || std::abs(wm) > wavenumber_absmin(Nx)) {
                memset(p, 0, solver->N*sizeof(p[0]));  // ...so we can zero,
                solver->apprx(false);                  // mark reuse moot, and
                continue;                              // then short circuit.
            }

            // Form complex-valued, wavenumber-dependent PA^TP^T
            static const char trans = 'T';
            if (solver->spec.in_place()) { // Pack for in-place LU
                SUZERAIN_TIMER_SCOPED("implicit operator assembly (packf)");
                suzerain_rholut_imexop_packf(
                        phi, km, kn, &s, &ref, &ld, cop.get(),
                        ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                        buf.data(), solver.get(), solver->LU.data());
            } else {                       // Pack for out-of-place LU
                SUZERAIN_TIMER_SCOPED("implicit operator assembly (packc)");
                suzerain_rholut_imexop_packc(
                        phi, km, kn, &s, &ref, &ld, cop.get(),
                        ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                        buf.data(), solver.get(), solver->PAPT.data());
            }
            // Apply boundary conditions to PA^TP^T
            {
                SUZERAIN_TIMER_SCOPED("implicit operator BCs");
                bc_enforcer.op(*solver, solver->PAPT.data(),
                                        solver->PAPT.colStride());
            }
            // Inform the solver about the new, unfactorized operator
            solver->supplied_PAPT();

            // Form right hand side, apply BCs, factorize, and solve.
            // Beware that much ugly, ugly magic is hidden just below.
            solver->supply_B(p);
            {
                SUZERAIN_TIMER_SCOPED("implicit right hand side BCs");
                bc_enforcer.rhs(solver->PB.data());
            }
            solver->solve(trans);
            solver->demand_X(p);

            // If necessary, solve any required integral constraints
            if (SUZERAIN_UNLIKELY(n == 0 && m == 0)) {
                for (std::size_t i = 0; i < nconstraints; ++i) {
                    SUZERAIN_TIMER_SCOPED("implicit constraint solution");
                    solver->supply_B(ic0->data() + i * solver->N);
                    bc_enforcer.rhs(solver->PB.data());
                    solver->solve(trans);
                    solver->demand_X(ic0->data() + i * solver->N);
                }
            }
        }
    }

    // State leaves method as coefficients in X, Y, and Z directions
}

suzerain_rholut_imexop_scenario
isothermal_hybrid_linear_operator::imexop_s() const
{
    suzerain_rholut_imexop_scenario retval;
    retval.Re    = scenario.Re;
    retval.Pr    = scenario.Pr;
    retval.Ma    = scenario.Ma;
    retval.alpha = scenario.alpha;
    retval.gamma = scenario.gamma;
    return retval;
}

hybrid_nonlinear_operator::hybrid_nonlinear_operator(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , msoln(msoln)
    , who("operator.N")
{
    // NOP
}

std::vector<real_t> hybrid_nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return apply_navier_stokes_spatial_operator<
                true,  linearize::rhome_xyz, slowgrowth::none
            >(this->scenario.alpha,
              this->scenario.beta,
              this->scenario.gamma,
              this->scenario.Ma,
              this->scenario.Pr,
              this->scenario.Re,
              *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return apply_navier_stokes_spatial_operator<
                false, linearize::rhome_xyz, slowgrowth::none
            >(this->scenario.alpha,
              this->scenario.beta,
              this->scenario.gamma,
              this->scenario.Ma,
              this->scenario.Pr,
              this->scenario.Re,
              *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} /* namespace perfect */ } /* namespace suzerain */
