//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc hybrid_operator.hpp
 */

#include "hybrid_operator.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/bsmbsm_solver.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/inorder.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/specification_isothermal.hpp>
#include <suzerain/specification_zgbsv.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/timers.h>

#include "nonlinear_operator_fwd.hpp"
#include "antioch_constitutive.hpp"
#include "definition_channel.hpp"
#include "reacting_ndx.hpp"

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

namespace suzerain {

namespace reacting {

isothermal_hybrid_linear_operator::isothermal_hybrid_linear_operator(
        const specification_zgbsv& spec,
        const antioch_constitutive &cmods,
        const specification_isothermal &isospec,
        const definition_channel &chdef,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : operator_base(grid, dgrid, cop, b)
    , flow_solver(bsmbsm_solver::build(suzerain_bsmbsm_construct(
                  5, dgrid.global_wave_extent.y(), cop.max_kl(), cop.max_ku()),
                  spec, 1))
    , isospec(isospec)
    , chdef(chdef)
    , cmods(cmods)
    , common(common)
    , who("operator.L")
{
    // Prepare species solvers
    species_solver.reserve(cmods.Ns()-1);
    
    for (unsigned int i=0; i<cmods.Ns()-1; ++i) {
        shared_ptr<bsmbsm_solver> tmp(
            bsmbsm_solver::build(suzerain_bsmbsm_construct(
                                     1, dgrid.global_wave_extent.y(), 
                                     cop.max_kl(), cop.max_ku()), spec, 1));

        species_solver.push_back(tmp);
    }
    
}

isothermal_hybrid_linear_operator::~isothermal_hybrid_linear_operator()
{
    // TODO Allreduce bsmbsm_solver statistics from all ranks
    if (flow_solver) {
        std::vector<std::string> summary = flow_solver->summarize_statistics();
        for (std::size_t i = 0; i < summary.size(); ++i) {
            INFO0(who, summary[i]);
        }
    }

    for (std::size_t i = 0; i < species_solver.size(); ++i) {
        if (species_solver[i]) {
            std::vector<std::string> 
                summary = species_solver[i]->summarize_statistics();
            for (std::size_t i = 0; i < summary.size(); ++i) {
                INFO0(who, summary[i]);
            }
        }
    }
}

void isothermal_hybrid_linear_operator::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("apply_mass_plus_scaled_operator");
    SUZERAIN_ENSURE(common.linearization != linearize::none);

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    SUZERAIN_UNUSED(substep_index);

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < flow_solver->S);
    assert(static_cast<int>(ndx::mx ) < flow_solver->S);
    assert(static_cast<int>(ndx::my ) < flow_solver->S);
    assert(static_cast<int>(ndx::mz ) < flow_solver->S);
    assert(static_cast<int>(ndx::rho) < flow_solver->S);

    // TODO: Assert species solver sized correctly

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
    SUZERAIN_ENSURE(state.shape()  [0] == (unsigned) (flow_solver->S+
                                                      species_solver.size()));
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned)             Ny);
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned)             Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==                         1);

    // Compute total state size for pencil
    int Ntot = flow_solver->N;
    if (species_solver.size()>0) {
        // All species have same size
        Ntot += species_solver.size()*species_solver[0]->N;
    }

    // Scratch for "in-place" suzerain_reacting_imexop_accumulate usage
    VectorXc tmp(Ntot);
    suzerain_reacting_imexop_scenario s(this->imexop_s());
    suzerain_reacting_imexop_ref   ref;
    suzerain_reacting_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Vector to save top boundary state in case we need it to
    // overwrite action of linear operator (for used in boundary layer
    // case)
    VectorXc top_state(state.shape()[0]);

    // Iterate across local wavenumbers and apply operator "in-place"
    // Short circuit "continues" occur for Nyquist and non-dealiased modes...
    // ...where the former will be zeroed during the later invert call.
    switch (common.linearization) {

    case linearize::rhome_y:
    {
        for (int n = dkbz; n < dkez; ++n) {
            const int wn = wavenumber(dNz, n);
            if (std::abs(wn) > wavenumber_absmin(Nz)) continue;

            for (int m = dkbx; m < dkex; ++m) {
                const int wm = wavenumber(dNx, m);
                if (std::abs(wm) > wavenumber_absmin(Nx)) continue;

                // Get pointer to (.,m,n)-th state pencil
                complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

                // if BL, save top state
                if (grid.one_sided())
                    for (std::size_t i=0; i<state.shape()[0]; ++i)
                        top_state[i] = state[i][Ny-1][m - dkbx][n - dkbz];

                // Copy pencil into temporary storage
                blas::copy(Ntot, p, 1, tmp.data(), 1);

                // Accumulate result back into state storage
                SUZERAIN_TIMER_SCOPED("suzerain_reacting_imexop_accumulate");

                // Flow equations
                suzerain_reacting_flow_imexop_accumulate(
                        phi, &s, &ref, &ld, cop.get(), 
                        0, // imagzero=false FIXME: update after removing from reacting_imexop
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

                // Species equations
                for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
                    suzerain_reacting_species_imexop_accumulate(
                        phi, &s, &ref, &ld, cop.get(),
                        0, // imagzero=false FIXME: update after removing from reacting_imexop
                        tmp.data() + (ndx::species+alfa) * Ny,
                        0.0,
                        p + (ndx::species+alfa) * Ny);
                }

                // if BL, replace top state
                if (grid.one_sided())
                    for (std::size_t i=0; i<state.shape()[0]; ++i)
                        state[i][Ny-1][m - dkbx][n - dkbz] = top_state[i];

            }
        }
        break;
    }

    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    }
}

void isothermal_hybrid_linear_operator::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    SUZERAIN_TIMER_SCOPED("accumulate_mass_plus_scaled_operator");
    SUZERAIN_ENSURE(common.linearization != linearize::none);

    // Shorthand
    using inorder::wavenumber;
    using inorder::wavenumber_absmin;
    SUZERAIN_UNUSED(substep_index);

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < flow_solver->S);
    assert(static_cast<int>(ndx::mx ) < flow_solver->S);
    assert(static_cast<int>(ndx::my ) < flow_solver->S);
    assert(static_cast<int>(ndx::mz ) < flow_solver->S);
    assert(static_cast<int>(ndx::rho) < flow_solver->S);

    // TODO: Assert species solver sized correctly

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
    SUZERAIN_ENSURE(input.shape()   [0] == (unsigned) (flow_solver->S +
                                                       species_solver.size()));
    SUZERAIN_ENSURE(input.strides() [0] == (unsigned)             Ny);
    SUZERAIN_ENSURE(input.shape()   [1] == (unsigned)             Ny);
    SUZERAIN_ENSURE(input.strides() [1] ==                         1);
    SUZERAIN_ENSURE(output.strides()[1] ==                         1);

    // Scratch for suzerain_reacting_imexop_accumulate usage
    suzerain_reacting_imexop_scenario s(this->imexop_s());
    suzerain_reacting_imexop_ref   ref;
    suzerain_reacting_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Vector to save top boundary state in case we need it to
    // overwrite action of linear operator (for used in boundary layer
    // case)
    VectorXc top_input (input.shape() [0]);
    VectorXc top_output(output.shape()[0]);

    // Iterate across local wavenumbers and apply operator "in-place"
    // Short circuit "continues" occur for Nyquist and non-dealiased modes...
    // ...where the former will be zeroed during the later invert call.
    switch (common.linearization) {

    case linearize::rhome_y:
    {
        for (int n = dkbz; n < dkez; ++n) {
            const int wn = wavenumber(dNz, n);
            if (std::abs(wn) > wavenumber_absmin(Nz)) continue;

            for (int m = dkbx; m < dkex; ++m) {
                const int wm = wavenumber(dNx, m);
                if (std::abs(wm) > wavenumber_absmin(Nx)) continue;


                // if BL, save top state
                if (grid.one_sided()) {
                    for (std::size_t i=0; i<input.shape()[0]; ++i) {
                        top_input[i]  = input [i][Ny-1][m - dkbx][n - dkbz];
                        top_output[i] = output[i][Ny-1][m - dkbx][n - dkbz];
                    }
                }

                // Accumulate result
                SUZERAIN_TIMER_SCOPED("suzerain_reacting_imexop_accumulate");

                // Flow equations
                suzerain_reacting_flow_imexop_accumulate(
                        phi, &s, &ref, &ld, cop.get(),
                        0, // imagzero=false FIXME: update after removing from reacting_imexop
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

                // Species equations
                for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
                    suzerain_reacting_species_imexop_accumulate(
                        phi, &s, &ref, &ld, cop.get(),
                        0, // imagzero=false FIXME: update after removing from reacting_imexop
                        &input [ndx::species+alfa][0][m - dkbx][n - dkbz],
                        beta,
                        &output[ndx::species+alfa][0][m - dkbx][n - dkbz]);
                }

                // if BL, replace top
                if (grid.one_sided()) {
                    // protect against 0*NaN = Nan
                    if (std::norm(beta)==0) {
                        for (std::size_t i=0; i<input.shape()[0]; ++i) {
                            output[i][Ny-1][m - dkbx][n - dkbz] = top_input[i];
                        }
                    }
                    else {
                        for (std::size_t i=0; i<input.shape()[0]; ++i) {
                            output[i][Ny-1][m - dkbx][n - dkbz] =
                                top_input[i] + beta*top_output[i];
                        }
                    }
                }


            }
        }
        break;
    }

    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
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
    enum { nsides = 2, nmomentum = 3 };

    // Indices within PA^TP^T at which to apply boundary conditions
    // Computed once within constructor and then repeatedly used
    int rho[nsides], noslip[nsides][nmomentum], e[nsides];

    // Number of wall BCs (must be 1 or 2 depending on channel or
    // boundary layer)
    int nwalls;

    // Precomputed coefficient based on the isothermal equation of state
    real_t e_tot;

    // Store velocity values in a vector to loop through components
    std::vector<real_t> velw;

    // To flip sign of v at upper wall for channel case
    std::vector<real_t> vsgn;
 
public:

    IsothermalNoSlipPATPTEnforcer(const suzerain_bsmbsm &A_T,
                                  const real_t &u_wall,
                                  const real_t &v_wall,
                                  const real_t &w_wall,                                  
                                  const real_t &e_tot, 
                                  const int nwalls)
        : // Pass internal energy and compute total energy to use in the object
          e_tot(e_tot+0.5*(u_wall*u_wall + v_wall*v_wall + w_wall*w_wall))
        , nwalls(nwalls)
         
    {
        // Starting offset to named scalars in interleaved_state pencil
        const int e0   = static_cast<int>(ndx::e  ) * A_T.n;
        const int mx0  = static_cast<int>(ndx::mx ) * A_T.n;
        const int my0  = static_cast<int>(ndx::my ) * A_T.n;
        const int mz0  = static_cast<int>(ndx::mz ) * A_T.n;
        const int rho0 = static_cast<int>(ndx::rho) * A_T.n;

        // Relative to foo0 what is the offset to lower, upper walls
        const int wall[nsides] = { 0, A_T.n - 1};

        // Prepare indices within PA^TP^T corresponding to the walls.
        // Uses that {A^T}_{i,j} maps to {PA^TP^T}_{{q^-1}(i),{q^(-1)}(j)}.
        for (int i = 0; i < nsides; ++i) {
            e     [i]    = suzerain_bsmbsm_qinv(A_T.S, A_T.n, e0   + wall[i]);
            noslip[i][0] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, mx0  + wall[i]);
            noslip[i][1] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, my0  + wall[i]);
            noslip[i][2] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, mz0  + wall[i]);
            rho   [i]    = suzerain_bsmbsm_qinv(A_T.S, A_T.n, rho0 + wall[i]);
        }

        // Initialize wall velocity vector
        velw.resize(3);
        velw[0] = u_wall;
        velw[1] = v_wall;
        velw[2] = w_wall;
        
        // Initialize sign constant for velocity components
        // NOTE: we need to declare the sign for the three components
        // Indices go by eqn+nmomentum*wall, with 
        // eqn {u,v,w}       = {1,2,3} and
        // wall{bottom, top} = {0,1}
        vsgn.resize(nsides*nmomentum);
        // u-velocity
        vsgn[0] =  1;
        vsgn[3] =  1;
        // v-velocity
        vsgn[1] =  1;
        vsgn[4] = -1;
        // w-velocity
        vsgn[2] =  1;
        vsgn[5] =  1;
    }

    /**
     * Zero RHS of momentum and energy equations at wall BCs.
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
                // Necessary to ensure constraint possible in degenerate case
                complex_t &rhocoeff = col[rho[wall]];
                if (rhocoeff == complex_t(0)) rhocoeff = 1;
                for (int i = begin; i < end; ++i) {
                    if (i == noslip[wall][eqn]) {
                        col[i] = rhocoeff;
                    } else if (i == rho[wall]) {
                        col[i] = - rhocoeff * vsgn[eqn+nmomentum*wall] * velw[eqn];
                    } else {
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
                    col[i] = -rhocoeff/e_tot; // TODO: Check me!
                } else if (i == rho[wall]) {
                    // NOP
                } else {
                    col[i] = 0;
                }
            }
        }
    }
};


/**
 * A functor for applying mass fraction for the species PA^TP^T-based
 * operators.  Encapsulated in a class so that the same logic may be
 * applied to real- and complex-valued operators and right hand sides.
 *
 */
class MassFractionPATPTEnforcer
{
    enum { nsides = 2, nmomentum = 3 };

    // Indices within PA^TP^T at which to apply boundary conditions
    // Computed once within constructor and then repeatedly used
    int rho_s[nsides];
    int Ny;

    // Number of wall BCs (must be 1 or 2 depending on channel or
    // boundary layer)
    int nwalls;

    // Precomputed coefficient based on the isothermal equation of state
    const std::vector<real_t>& wall_mass_fractions;

public:

    MassFractionPATPTEnforcer(const suzerain_bsmbsm &A_T,
                              const std::vector<real_t>& wall_mass_fractions,
                              const int nwalls)
        : Ny(A_T.n)
        , wall_mass_fractions(wall_mass_fractions)
        , nwalls(nwalls)
    {
        // Starting offset to named scalars in interleaved_state pencil
        const int rho_s0 = 0;

        // Relative to foo0 what is the offset to lower, upper walls
        const int wall[nsides] = { 0, A_T.n - 1};

        // Prepare indices within PA^TP^T corresponding to the walls.
        // Uses that {A^T}_{i,j} maps to {PA^TP^T}_{{q^-1}(i),{q^(-1)}(j)}.
        for (int i = 0; i < nsides; ++i) {
            rho_s[i] = suzerain_bsmbsm_qinv(A_T.S, A_T.n, rho_s0 + wall[i]);
        }
    }

    /**
     * Set RHS to rho*wall_mass_fractions
     */
    template<typename T>
    void rhs(T * const bspec, T * const xflow, const int alfa)
    {
        // Note two assumptions...
        // 1.) alfa+1 implies wall_mass_fractions[0] is diluter
        // 2.) incoming xflow is non-permuted flow solution

        for (int wall = 0; wall < nwalls; ++wall) {
            bspec[rho_s[wall]] = 
                xflow[ndx::rho*Ny+wall*(Ny-1)]*wall_mass_fractions[alfa+1];
        }
    }

    /**
     * Set RHS to rho*wall_mass_fractions, non-permuted version
     */
    template<typename T>
    void rhs_non_permuted(T * const bspec, T * const xflow, const int alfa)
    {
        // Note two assumptions...
        // 1.) alfa+1 implies wall_mass_fractions[0] is diluter
        // 2.) incoming xflow is non-permuted flow solution

        for (int wall = 0; wall < nwalls; ++wall) {
            bspec[wall*(Ny-1)] = 
                xflow[ndx::rho*Ny+wall*(Ny-1)]*wall_mass_fractions[alfa+1];
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
        // No attempt made to not unnecessarily disturb matrix conditioning.

        for (int wall = 0; wall < nwalls; ++wall) {
            int begin, end;

            T * const col = (T *) suzerain_gbmatrix_col(
                    A_T.N, A_T.N, A_T.KL, A_T.KU, (void *) patpt, patpt_ld,
                    sizeof(T), rho_s[wall], &begin, &end);
            // Scan row and adjust coefficients for constraint
            for (int i = begin; i < end; ++i) {
                if (i == rho_s[wall]) {
                    col[i] = 1.0;
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
        const lowstorage::method_interface<complex_t> &method,
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

    // number of flow equations
    enum {Nf=5};

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    assert(static_cast<int>(ndx::e  ) < flow_solver->S);
    assert(static_cast<int>(ndx::mx ) < flow_solver->S);
    assert(static_cast<int>(ndx::my ) < flow_solver->S);
    assert(static_cast<int>(ndx::mz ) < flow_solver->S);
    assert(static_cast<int>(ndx::rho) < flow_solver->S);

    // TODO: Assert species solver sized correctly

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
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned)             Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==                         1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned)             Ny);
    SUZERAIN_ENSURE(state.shape()  [0] == (unsigned) (flow_solver->S +
                                                      species_solver.size()));

    // Compute how many additional mean constraints we must solve
    // Ensure conformant, mean constraints are arriving on the correct rank
    const std::size_t nconstraints = ic0 ? ic0->shape()[2]*ic0->shape()[3] : 0;
    if (nconstraints) {
        SUZERAIN_ENSURE(dgrid.has_zero_zero_modes());
        SUZERAIN_ENSURE(ic0->shape()  [1] == (unsigned)             Ny); 
        SUZERAIN_ENSURE(ic0->strides()[1] ==                         1);
        SUZERAIN_ENSURE(ic0->strides()[0] == (unsigned)             Ny);
        SUZERAIN_ENSURE(ic0->shape()  [0] == (unsigned) (flow_solver->S +
                                                         species_solver.size()));
    }

    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions

    // Pack reference details for suzerain_reacting_imexop routines
    suzerain_reacting_imexop_scenario s(this->imexop_s());
    suzerain_reacting_imexop_ref   ref;
    suzerain_reacting_imexop_refld ld;
    common.imexop_ref(ref, ld);

    // Prepare an almost functor mutating RHS and PA^TP^T to enforce BCs.
    // FIXME: permit upper_T, upper_cs, and upper_{u,v,w}
    int nwalls=2;
    if (grid.one_sided()) nwalls=1;

    IsothermalNoSlipPATPTEnforcer flow_bc_enforcer(
        *flow_solver, 
        isospec.lower_u, isospec.lower_v, isospec.lower_w, 
        cmods.e_from_T(isospec.lower_T, isospec.lower_cs), nwalls);

    shared_ptr<MassFractionPATPTEnforcer> species_bc_enforcer;
    if (species_solver.size()>0) {
        // FIXME: permit upper_T and upper_cs
        species_bc_enforcer = 
            make_shared<MassFractionPATPTEnforcer>(
                *(species_solver[0]), isospec.lower_cs, nwalls);
    }

    // Prepare a scratch buffer for packc/packf usage
    //
    // NOTE: buf sized for usage with flow operator is always big
    // enough for usage with decoupled species operators... beware if
    // you do anything else
    ArrayXXc buf(flow_solver->ld, flow_solver->n);

    // Iterate across local wavenumbers and "invert" operator "in-place"
    //
    // The loop structure changes markedly on rhome_xyz vs rhome_y,
    // but the basic sequence of steps remains identical in each case.
    switch (common.linearization) {

    case linearize::rhome_y:
    {
        // Form complex-valued, wavenumber-independent PA^TP^T
        static const char trans = 'T';
        if (flow_solver->spec.in_place()) { // Pack for in-place LU
            SUZERAIN_TIMER_SCOPED("implicit operator assembly (packf)");
            suzerain_reacting_flow_imexop_packf(
                    phi, &s, &ref, &ld, cop.get(),
                    ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                    buf.data(), flow_solver.get(), flow_solver->LU.data());
        } else {                       // Pack for out-of-place LU
            SUZERAIN_TIMER_SCOPED("implicit operator assembly (packc)");
            suzerain_reacting_flow_imexop_packc(
                    phi, &s, &ref, &ld, cop.get(),
                    ndx::e, ndx::mx, ndx::my, ndx::mz, ndx::rho,
                    buf.data(), flow_solver.get(), flow_solver->PAPT.data());
        }

        for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
            if (species_solver[alfa]->spec.in_place()) { // Pack for in-place LU
                SUZERAIN_TIMER_SCOPED("implicit operator assembly (packf)");
                suzerain_reacting_species_imexop_packf(
                    phi, &s, &ref, &ld, cop.get(),
                    buf.data(), species_solver[alfa].get(), 
                    species_solver[alfa]->LU.data());
            } else {                       // Pack for out-of-place LU
                SUZERAIN_TIMER_SCOPED("implicit operator assembly (packc)");
                suzerain_reacting_species_imexop_packc(
                    phi, &s, &ref, &ld, cop.get(),
                    buf.data(), species_solver[alfa].get(), 
                    species_solver[alfa]->PAPT.data());
            }
        }

        // Apply boundary conditions to PA^TP^T
        {
            SUZERAIN_TIMER_SCOPED("implicit operator BCs");
            flow_bc_enforcer.op(*flow_solver, flow_solver->PAPT.data(),
                                         flow_solver->PAPT.colStride());
        }

        // If one-sided, set 5x5 PAPT block corresponding to
        // freestream boundary (i.e., the top) to identity matrix
        if (grid.one_sided()) {
            const int nr  = flow_solver->KL+flow_solver->KU+1;
            const int ncb = flow_solver->S*(flow_solver->n-1);

            // zero the whole block
            flow_solver->PAPT.block(0              , ncb, nr, Nf).setZero();

            // set diagonal to 1
            flow_solver->PAPT.block(flow_solver->KL, ncb, 1 , Nf).setOnes();
        }

        // Inform the flow_solver about the new, unfactorized operator
        flow_solver->supplied_PAPT();


        for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
            
            {
                // Apply species BCs to operator
                SUZERAIN_TIMER_SCOPED("implicit operator BCs"); 
                species_bc_enforcer->op(*(species_solver[alfa]), 
                                        species_solver[alfa]->PAPT.data(),
                                        species_solver[alfa]->PAPT.colStride());
            }

            // If one-sided, set 1x1 PAPT block corresponding to
            // freestream boundary (i.e., the top) to identity matrix
            if (grid.one_sided()) {
                const int KL = species_solver[alfa]->KL;
                const int KU = species_solver[alfa]->KU;

                const int nr  = KL+KU+1;
                const int ncb = 
                    species_solver[alfa]->S*(species_solver[alfa]->n-1);

                // zero the whole block
                species_solver[alfa]->PAPT.block(0 , ncb, nr, 1).setZero();
                // set diagonal to 1
                species_solver[alfa]->PAPT.block(KL, ncb, 1 , 1).setOnes();
            }

            // Inform species_solver about new, unfactorized operator
            species_solver[alfa]->supplied_PAPT();
        }

        // Species solver supplied call

        // Solve using a new right hand side for each wavenumber pair km, kn
        for (int n = dkbz; n < dkez; ++n) {
            const int wn = wavenumber(dNz, n);

            for (int m = dkbx; m < dkex; ++m) {
                const int wm = wavenumber(dNx, m);

                // Get pointer to (.,m,n)-th state pencil
                complex_t * const p = &state[0][0][m - dkbx][n - dkbz];

                // Continue didn't yet occur for Nyquist/dealiasing modes...
                if (   std::abs(wn) > wavenumber_absmin(Nz)
                    || std::abs(wm) > wavenumber_absmin(Nx)) {
                    memset(p, 0, flow_solver->N*sizeof(p[0]));  // ...so we can zero,
                    continue;                              // and short circuit.
                }

                // Form right hand side, apply BCs, factorize, and solve.
                // Beware that much ugly, ugly magic is hidden just below.

                // Flow RHS, BCs, solve
                flow_solver->supply_B(p);

                // Boundary conditions
                {
                    SUZERAIN_TIMER_SCOPED("implicit right hand side BCs");
                    flow_bc_enforcer.rhs(flow_solver->PB.data());
                }

                flow_solver->solve(trans);
                flow_solver->demand_X(p);


                // Species RHSs, BCs, solves
                for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
                    species_solver[alfa]->supply_B(p + (ndx::species+alfa)*Ny);

                    {
                        SUZERAIN_TIMER_SCOPED("implicit right hand side BCs");
                        species_bc_enforcer->rhs(species_solver[alfa]->PB.data(),
                                                 p, alfa);
                    }

                    species_solver[alfa]->solve(trans);
                    species_solver[alfa]->demand_X(p + (ndx::species+alfa)*Ny);
                }

            }
        }

        // Compute total state size for pencil
        int Ntot = flow_solver->N;
        if (species_solver.size()>0) {
            // All species have same size
            Ntot += species_solver.size()*species_solver[0]->N;
        }

        // If necessary, solve any required integral constraints
        for (std::size_t i = 0; i < nconstraints; ++i) {
            SUZERAIN_TIMER_SCOPED("implicit constraint solution");
            flow_solver->supply_B(ic0->data() + i * Ntot);
            flow_bc_enforcer.rhs(flow_solver->PB.data());
            flow_solver->solve(trans);
            flow_solver->demand_X(ic0->data() + i * Ntot);


            // Handle effects of constraints on species wall BCs
            for (std::size_t alfa=0; alfa<species_solver.size(); ++alfa) {
                species_bc_enforcer->rhs_non_permuted(
                    ic0->data() + i * Ntot + (ndx::species+alfa)*Ny, 
                    ic0->data() + i * Ntot, 
                    alfa);
            }
        }

        break;
    }


    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

    // State leaves method as coefficients in X, Y, and Z directions
}

suzerain_reacting_imexop_scenario
isothermal_hybrid_linear_operator::imexop_s() const
{
    suzerain_reacting_imexop_scenario retval;
    retval.alpha = cmods.alpha;
    return retval;
}

} // namespace reacting

} // namespace suzerain
