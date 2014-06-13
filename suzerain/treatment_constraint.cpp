//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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
 * @copydoc treatment_constraint.hpp
 */

#include <suzerain/treatment_constraint.hpp>

#include <suzerain/bspline.h>
#include <suzerain/constraint.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>

namespace suzerain {

namespace constraint {

treatment::inputs::~inputs()
{
}

treatment::outputs::~outputs()
{
}

treatment::treatment(const real_t& Ma,
                     const pencil_grid& dgrid,
                     bspline& b,
                     inputs& inp,
                     outputs& out)
    : none(new constraint::disabled(b))
    , Ma(Ma)
    , rank_has_zero_zero_modes(dgrid.has_zero_zero_modes())
    , inp(inp)
    , out(out)
    , jacobiSvd(0, 0, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    using std::fill;
    fill(physical.begin(),  physical.end(),  none);
    fill(numerical.begin(), numerical.end(), none);
}

void
treatment::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    return L->apply_mass_plus_scaled_operator(
            phi, state, substep_index);
}

void
treatment::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    return L->accumulate_mass_plus_scaled_operator(
            phi, input, beta, output, substep_index);
}

void
treatment::invert_mass_plus_scaled_operator(
        const complex_t& phi,
        multi_array::ref<complex_t,4>& state,
        const lowstorage::method_interface<complex_t>& method,
        const real_t delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("constraint::treatment");

    // Shorthand
    using ndx::e;
    using ndx::mx;
    using ndx::my;
    using ndx::mz;
    using ndx::rho;
    const int rho_1 = ndx::rho_(1);

    // Sidesteps assertions when local rank contains no wavespace information
    const int N = (int) state.shape()[1];
    if (SUZERAIN_UNLIKELY(0 == N)) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    // Any amount of incoming state is valid so long as there's enough there
    SUZERAIN_ENSURE(state.strides()[1] ==            1  );
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) N  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) e  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) rho);

    // See treatment_channel writeup (redux) for information on steps below.

    // Have a tantrum if caller expects us to compute any constraints
    SUZERAIN_ENSURE(ic0 == NULL);

    // On zero-zero rank, re-use ic0 to wrap cdata for BaseClass invocation
    if (rank_has_zero_zero_modes) {

        // Prepare RHS data for density, momentum, and total energy constraints
        //
        // Notice scaling by Mach^2 to cause mx-, my-, and mz-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        //
        // Notice ndx::rho constraint ignores other equations BY DESIGN.
        // This choice is reflected in the later update of out.Crho.
        cdata.setZero(state.shape()[0]*N, cdata.cols());
        cdata.col(e  ).segment(N*e  , N).real() = physical[e  ]->shape;
        cdata.col(mx ).segment(N*e  , N).real() = physical[mx ]->shape
                                                * inp.u() * (Ma * Ma);
        cdata.col(mx ).segment(N*mx , N).real() = physical[mx ]->shape;
        cdata.col(my ).segment(N*e  , N).real() = physical[my ]->shape
                                                * inp.v() * (Ma * Ma);
        cdata.col(my ).segment(N*my , N).real() = physical[my ]->shape;
        cdata.col(mz ).segment(N*e  , N).real() = physical[mz ]->shape
                                                * inp.w() * (Ma * Ma);
        cdata.col(mz ).segment(N*mz , N).real() = physical[mz ]->shape;
        cdata.col(rho).segment(N*rho, N).real() = physical[rho]->shape;

        // Prepare RHS for one decoupled numerical constraint per equation
        cdata.col(rho_1+e  ).segment(N*e  , N).real() = numerical[e  ]->shape;
        cdata.col(rho_1+mx ).segment(N*mx , N).real() = numerical[mx ]->shape;
        cdata.col(rho_1+my ).segment(N*my , N).real() = numerical[my ]->shape;
        cdata.col(rho_1+mz ).segment(N*mz , N).real() = numerical[mz ]->shape;
        cdata.col(rho_1+rho).segment(N*rho, N).real() = numerical[rho]->shape;

        // Wrap data into appropriately digestible format
        using std::size_t;
        const array<size_t,4> sizes = {{
                state.shape()[0], (size_t) N, (size_t) cdata.cols(), 1
        }};
        ic0 = new multi_array::ref<complex_t,4>(
                cdata.data(), sizes, storage::interleaved<4>());
    }

    // Delegate to wrapped operator for the solution procedure
    L->invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // Clean up any integral constraint data wrapper we may have employed
    delete ic0;

    // Only the rank with zero-zero modes proceeds!
    if (!rank_has_zero_zero_modes) return;

    // Get an Eigen-friendly map of the zero-zero mode coefficients
    Map<VectorXc> mean(state.origin(), state.shape()[0]*N);

    // Solve the requested, possibly simultaneous constraint problem.  A fancy
    // decomposition is used for a simple 10x10 solve or to permit one or more
    // inactive constraints via least squares.  Least squares also adds
    // robustness whenever the constraints are somehow incompatible.
    //
    // 1) Assemble the matrix problem for simultaneous integral constraints.
    //    Nonzero only when both the i-th and j-th constraints are enabled.
    MatrixAr cmat(MatrixAr::Zero());
    for (int j = 0; j < cdata.cols(); ++j) {
        if (j<rho_1 ? physical[j]->enabled() : numerical[j-rho_1]->enabled()) {
            for (int i = 0; i < physical.static_size;  ++i) {
                if (physical[i]->enabled()) {
                    cmat(i,j) = physical[i]->coeff.dot(
                            cdata.col(j).segment(N*i,N).real());
                }
            }
            for (int i = 0; i < numerical.static_size; ++i) {
                if (numerical[i]->enabled()) {
                    cmat(rho_1+i,j) = numerical[i]->coeff.dot(
                            cdata.col(j).segment(N*i,N).real());
                }
            }
        }
    }
    // 2) Prepare initial values relative to desired targets when active.
    //    Otherwise, zero the associated row and column in matrix.
    VectorAr crhs(VectorAr::Zero());
    for (int i = 0; i < physical.static_size;  ++i) {
        if (physical[i]->enabled()) {
            crhs[i] = physical[i]->target()
                    - physical[i]->coeff.dot(
                            mean.segment(N*i,N).real());
        }
    }
    for (int i = 0; i < numerical.static_size;  ++i) {
        if (numerical[i]->enabled()) {
            crhs[rho_1+i] = numerical[i]->target()
                          - numerical[i]->coeff.dot(
                                  mean.segment(N*i,N).real());
        }
    }
    // 3) Prepare matrix decomposition and solve using least squares
    VectorAr cphi = jacobiSvd.compute(cmat).solve(crhs);
    // 4) Add correctly scaled results to the mean state to satisfy constraints
    mean += (cdata * cphi.asDiagonal()).rowwise().sum();

    // The implicitly applied integral constraints, as point values, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations.
    //
    // Per discussion with M. K. Lee, the implicit operator used to compute cphi
    // already accounts for variable substep effects and so only a simple
    // in-place running average is necessary:
    //
    //    mean += (nsubsteps*(sample / delta_t) - mean) / (substep_index + 1)
    //
    // The delta_t accounts for step sizes already implicitly included in cphi.
    // The correctness of this approach can be confirmed by running a 1D channel
    // known a prior to be stationary followed by checking the total stress.
    //
    // Notice mx-related forcing is NOT scaled by Mach^2 when tracked
    // because our post-processing routines will account for Mach^2 factor.
    //
    // Notice physical[ndx::rho] constraint lumped into Crho BY DESIGN.
    cphi               *= method.substeps() / delta_t;
    const real_t iota   = static_cast<real_t>(1) / (substep_index + 1);
    out.fx()           += iota*(
                              cphi[mx]*physical[mx]->shape
                            - out.fx()
                          );
    out.fy()           += iota*(
                              cphi[my]*physical[my]->shape
                            - out.fy()
                          );
    out.fz()           += iota*(
                              cphi[mz]*physical[mz]->shape
                            - out.fz()
                          );
    out.f_dot_u()      += iota*(
                              cphi[mx]*physical[mx]->shape*inp.u()
                            + cphi[my]*physical[my]->shape*inp.v()
                            + cphi[mz]*physical[mz]->shape*inp.w()
                            - out.f_dot_u()
                          );
    out.qb()           += iota*(
                              cphi[e]*physical[e]->shape
                            - out.qb()
                          );
    out.CrhoE()        += iota*(
                              cphi[rho_1+e  ]*numerical[e  ]->shape
                            - out.CrhoE()
                          );
    out.C2rhoE()       += iota*(
                              (cphi[rho_1+e  ]*numerical[e  ]->shape).square()
                            - out.C2rhoE()
                          );
    out.Crhou()        += iota*(
                              cphi[rho_1+mx ]*numerical[mx ]->shape
                            - out.Crhou()
                          );
    out.C2rhou()       += iota*(
                              (cphi[rho_1+mx ]*numerical[mx ]->shape).square()
                            - out.C2rhou()
                          );
    out.Crhov()        += iota*(
                              cphi[rho_1+my ]*numerical[my ]->shape
                            - out.Crhov()
                          );
    out.C2rhov()       += iota*(
                              (cphi[rho_1+my ]*numerical[my ]->shape).square()
                            - out.C2rhov()
                          );
    out.Crhow()        += iota*(
                              cphi[rho_1+mz ]*numerical[mz ]->shape
                            - out.Crhow()
                          );
    out.C2rhow()       += iota*(
                              (cphi[rho_1+mz ]*numerical[mz ]->shape).square()
                            - out.C2rhow()
                          );
    out.Crho()         += iota*(
                              cphi[rho]      *physical [rho]->shape
                            + cphi[rho_1+rho]*numerical[rho]->shape
                            - out.Crho()
                          );
    out.C2rho()        += iota*( // Terms like aa + 2ab + bb
                                (cphi[rho]      *physical [rho]->shape).square()
                            + 2*(cphi[rho]      *physical [rho]->shape)
                               *(cphi[rho_1+rho]*numerical[rho]->shape)
                            +   (cphi[rho_1+rho]*numerical[rho]->shape).square()
                            - out.C2rho()
                          );
    out.Crhou_dot_u()  += iota*(
                              cphi[rho_1+mx]*numerical[mx]->shape*inp.u()
                            + cphi[rho_1+my]*numerical[my]->shape*inp.v()
                            + cphi[rho_1+mz]*numerical[mz]->shape*inp.w()
                            - out.Crhou_dot_u()
                          );
    out.C2rhou_dot_u() += iota*( // Terms like uu + 2uv + 2uw + vv + 2vw + ww
                                (cphi[rho_1+mx]*numerical[mx]->shape).square()
                               *inp.uu()
                            + 2*(cphi[rho_1+mx]*cphi[rho_1+my])
                               *numerical[mx]->shape
                               *numerical[my]->shape
                               *inp.uv()
                            + 2*(cphi[rho_1+mx]*cphi[rho_1+mz])
                               *numerical[mx]->shape
                               *numerical[mz]->shape
                               *inp.uw()
                            +   (cphi[rho_1+my]*numerical[my]->shape).square()
                               *inp.vv()
                            + 2*(cphi[rho_1+my]*cphi[rho_1+mz])
                               *numerical[my]->shape
                               *numerical[mz]->shape
                               *inp.vw()
                            +   (cphi[rho_1+mz]*numerical[mz]->shape).square()
                               *inp.ww()
                            - out.C2rhou_dot_u()
                          );

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace constraint

} // namespace suzerain
