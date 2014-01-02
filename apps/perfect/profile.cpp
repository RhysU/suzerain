//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc profile.hpp
 */

#include "profile.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bl.h>
#include <suzerain/bspline.hpp>
#include <suzerain/channel.h>
#include <suzerain/error.h>
#include <suzerain/grid_specification.hpp>
#include <suzerain/largo_specification.hpp>
#include <suzerain/largo_state.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>

#include "quantities.hpp"
#include "scenario_definition.hpp"

using boost::numeric_cast;

namespace suzerain {

namespace perfect {

profile::profile()
{
}

profile::profile(
        profile::storage_type::Index Ny)
    : storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
{
}

profile& profile::operator=(const quantities &q)
{
    // Resize our storage and defensively NaN in case we miss something
    this->storage.setConstant(q.storage.rows(),
            storage_type::ColsAtCompileTime,
            std::numeric_limits<storage_type::Scalar>::quiet_NaN());

    // Copy boundary layer profiles of interest from q
    this->rho()          = q.rho();
    this->rho_u().col(0) = q.rho_u().col(0);
    this->rho_u().col(1) = q.rho_u().col(1);
    this->a()            = q.a();
    this->H0()           = q.H0();
    this->ke()           = q.ke();
    this->T()            = q.T();
    this->mu()           = q.mu();
    this->u().col(0)     = q.u().col(0);
    this->u().col(1)     = q.u().col(1);

    return *this;
}

// This logic is a trimmed down version of quantities.cpp. See comments there.
profile sample_profile(
        const scenario_definition &scenario,
        const operator_tools& otool,
        contiguous_state<4,complex_t> &swave)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("sample_profile");

    // Shorthand for the operator_tools member(s) commonly used
    const pencil_grid &dgrid = otool.dgrid;

    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };
    assert(static_cast<int>(ndx::e  ) < state_count);
    assert(static_cast<int>(ndx::mx ) < state_count);
    assert(static_cast<int>(ndx::my ) < state_count);
    assert(static_cast<int>(ndx::mz ) < state_count);
    assert(static_cast<int>(ndx::rho) < state_count);

    // Sanity check incoming swave's shape and contiguity
    const std::size_t Ny = swave.shape()[1];
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);

    // Rank-specific details accumulated in ret to be MPI_Allreduce-d later
    profile ret(Ny);

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret.rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret.rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret.rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Transform to obtain physical space view of state on collocation points
    physical_view<state_count>sphys(dgrid, swave);
    for (std::size_t i = 0; i < state_count; ++i) {
        otool.zero_dealiasing_modes(swave, i);
        otool.bop_apply(0, 1., swave, i);
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The two loop structure is used as positions x(i), and z(k) unnecessary.
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Type used to accumulate means versus wall-normal position
        typedef boost::accumulators::accumulator_set<
                real_t, boost::accumulators::stats<
                    boost::accumulators::tag::sum_kahan
                >
            > accumulator_type;

        // Accumulators for each mean quantity computed in physical space.
        // For example, quantity "foo" has accumulator "sum_foo".
        // Declared within Y loop so they are reset on each Y iteration.
#define DECLARE(r, data, tuple)                                              \
        accumulator_type BOOST_PP_CAT(sum_,BOOST_PP_TUPLE_ELEM(2, 0, tuple)) \
                [BOOST_PP_TUPLE_ELEM(2, 1, tuple)];
        BOOST_PP_SEQ_FOR_EACH(DECLARE,,
                SUZERAIN_PERFECT_PROFILE_PHYSICAL)
#undef DECLARE

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack conserved state
            const real_t   e(sphys(ndx::e,   offset));
            const Vector3r m(sphys(ndx::mx,  offset),
                             sphys(ndx::my,  offset),
                             sphys(ndx::mz,  offset));
            const real_t rho(sphys(ndx::rho, offset));

            // Compute quantities of interest
            const Vector3r u = rholut::u(rho, m);
            real_t p, T, mu, lambda;
            rholut::p_T_mu_lambda(scenario.alpha, scenario.beta,
                                  scenario.gamma, scenario.Ma,
                                  rho, m, e, p, T, mu, lambda);

            // Accumulate into sum_XXX using function syntax.
            sum_a [0](std::sqrt(T));
            sum_H0[0]((e + p) / rho);
            sum_ke[0](u.squaredNorm() / 2);
            sum_mu[0](mu);
            sum_T [0](T);
            sum_u [0](u.x());
            sum_u [1](u.y());

        } // end X // end Z

        // Move y-specific sums into MPI-reduction-ready storage for y(j) using
        // Eigen comma initialization syntax.  Yes, this is not stride 1.
#define EXTRACT_SUM(z, n, q) boost::accumulators::sum(sum_##q[n])
#define MOVE_SUM_INTO_TMP(r, data, tuple)                                 \
        ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().row(j) <<                  \
            BOOST_PP_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, tuple),               \
                          EXTRACT_SUM, BOOST_PP_TUPLE_ELEM(2, 0, tuple));
        BOOST_PP_SEQ_FOR_EACH(MOVE_SUM_INTO_TMP,,
                SUZERAIN_PERFECT_PROFILE_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // profile while other ranks have zeros in those locations.

    // Allreduce to obtain global sums on every rank
    SUZERAIN_MPICHKR(MPI_Allreduce(
            MPI_IN_PLACE, ret.storage.data(), ret.storage.size(),
            mpi::datatype<profile::storage_type::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));

    // Physical space sums, which are at collocation points, need to be
    // scaled by the dealiased extents and converted to coefficients.
    ret.physical() *= dgrid.chi();
    otool.masslu()->solve(ret.physical().outerSize(),
                          ret.physical().data(),
                          ret.physical().innerStride(),
                          ret.physical().outerStride());

    return ret;
}

void summarize_boundary_layer_nature(
        const profile &prof,
        const scenario_definition &scenario,
        const shared_ptr<largo_specification> &sg,
        const bsplineop_lu &masslu,
        bspline &b,
        suzerain_bl_local       &wall,
        suzerain_bl_viscous     &viscous,
        suzerain_bl_thicknesses &thick,
        suzerain_bl_local       &edge,
        suzerain_bl_reynolds    &reynolds,
        suzerain_bl_qoi         &qoi,
        suzerain_bl_pg          &pg)
{
    // Care taken to avoid any B-spline evaluation using NaN thicknesses
    // per Redmine #2940 as this can bring down a simulation from GSL_ERROR.
    using boost::math::isnan;

    // Prepare local state at the wall (y=0)
    // Uses B-spline 0th coefficient being wall value to reduce costs.
    std::fill_n(reinterpret_cast<double *>(&wall),
                sizeof(wall)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    wall.a     = prof.a()[0];
    wall.gamma = scenario.gamma;
    wall.mu    = prof.mu()[0];
    wall.Pr    = scenario.Pr;
    wall.rho   = prof.rho()[0];
    assert(&(wall.T) + 1 == &(wall.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), 0.0, &(wall.T), 1);
    assert(&(wall.u) + 1 == &(wall.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), 0.0, &(wall.u), 1);
    wall.v     = prof.u().col(1)[0];
    wall.y     = 0.0;

    // Compute viscous quantities based only on wall information
    suzerain_bl_compute_viscous(scenario.Re, &wall, &viscous);

    // Prepare inviscid baseflow coefficients if necessary
    ArrayX4r coeffs_inviscid; // Columns are H0, rhou, u, v
    if (sg && sg->formulation.enabled() && sg->baseflow) {
        coeffs_inviscid.resize(b.n(), NoChange);

        // Retrieve values on collocation points
        largo_state state, dontcare;
        baseflow_interface &baseflow = *sg->baseflow;
        for (int j = 0; j < b.n(); ++j) {
            const real_t y_j = b.collocation_point(j);
            baseflow.conserved(y_j, state.as_is(),
                               dontcare.as_is(), dontcare.as_is());
            baseflow.pressure (y_j, state.p,
                               dontcare.p, dontcare.p);
            coeffs_inviscid(j, 0) = state.H0();
            coeffs_inviscid(j, 1) = state.mx;
            coeffs_inviscid(j, 2) = state.u();
            coeffs_inviscid(j, 3) = state.v();
        }

        // Convert to coefficients using mass matrix
        masslu.solve(coeffs_inviscid.outerSize(),
                     coeffs_inviscid.data(),
                     coeffs_inviscid.innerStride(),
                     coeffs_inviscid.outerStride());
    }

    // Compute boundary layer thicknesses, including delta
    // TODO WARN on non-SUCCESS return
    if (0 == coeffs_inviscid.size()) {
        suzerain_bl_compute_thicknesses(scenario.Ma,
                                        prof.H0().data(),
                                        prof.ke().data(),
                                        prof.rho_u().col(0).data(),
                                        prof.u().col(0).data(),
                                        &thick, b.bw, b.dbw);
    } else {
        suzerain_bl_compute_thicknesses_baseflow(scenario.Ma,
                                                 prof.H0().data(),
                                                 prof.ke().data(),
                                                 prof.rho_u().col(0).data(),
                                                 prof.u().col(0).data(),
                                                 coeffs_inviscid.innerStride(),
                                                 coeffs_inviscid.col(0).data(),
                                                 coeffs_inviscid.col(1).data(),
                                                 coeffs_inviscid.col(2).data(),
                                                 coeffs_inviscid.col(3).data(),
                                                 &thick, b.bw, b.dbw);
    }

    // Evaluate state at the edge (y=thick.delta) from B-spline coefficients
    std::fill_n(reinterpret_cast<double *>(&edge),
                sizeof(edge)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    edge.gamma = scenario.gamma;
    edge.Pr    = scenario.Pr;
    edge.y     = thick.delta;
    if (SUZERAIN_UNLIKELY((isnan)(edge.y))) {
        // NOP as fill_n above ensured NaN propagated correctly
    } else {
        // Later, could more quickly evaluate basis once and then re-use that
        // result to repeatedly form the necessary linear combinations.
        b.linear_combination(0, prof.a().data(),        edge.y, &(edge.a));
        b.linear_combination(0, prof.mu().data(),       edge.y, &(edge.mu));
        b.linear_combination(0, prof.rho().data(),      edge.y, &(edge.rho));
        assert(&(edge.T) + 1 == &(edge.T__y)); // Next, compute both T and T__y
        b.linear_combination(1, prof.T().col(0).data(), edge.y, &(edge.T), 1);
        assert(&(edge.u) + 1 == &(edge.u__y)); // Next, compute both u and u__y
        b.linear_combination(1, prof.u().col(0).data(), edge.y, &(edge.u), 1);
        b.linear_combination(0, prof.u().col(1).data(), edge.y, &(edge.v));
    }

    // Compute Reynolds numbers
    // TODO WARN on non-SUCCESS return
    if (0 == coeffs_inviscid.size()) {
        suzerain_bl_compute_reynolds(scenario.Re, &edge, &thick, &reynolds);
    } else {
        suzerain_bl_compute_reynolds_baseflow(scenario.Ma,
                                              scenario.Re,
                                              prof.H0().data(),
                                              prof.ke().data(),
                                              prof.rho_u().col(0).data(),
                                              prof.u().col(0).data(),
                                              1,
                                              coeffs_inviscid.col(0).data(),
                                              coeffs_inviscid.col(1).data(),
                                              coeffs_inviscid.col(2).data(),
                                              coeffs_inviscid.col(3).data(),
                                              &edge, &reynolds, b.bw);
    }

    // Compute general quantities of interest
    suzerain_bl_compute_qoi(scenario.Ma, scenario.Re,
                            &wall, &viscous, &edge, &thick, &qoi);

    // Mean pressure and streamwise velocity gradients come from slow growth
    double edge_p__x = 0, edge_u__x = 0;
    if (sg && sg->formulation.enabled() && sg->baseflow) {
        const double delta = thick.delta;
        if (SUZERAIN_UNLIKELY((isnan)(delta))) {
            edge_p__x = edge_u__x = std::numeric_limits<real_t>::quiet_NaN();
        } else {
            largo_state base, dy, dx; // as_is()
            sg->baseflow->conserved(delta, base.as_is(), dy.as_is(), dx.as_is());
            sg->baseflow->pressure(delta, base.p, dy.p, dx.p);
            edge_p__x = dx.p;                                       // Direct
            edge_u__x = (dx.mx - dx.mx/base.rho*dx.rho) / base.rho; // Chained
        }
    }
    suzerain_bl_compute_pg(scenario.Ma, scenario.Re, &wall, &viscous, &edge,
                           edge_p__x, edge_u__x, &thick, &pg);
}

void summarize_channel_nature(
        const profile &prof,
        const scenario_definition &scenario,
        bspline &b,
        suzerain_channel_local   &wall,
        suzerain_channel_viscous &viscous,
        suzerain_channel_local   &center,
        suzerain_channel_qoi     &qoi)
{
    // Prepare local state at the wall (y=0)
    // Uses B-spline 0th coefficient being wall value to reduce costs.
    std::fill_n(reinterpret_cast<double *>(&wall),
                sizeof(wall)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    wall.y     = b.collocation_point(0);
    wall.a     = prof.a()[0];
    wall.gamma = scenario.gamma;
    wall.mu    = prof.mu()[0];
    wall.Pr    = scenario.Pr;
    wall.rho   = prof.rho()[0];
    assert(&(wall.T) + 1 == &(wall.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), wall.y, &(wall.T), 1);
    assert(&(wall.u) + 1 == &(wall.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), wall.y, &(wall.u), 1);
    wall.v     = prof.u().col(1)[0];

    // Compute viscous quantities based only on wall information
    suzerain_channel_compute_viscous(scenario.Re, &wall, &viscous);

    // Evaluate state at the centerline from B-spline basis
    std::fill_n(reinterpret_cast<double *>(&center),
                sizeof(center)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    center.y = b.collocation_point(b.n()-1) / 2;
    b.linear_combination(0, prof.a().data(),        center.y, &(center.a));
    center.gamma = scenario.gamma;
    b.linear_combination(0, prof.mu().data(),       center.y, &(center.mu));
    center.Pr    = scenario.Pr;
    b.linear_combination(0, prof.rho().data(),      center.y, &(center.rho));
    assert(&(center.T) + 1 == &(center.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), center.y, &(center.T), 1);
    assert(&(center.u) + 1 == &(center.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), center.y, &(center.u), 1);
    b.linear_combination(0, prof.u().col(1).data(), center.y, &(center.v));

    // Compute general quantities of interest
    suzerain_channel_compute_qoi(scenario.Ma, scenario.Re,
                                 &wall, &viscous, &center, &qoi);
}

} // namespace perfect

} // namespace suzerain
