//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc layers.hpp
 */

#include "layers.hpp"

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bl.h>
#include <suzerain/error.h>
#include <suzerain/grid_specification.hpp>
#include <suzerain/largo_state.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/largo_definition.hpp>

#include "quantities.hpp"
#include "scenario_definition.hpp"

using boost::numeric_cast;

namespace suzerain {

namespace perfect {

layers::layers()
{
    // NOP
}

layers::layers(
        layers::storage_type::Index Ny)
    : storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
{
    // NOP
}

layers& layers::operator=(const quantities &q)
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
    this->T()            = q.T();
    this->mu()           = q.mu();
    this->u().col(0)     = q.u().col(0);
    this->u().col(1)     = q.u().col(1);

    return *this;
}

// This logic is a trimmed down version of quantities.cpp. See comments there.
layers sample_layers(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        contiguous_state<4,complex_t> &swave)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("sample_layers");

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

    // Rank-specific details accumulated in ret to be MPI_Reduce-d later
    layers ret(Ny);

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret.rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret.rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret.rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Transform to obtain physical space view of state on collocation points
    physical_view<state_count>sphys(dgrid, swave);
    operator_tools otool(grid, dgrid, cop);
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
                SUZERAIN_PERFECT_LAYERS_PHYSICAL)
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
                SUZERAIN_PERFECT_LAYERS_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // layers while other ranks have zeros in those locations.

    // Reduce sums onto rank zero and then return garbage from non-zero ranks
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {

        // Reduce operation requires no additional storage on rank-zero
        SUZERAIN_MPICHKR(MPI_Reduce(
                MPI_IN_PLACE, ret.storage.data(), ret.storage.size(),
                mpi::datatype<layers::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

    } else {

        // Reduce operation requires temporary storage on non-zero ranks
        ArrayXXr tmp;
        tmp.resizeLike(ret.storage);
        tmp.setZero();
        SUZERAIN_MPICHKR(MPI_Reduce(ret.storage.data(), tmp.data(), tmp.size(),
                mpi::datatype<layers::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

        // Force non-zero ranks contain all NaNs to help detect usage errors
        ret.storage.fill(std::numeric_limits<
                layers::storage_type::Scalar>::quiet_NaN());

        // Return from all non-zero ranks
        return ret;

    }

    // Only rank zero reaches this logic because of return statement just above
    assert(mpi::comm_rank(MPI_COMM_WORLD) == 0);

    // Physical space sums, which are at collocation points, need to be
    // divided by the dealiased extents and converted to coefficients.
    const real_t scale_factor = 1 / dgrid.chi();
    bsplineop_lu scaled_mass(cop);
    scaled_mass.opform(1, &scale_factor, cop);
    scaled_mass.factor();
    scaled_mass.solve(layers::nscalars::physical,
            ret.storage.middleCols<layers::nscalars::physical>(
                layers::start::physical).data(),
            ret.storage.innerStride(), ret.storage.outerStride());

    return ret;
}

void summarize_boundary_layer_nature(
        const layers &lay,
        const scenario_definition &scenario,
        const shared_ptr<support::largo_definition> &sg,
        bspline &b,
        suzerain_bl_local       &wall,
        suzerain_bl_viscous     &viscous,
        suzerain_bl_local       &edge,
        suzerain_bl_thicknesses &thick,
        suzerain_bl_qoi         &qoi,
        suzerain_bl_pg          &pg)
{
    // Care taken to avoid any B-spline evaluation using NaN thicknesses
    // per Redmine #2940 as this can bring down a simulation from GSL_ERROR.
    using boost::math::isnan;

    // Compute boundary layer thicknesses, including delta
    suzerain_bl_compute_thicknesses(lay.H0().data(),
                                    lay.rho_u().col(0).data(),
                                    lay.u().col(0).data(),
                                    &thick, b.bw, b.dbw);

    // Prepare local state at the wall (y=0)
    // Uses B-spline 0th coefficient being wall value to reduce costs.
    std::fill_n(reinterpret_cast<double *>(&wall),
                sizeof(wall)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    wall.a     = lay.a()[0];
    wall.gamma = scenario.gamma;
    wall.mu    = lay.mu()[0];
    wall.Pr    = scenario.Pr;
    wall.rho   = lay.rho()[0];
    wall.T     = lay.T()[0];
    assert(&(wall.u) + 1 == &(wall.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, lay.u().col(0).data(), 0.0, &(wall.u), 1);
    wall.v     = lay.u().col(1)[0];

    // Compute viscous quantities based only on wall information
    suzerain_bl_compute_viscous(&wall, &viscous);

    // Evaluate state at the edge (y=thick.delta) from B-spline coefficients
    std::fill_n(reinterpret_cast<double *>(&edge),
                sizeof(edge)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    edge.gamma = scenario.gamma;
    edge.Pr    = scenario.Pr;
    if (SUZERAIN_UNLIKELY((isnan)(thick.delta))) {
        // NOP as fill_n above ensured NaN propagated correctly
    } else {
        // Later, could more quickly evaluate basis once and then re-use that
        // result to repeatedly form the necessary linear combinations.
        const double delta = thick.delta;
        b.linear_combination(0, lay.a().data(),        delta, &(edge.a));
        b.linear_combination(0, lay.mu().data(),       delta, &(edge.mu));
        b.linear_combination(0, lay.rho().data(),      delta, &(edge.rho));
        b.linear_combination(0, lay.T().data(),        delta, &(edge.T));
        assert(&(edge.u) + 1 == &(edge.u__y)); // Next, compute both u and u__y
        b.linear_combination(1, lay.u().col(0).data(), delta, &(edge.u), 1);
        b.linear_combination(0, lay.u().col(1).data(), delta, &(edge.v));
    }

    // Compute general quantities of interest
    suzerain_bl_compute_qoi(scenario.Ma, scenario.Re,
                            &wall, &viscous, &edge, & thick, &qoi);

    // Mean pressure and streamwise velocity gradients come from slow growth
    double edge_p__x = 0, edge_u__x = 0;
    if (sg && sg->formulation.enabled()) {
        const double delta = thick.delta;
        if (SUZERAIN_UNLIKELY((isnan)(delta))) {
            edge_p__x = edge_u__x = std::numeric_limits<real_t>::quiet_NaN();
        } else {
            largo_state base, dy, dx; // as_is()
            sg->get_baseflow(delta, base.as_is(), dy.as_is(), dx.as_is());
            sg->get_baseflow_pressure(delta, base.p, dy.p, dx.p);
            edge_p__x = dx.p;                                       // Direct
            edge_u__x = (dx.mx - dx.mx/base.rho*dx.rho) / base.rho; // Chained
        }
    }
    suzerain_bl_compute_pg(scenario.Ma, scenario.Re, &wall, &viscous, &edge,
                           edge_p__x, edge_u__x, &thick, &pg);
}

} // namespace perfect

} // namespace suzerain
