//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_HPP
#define SUZERAIN_PERFECT_HPP

/** @file
 * Support logic for the perfect gas application.
 */

#include <suzerain/common.hpp>
#include <suzerain/state_fwd.hpp>

// Forward declarations
struct suzerain_bl_local;
struct suzerain_bl_pg;
struct suzerain_bl_qoi;
struct suzerain_bl_reynolds;
struct suzerain_bl_thicknesses;
struct suzerain_bl_viscous;
struct suzerain_channel_local;
struct suzerain_channel_qoi;
struct suzerain_channel_viscous;

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class bsplineop_lu;
class operator_tools;
class pencil_grid;
class profile;
class samples;
class specification_grid;
class specification_largo;
class specification_noise;
namespace support { class field; }
template <int> struct physical_view;

/**
 * Functionality for Suzerain perfect gas applications.
 */
namespace perfect {

// Forward declarations
class definition_scenario;
class instantaneous;
class references;

/**
 * Type of Boost.Accumulator to use for summation processes.
 *
 * Kahan summation preferred when available as incremental cost is small and we
 * may add many small numbers to a large magnitude sum. During debugging, the
 * number of samples collected is also available.
 */
typedef boost::accumulators::accumulator_set<
            real_t,
            boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
#ifndef NDEBUG
                , boost::accumulators::tag::count
#endif
            >
        > summing_accumulator_type;

/** Return default nondimensional field information per \ref suzerain::ndx. */
std::vector<support::field>
default_fields();

/**
 * Hold temperature and density constant while changing the Mach number and
 * ratio of specific heats.  On input, \c state should contain total energy
 * fields using \c old_Ma and \c old_gamma.  On output \c state will contain
 * total energy fields using <tt>scenario.Ma</tt> and <tt>scenario.gamma</tt>.
 */
void
adjust_scenario(contiguous_state<4,complex_t> &swave,
                const definition_scenario& scenario,
                const specification_grid& grid,
                const pencil_grid& dgrid,
                const bsplineop& cop,
                const real_t old_Ma,
                const real_t old_gamma);

/**
 * Add random momentum field perturbations ("noise") according to
 * the provided definition_noise.
 */
void
add_noise(contiguous_state<4,complex_t> &state,
          const specification_noise& noise,
          const definition_scenario& scenario,
          const specification_grid& grid,
          const pencil_grid& dgrid,
          const bsplineop& cop,
          bspline &b);

/**
 * Using the provided state, sample the quantities declared in \ref samples
 * with the notable exceptions of those listed in \ref
 * SUZERAIN_SAMPLES_IMPLICIT.  This is an expensive, collective method.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     sg       Slow growth specification in use.
 * @param[in]     otool    Operator definitions in use.
 * @param[in]     b        B-spline basis for wall-normal direction.
 * @param[in,out] swave    Destroyed in the computation
 * @param[in]     t        Current simulation time
 *
 * @return Mean quantities as B-spline coefficients.
 */
std::auto_ptr<samples>
take_samples(const definition_scenario &scenario,
             const specification_largo &sg,
             const operator_tools &otool,
             const bspline &b,
             contiguous_state<4,complex_t> &swave,
             const real_t t);

/**
 * Using the provided state, sample the mean solution profiles declared
 * in \ref profile.  This is a mildly expensive, collective method.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     otool    Operator definitions in use.
 * @param[in,out] swave    Destroyed in the computation
 *
 * @return Mean quantity profiles as B-spline coefficients.
 */
std::auto_ptr<profile>
take_profile(const definition_scenario &scenario,
             const operator_tools &otool,
             contiguous_state<4,complex_t> &swave);

/**
 * Using the provided state, find the mean reference quantities declared in \ref
 * references.  This is a mildly expensive, collective method.
 *
 * @param[in]  scenario Scenario parameters
 * @param[in]  grid     Grid definitions in use
 * @param[in]  dgrid    Parallel, dealiased grid in use.
 * @param[in]  sphys    Physical space state on collocation points.
 * @param[out] refs     To be populated during the method invocation.
 */
void
collect_references(const definition_scenario &scenario,
                   const specification_grid& grid,
                   const pencil_grid &dgrid,
                   const physical_view<5> &sphys,
                   references &refs);

/**
 * Using the provided state, find the mean instantaneous quantities declared in
 * \ref instantaneous.  This is a mildly expensive, collective method.
 *
 * @param[in]  scenario Scenario parameters
 * @param[in]  grid     Grid definitions in use
 * @param[in]  dgrid    Parallel, dealiased grid in use.
 * @param[in]  sphys    Physical space state on collocation points.
 * @param[out] inst     To be populated during the method invocation.
 */
void
collect_instantaneous(const definition_scenario &scenario,
                      const specification_grid &grid,
                      const pencil_grid &dgrid,
                      const physical_view<5> &sphys,
                      instantaneous &inst);

/**
 * Use the boundary layer information in \c prof and possibly base flow
 * information in \c sg to compute many quantities of interest.  This is
 * a purely local computation requiring no communication.
 *
 * @param[in]  prof     Profile information from \ref sample_profile().
 * @param[in]  scenario Scenario of interest.
 * @param[in]  sg       Slow growth definition optionally in use
 *                      which provides base flow details for
 *                      streamwise pressure and velocity gradients.
 * @param[in]  masslu   A factored mass matrix corresponding to \c b.
 * @param[in]  b        The B-spline basis in use.
 * @param[out] wall     Populated on return.
 * @param[out] viscous  Populated on return.
 * @param[out] thick    Populated on return.
 * @param[out] edge     Populated on return.
 * @param[out] edge99   Populated on return.
 * @param[out] reynolds Populated on return.
 * @param[out] qoi      Populated on return.
 * @param[out] pg       Populated on return.
 */
void summarize_boundary_layer_nature(
        const profile &prof,
        const definition_scenario &scenario,
        const shared_ptr<specification_largo> &sg,
        const bsplineop_lu &masslu,
        bspline &b,
        suzerain_bl_local       &wall,
        suzerain_bl_viscous     &viscous,
        suzerain_bl_thicknesses &thick,
        suzerain_bl_local       &edge,
        suzerain_bl_local       &edge99,
        suzerain_bl_reynolds    &reynolds,
        suzerain_bl_qoi         &qoi,
        suzerain_bl_pg          &pg);

/**
 * Use the boundary layer information in \c prof to compute many quantities of
 * interest.  This is a purely local computation requiring no communication.
 *
 * @param[in]  prof     Profile information from \ref sample_profile().
 * @param[in]  scenario Scenario of interest.
 * @param[in]  b        The B-spline basis in use.
 * @param[out] wall     Populated on return.
 * @param[out] viscous  Populated on return.
 * @param[out] center   Populated on return.
 * @param[out] qoi      Populated on return.
 */
void summarize_channel_nature(
        const profile &prof,
        const definition_scenario &scenario,
        bspline &b,
        suzerain_channel_local   &wall,
        suzerain_channel_viscous &viscous,
        suzerain_channel_local   &center,
        suzerain_channel_qoi     &qoi);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_HPP
