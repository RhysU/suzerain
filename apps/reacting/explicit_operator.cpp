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
 * @copydoc explicit_operator.hpp
 */

#include "explicit_operator.hpp"

#include <suzerain/common.hpp>
#include <suzerain/largo_specification.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/state.hpp>

#include "nonlinear_operator.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace reacting {

isothermal_mass_operator::isothermal_mass_operator(
        const antioch_constitutive& cmods,
        const isothermal_specification &isospec,
        const channel_definition &chdef,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : suzerain::isothermal_mass_operator(isospec, grid, dgrid, cop, b)
    , cmods(cmods)
    , chdef(chdef)
    , common(common)
    , who("operator.L")
{
    // NOP
}

explicit_nonlinear_operator::explicit_nonlinear_operator(
        const antioch_constitutive& cmods,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const filter_definition &fsdef,
        const largo_specification &sgdef,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , cmods(cmods)
    , common(common)
    , msoln(msoln)
    , fsdef(fsdef)
    , sgdef(sgdef)
    , who("operator.N")
{
    // Ensure cached mass matrix factorized prior to first use
    // Strictly speaking unnecessary, but reduces timing variability
    this->masslu();
    this->massluz();
}

std::vector<real_t> explicit_nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{

#define ARGUMENTS *this, common, fsdef, sgdef, msoln, cmods, *massluz(),      \
                  time, swave, method

    // Dispatch to an optimized implementation depending on case:
    switch (common.filter_treatment) {

    case filter::none:

        switch (common.linearization) {
        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::none, filter::none>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::none, filter::none>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_y, filter::none>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_y, filter::none>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }
        break;

    case filter::cook:

        switch (common.linearization) {
        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::none, filter::cook>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::none, filter::cook>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_y, filter::cook>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_y, filter::cook>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }
        break;

    case filter::viscous:

        switch (common.linearization) {
        case linearize::none:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::none, filter::viscous>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::none, filter::viscous>(ARGUMENTS);

        case linearize::rhome_y:
            return (substep_index == 0)
                 ? apply_navier_stokes_spatial_operator<true,
                        linearize::rhome_y, filter::viscous>(ARGUMENTS)
                 : apply_navier_stokes_spatial_operator<false,
                        linearize::rhome_y, filter::viscous>(ARGUMENTS);

        default:
            SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
            break;
        }
        break;

    default:
        SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
        break;

    }


#undef ARGUMENTS
}

real_t isothermal_mass_operator::lower_E(
        const real_t lower_T,
        const real_t lower_u,
        const real_t lower_v,
        const real_t lower_w,
        const std::vector<real_t> lower_cs) const
{
    // FIXME Gather and use instantaneous averages when isnan(lower_[uvw]).

    // compute internal and kinetic energies
    const real_t E_internal = cmods.e_from_T(lower_T, lower_cs);
    const real_t E_kinetic  = 0.5 * (  lower_u*lower_u
                                     + lower_v*lower_v
                                     + lower_w*lower_w);

    return E_internal + E_kinetic;
}

real_t isothermal_mass_operator::upper_E(
        const real_t upper_T,
        const real_t upper_u,
        const real_t upper_v,
        const real_t upper_w,
        const std::vector<real_t> upper_cs) const
{
    // FIXME Gather and use instantaneous averages when isnan(lower_[uvw]).
    return lower_E(upper_T, upper_u, upper_v, upper_w, upper_cs);
}



} // namespace reacting

} // namespace suzerain
