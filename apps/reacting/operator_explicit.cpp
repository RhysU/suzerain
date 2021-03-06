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
 * @copydoc operator_explicit.hpp
 */

#include "operator_explicit.hpp"

#include <suzerain/common.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>

#include "operator_nonlinear.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace reacting {

operator_mass_isothermal::operator_mass_isothermal(
        const antioch_constitutive& cmods,
        const specification_isothermal &isospec,
        const definition_channel &chdef,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : suzerain::operator_mass_isothermal(isospec, grid, dgrid, cop, b)
    , cmods(cmods)
    , chdef(chdef)
    , common(common)
    , who("operator.L")
{
}

explicit_operator_nonlinear::explicit_operator_nonlinear(
        const antioch_constitutive& cmods,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const definition_filter &fsdef,
        const specification_largo &sgdef,
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

std::vector<real_t> explicit_operator_nonlinear::apply_operator(
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

real_t operator_mass_isothermal::lower_E(
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

real_t operator_mass_isothermal::upper_E(
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
