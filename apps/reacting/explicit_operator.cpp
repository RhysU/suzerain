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
#include <suzerain/ndx.hpp>
#include <suzerain/state.hpp>

#include "nonlinear_operator.hpp"

#pragma warning(disable:383 1572)

namespace suzerain {

namespace reacting {

// A helper class for implementing isothermal, no-slip boundary conditions
class IsothermalNoSlipFunctor
{
private:
    const ptrdiff_t field_stride;
    const real_t    e_tot;
    const std::vector<real_t>& wall_mass_fractions;

public:
    IsothermalNoSlipFunctor(ptrdiff_t field_stride,
                            const real_t e_tot,
                            const std::vector<real_t>& wall_mass_fractions)
        : field_stride(field_stride)
        , e_tot(e_tot)
        , wall_mass_fractions(wall_mass_fractions)
    {}

    void operator()(complex_t &rho) const
    {
        // TODO Reorder to linearly access e, mx, my, mz, species

        for (size_t s=1; s<wall_mass_fractions.size(); ++s) {
            (&rho)[s*field_stride] = rho*wall_mass_fractions[s];
        }

        (&rho)[(ndx::mx - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::my - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::mz - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::e  - ndx::rho)*field_stride] = rho*e_tot;
    }
};

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

void isothermal_mass_operator::invert_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::method_interface<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using boost::indices;
    typedef boost::multi_array_types::index_range range;

    // Indexes only the first and last collocation point
    const std::size_t Ny         = state.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;
    range walls(wall_lower, wall_upper + 1, wall_upper - wall_lower);

    // Prepare a state view of density locations at lower and upper walls
    multi_array::ref<complex_t,4>::array_view<3>::type state_view
            = state[indices[ndx::rho][walls][range()][range()]];

    // Prepare functor setting pointwise BCs given density locations
    const IsothermalNoSlipFunctor bc_functor(
        state.strides()[0],
        cmods.e_from_T(chdef.T_wall, chdef.wall_mass_fractions),
        chdef.wall_mass_fractions);

    // Apply the functor to all wall-only density locations
    multi_array::for_each(state_view, bc_functor);

    // Apply boundary conditions to any requested constraint problems
    if (ic0) {
        multi_array::ref<complex_t,4>::array_view<3>::type ic0_view
                = (*ic0)[indices[ndx::rho][walls][range()][range()]];
        SUZERAIN_ENSURE(state.strides()[0] == ic0->strides()[0]); // NB!
        multi_array::for_each(ic0_view, bc_functor);
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // State leaves method as coefficients in X, Y, and Z directions
}

explicit_nonlinear_operator::explicit_nonlinear_operator(
        const antioch_constitutive& cmods,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common,
        const filter_definition &fsdef,
        const shared_ptr<const manufactured_solution>& msoln)
    : operator_base(grid, dgrid, cop, b)
    , cmods(cmods)
    , common(common)
    , msoln(msoln)
    , fsdef(fsdef)
    , who("operator.N")
{
    // Ensure cached mass matrix factorized prior to first use
    // Strictly speaking unnecessary, but reduces timing variability
    this->massluz();
}

std::vector<real_t> explicit_nonlinear_operator::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{

#define ARGUMENTS *this, common, fsdef, msoln, cmods, *massluz(), \
                  time, swave, evmaxmag_real, evmaxmag_imag

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

    // The diluter mass fraction is not present in lower_cs, 
    // and it's required by antioch to compute total (internal) energy.
    // Compute the diluter mass fraction in cs[0]
    const unsigned int Ns = lower_cs.size()+1;
    std::vector<real_t>  cs; // species mass fractions
    cs.resize(Ns);
    cs[0] = 1.0;
    for (unsigned int is=1; is<Ns; ++is) {
        cs[is]  = lower_cs[is-1];
        cs[0 ] -= cs[is];
    }

    // compute internal and kinetic energies
    const real_t E_internal = cmods.e_from_T(lower_T, cs);
//     std::cout << "E_internal = " << E_internal << ", " << 717.5 * lower_T << std::endl;
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
