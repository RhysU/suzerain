//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// explicit_n.hpp: Nonlinear operators for channel_explicit
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef EXPLICIT_N_HPP
#define EXPLICIT_N_HPP

#include <suzerain/bspline_operators.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state.hpp>

#include "precision.hpp"
#include "channel.hpp"

#pragma warning(disable:383 1572)

namespace channel {

class NonlinearOperator
    : public suzerain::timestepper::INonlinearOperator<
            suzerain::NoninterleavedState<4,complex_t>
      >
{

public:
    typedef suzerain::NoninterleavedState<4,complex_t> state_type;

    NonlinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            const suzerain::bsplineop_luz &bopluz);

    virtual real_t applyOperator(
            suzerain::NoninterleavedState<4,complex_t> &state,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const bool delta_t_requested = false) const;

protected:

    const suzerain::problem::ScenarioDefinition<real_t> &scenario;
    const suzerain::problem::GridDefinition &grid;
    const suzerain::pencil_grid &dgrid;
    const suzerain::bsplineop &bop;
    const suzerain::bsplineop_luz &bopluz;

    Eigen::ArrayXr one_over_delta_y;

    // Details for suzerain::diffwave::* calls
    const int Ny;
    const int Nx;
    const int dNx;
    const int dkbx;
    const int dkex;
    const int Nz;
    const int dNz;
    const int dkbz;
    const int dkez;

private:
    mutable suzerain::pencil<> rho_x, rho_y, rho_z;
    mutable suzerain::pencil<> rho_xx, rho_xy, rho_xz, rho_yy, rho_yz, rho_zz;

    mutable suzerain::pencil<> mx_x, mx_y, mx_z;
    mutable suzerain::pencil<> mx_xx, mx_xy, mx_xz, mx_yy, mx_yz, mx_zz;

    mutable suzerain::pencil<> my_x, my_y, my_z;
    mutable suzerain::pencil<> my_xx, my_xy, my_xz, my_yy, my_yz, my_zz;

    mutable suzerain::pencil<> mz_x, mz_y, mz_z;
    mutable suzerain::pencil<> mz_xx, mz_xy, mz_xz, mz_yy, mz_yz, mz_zz;

    mutable suzerain::pencil<> e_x, e_y, e_z, div_grad_e;
};

class NonlinearOperatorWithBoundaryConditions : public NonlinearOperator
{

public:

    typedef NonlinearOperator base;

    NonlinearOperatorWithBoundaryConditions(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            const suzerain::bsplineop_luz &bopluz);

    virtual real_t applyOperator(
            suzerain::NoninterleavedState<4,complex_t> &state,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const bool delta_t_requested = false) const;

protected:
    Eigen::ArrayXr bintcoeff;

};

} // namespace channel

#endif  /* EXPLICIT_N_HPP */
