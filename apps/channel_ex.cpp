//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
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
// channel_ex.cpp: A fully explicit channel calculation using Suzerain
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <esio/esio.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/math.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/utility.hpp>

#include "logger.hpp"
#include "channel_common.hpp"

#pragma warning(disable:383 1572)

using boost::make_shared;
using boost::numeric_cast;
using boost::shared_ptr;
using std::numeric_limits;

// Explicit timestepping scheme uses only complex_t 4D NoninterleavedState
// State indices range over (scalar field, Y, X, Z) in wave space
// Establish some shorthand for some overly-templated operator and state types
typedef suzerain::storage::noninterleaved<4> storage_type;
typedef suzerain::IState<
            storage_type::dimensionality, complex_t, storage_type
        > istate_type;
typedef suzerain::timestepper::lowstorage::ILinearOperator<
            storage_type::dimensionality, complex_t, storage_type
        > ilinearoperator_type;
typedef suzerain::timestepper::INonlinearOperator<
            storage_type::dimensionality, complex_t, storage_type
        > inonlinearoperator_type;
typedef suzerain::NoninterleavedState<
            storage_type::dimensionality, complex_t
        > state_type;

// Global scenario parameters initialized in main().  These are declared const
// to avoid accidental modification but have their const-ness const_cast away
// where necessary to load settings.
using suzerain::problem::ScenarioDefinition;
using suzerain::problem::GridDefinition;
using suzerain::problem::RestartDefinition;
using suzerain::problem::TimeDefinition;
static const ScenarioDefinition<real_t> scenario(0, 0, 0, 0, 0, 0, 0);
static const GridDefinition<real_t> grid(0, 0, 0, 0, 0, 0);
static const RestartDefinition<> restart(/* load         */ "",
                                         /* metadata     */ "metadata.h5",
                                         /* uncommitted  */ "uncommitted.h5",
                                         /* desttemplate */ "restart#.h5",
                                         /* retain       */ 1,
                                         /* every_dt     */ 0,
                                         /* every_nt     */ 100);
static const TimeDefinition<real_t> timedef(0, 1, 0, 1);

// Global grid-details initialized in main()
static shared_ptr<suzerain::bspline>     bspw;
static shared_ptr<suzerain::bspline_luz> bspluzw;
static shared_ptr<suzerain::pencil_grid> pg;

// State details specific to this rank initialized in main()
static shared_ptr<state_type> state_linear;
static shared_ptr<state_type> state_nonlinear;
static boost::array<suzerain::pencil_grid::index,3> state_start;
static boost::array<suzerain::pencil_grid::index,3> state_end;
static boost::array<suzerain::pencil_grid::index,3> state_extent;

// Bulk density is computed in main() and used within NonlinearOperator
static real_t bulk_density;

// TODO Incorporate IOperatorLifecycle semantics
// TODO Refactor MassOperator into templated BsplineMassOperator

class MassOperator : public ilinearoperator_type
{
public:
    explicit MassOperator(real_t scaling = 1)
        : opscaling_(scaling), luzw_(*bspw)
    {
        complex_t coefficient;
        suzerain::complex::assign_complex(coefficient, scaling);
        luzw_.form_general(1, &coefficient, *bspw);
    }

    virtual void applyMassPlusScaledOperator(
            const complex_t scale,
            istate_type &istate) const throw (std::exception)
    {
        SUZERAIN_UNUSED(scale);
        state_type &state = dynamic_cast<state_type&>(istate);

        const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
        assert(1 == state.strides()[1]);
        assert(static_cast<unsigned>(luzw_.ndof()) == state.shape()[1]);
        bspw->apply_operator(0, nrhs, opscaling_,
                state.memory_begin(), 1, state.strides()[2]);
    }

    virtual void accumulateMassPlusScaledOperator(
            const complex_t scale,
            const istate_type &iinput,
            istate_type &ioutput) const throw (std::exception)
    {
        SUZERAIN_UNUSED(scale);
        const state_type &x = dynamic_cast<const state_type&>(iinput);
        state_type &y = dynamic_cast<state_type&>(ioutput);
        assert(std::equal(x.shape(), x.shape() + state_type::dimensionality,
                          y.shape()));

        typedef state_type::index index;
        for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
            ix < static_cast<index>(x.index_bases()[0] + x.shape()[0]);
            ++ix, ++iy) {

            for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
                lx < static_cast<index>(x.index_bases()[3] + x.shape()[3]);
                ++lx, ++ly) {

                bspw->accumulate_operator(0, x.shape()[2], opscaling_,
                        &(x[ix][x.index_bases()[1]][x.index_bases()[2]][lx]),
                        x.strides()[1], x.strides()[2],
                        1.0, &(y[iy][y.index_bases()[1]][y.index_bases()[2]][ly]),
                        y.strides()[1], y.strides()[2]);
            }
        }
    }

    virtual void invertMassPlusScaledOperator(
            const complex_t scale,
            istate_type &istate) const throw (std::exception)
    {
        SUZERAIN_UNUSED(scale);
        state_type &state = dynamic_cast<state_type&>(istate);

        const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
        assert(1 == state.strides()[1]);
        assert(static_cast<unsigned>(luzw_.ndof()) == state.shape()[1]);
        luzw_.solve(nrhs, state.memory_begin(), 1, state.strides()[2]);
    }

private:
    const real_t opscaling_;
    suzerain::bspline_luz luzw_;
};

// TODO Incorporate IOperatorLifecycle semantics for NonlinearOperator

class NonlinearOperator : public inonlinearoperator_type
{

public:

    NonlinearOperator() :
        Ny(pg->global_extents()[1]),
        dNx(pg->global_extents()[0]),
        dkbx(pg->local_wave_start()[0]),
        dkex(pg->local_wave_end()[0]),
        dNz(pg->global_extents()[2]),
        dkbz(pg->local_wave_start()[2]),
        dkez(pg->local_wave_end()[2]),
        rho_x(*pg), rho_y(*pg), rho_z(*pg),
        rho_xx(*pg), rho_xy(*pg), rho_xz(*pg),
                     rho_yy(*pg), rho_yz(*pg),
                                  rho_zz(*pg),
        mx_x(*pg), mx_y(*pg), mx_z(*pg),
        mx_xx(*pg), mx_xy(*pg), mx_xz(*pg),
                    mx_yy(*pg), mx_yz(*pg),
                                mx_zz(*pg),
        my_x(*pg), my_y(*pg), my_z(*pg),
        my_xx(*pg), my_xy(*pg), my_xz(*pg),
                    my_yy(*pg), my_yz(*pg),
                                my_zz(*pg),
        mz_x(*pg), mz_y(*pg), mz_z(*pg),
        mz_xx(*pg), mz_xy(*pg), mz_xz(*pg),
                    mz_yy(*pg), mz_yz(*pg),
                                mz_zz(*pg),
        e_x(*pg), e_y(*pg), e_z(*pg), div_grad_e(*pg)
    {
        // NOP
    }

    real_t applyOperator(
        istate_type &istate,
        const bool delta_t_requested = false)
        const
        throw(std::exception) {

        SUZERAIN_UNUSED(delta_t_requested);
        real_t convective_delta_t = numeric_limits<real_t>::max();
        real_t diffusive_delta_t  = numeric_limits<real_t>::max();
        const real_t one_over_delta_x = scenario.Lx / grid.Nx;
        const real_t one_over_delta_y = scenario.Ly / grid.Ny;
        const real_t one_over_delta_z = scenario.Lz / grid.Nz;

        // Get state information with appropriate type
        state_type &state = dynamic_cast<state_type&>(istate);

        // Create 3D views of 4D state information
        using boost::array_view_gen;
        using boost::indices;
        using boost::multi_array_types::index_range;
        array_view_gen<state_type,3>::type state_rho
            = state[indices[0][index_range()][index_range()][index_range()]];
        array_view_gen<state_type,3>::type state_rhou
            = state[indices[1][index_range()][index_range()][index_range()]];
        array_view_gen<state_type,3>::type state_rhov
            = state[indices[2][index_range()][index_range()][index_range()]];
        array_view_gen<state_type,3>::type state_rhow
            = state[indices[3][index_range()][index_range()][index_range()]];
        array_view_gen<state_type,3>::type state_rhoe
            = state[indices[4][index_range()][index_range()][index_range()]];

        const real_t complex_one[2]  = { 1.0, 0.0 };
        const real_t complex_zero[2] = { 0.0, 0.0 };

        // All state enters routine as coefficients in X, Y, and Z directions

        // On "zero-zero" rank save a copy of the constant x-momentum modes
        // Need these later to apply the momentum forcing's energy contribution
        boost::scoped_array<real_t> original_state_mx;
        if (dkbx == 0 && dkbz == 0) {
            original_state_mx.reset(new real_t[state.shape()[1]]);
            for (std::size_t i = 0; i < state.shape()[1]; ++i) {
                original_state_mx[i]
                    = suzerain::complex::real(state_rhou[i][0][0]);
            }
        }

        // Compute Y derivatives of density at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, state_rho.origin(), 1, state.shape()[1],
            complex_zero, rho_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, state_rho.origin(), 1, state.shape()[1],
            complex_zero, rho_yy.wave.begin(), 1, state.shape()[1]);
        // Compute density at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, state_rho.origin(), 1, state.shape()[1]);

        // Compute X-related derivatives of density at collocation points
        suzerain::diffwave::accumulate(1, 0, // dx
            complex_one, state_rho.origin(),
            complex_zero, rho_x.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rho.origin(),
            complex_zero, rho_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rho.origin(),
            complex_zero, rho_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of density at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rho.origin(),
            complex_zero, rho_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rho.origin(),
            complex_zero, rho_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);


        // Compute Y derivatives of X momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhou.origin(), 1, state.shape()[1],
            complex_zero, mx_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhou.origin(), 1, state.shape()[1],
            complex_zero, mx_yy.wave.begin(), 1, state.shape()[1]);
        // Compute X momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, state_rhou.origin(), 1, state.shape()[1]);

        // Compute X-related derivatives of X momentum at collocation points
        suzerain::diffwave::accumulate(1, 0, // dx
            complex_one, state_rhou.origin(),
            complex_zero, mx_x.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhou.origin(),
            complex_zero, mx_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhou.origin(),
            complex_zero, mx_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of X momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhou.origin(),
            complex_zero, mx_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhou.origin(),
            complex_zero, mx_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);


        // Compute Y derivatives of Y momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhov.origin(), 1, state.shape()[1],
            complex_zero, my_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhov.origin(), 1, state.shape()[1],
            complex_zero, my_yy.wave.begin(), 1, state.shape()[1]);
        // Compute Y momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, state_rhov.origin(), 1, state.shape()[1]);

        // Compute X-related derivatives of Y momentum at collocation points
        suzerain::diffwave::accumulate(1, 0, // dx
            complex_one, state_rhov.origin(),
            complex_zero, my_x.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhov.origin(),
            complex_zero, my_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, my_y.wave.begin(),
            complex_zero, my_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhov.origin(),
            complex_zero, my_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of Y momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhov.origin(),
            complex_zero, my_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, my_y.wave.begin(),
            complex_zero, my_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhov.origin(),
            complex_zero, my_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);


        // Compute Y derivatives of Z momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhow.origin(), 1, state.shape()[1],
            complex_zero, mz_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhow.origin(), 1, state.shape()[1],
            complex_zero, mz_yy.wave.begin(), 1, state.shape()[1]);
        // Compute Y momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, state_rhow.origin(), 1, state.shape()[1]);

        // Compute X-related derivatives of Z momentum at collocation points
        suzerain::diffwave::accumulate(1, 0, // dx
            complex_one, state_rhow.origin(),
            complex_zero, mz_x.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhow.origin(),
            complex_zero, mz_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhow.origin(),
            complex_zero, mz_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of Z momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhow.origin(),
            complex_zero, mz_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhow.origin(),
            complex_zero, mz_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);


        // Compute Y derivatives of total energy at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhoe.origin(), 1, state.shape()[1],
            complex_zero, e_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, state_rhoe.origin(), 1, state.shape()[1],
            complex_zero, div_grad_e.wave.begin(), 1, state.shape()[1]);
        // Compute total energy at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, state_rhoe.origin(), 1, state.shape()[1]);

        // Compute X-related derivatives of total energy at collocation points
        suzerain::diffwave::accumulate(1, 0, // dx
            complex_one, state_rhoe.origin(),
            complex_zero, e_x.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhoe.origin(),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of total energy at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhoe.origin(),
            complex_zero, e_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhoe.origin(),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            scenario.Lx, scenario.Lz,
            Ny, grid.Nx, dNx, dkbx, dkex, grid.Nz, dNz, dkbz, dkez);

        // Collectively convert state to physical space
        pg->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rho.origin()));
        pg->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhou.origin()));
        pg->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhov.origin()));
        pg->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhow.origin()));
        pg->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhoe.origin()));

        // Collectively convert state derivatives to physical space
        pg->transform_wave_to_physical(rho_x.data());   // density
        pg->transform_wave_to_physical(rho_y.data());
        pg->transform_wave_to_physical(rho_z.data());
        pg->transform_wave_to_physical(rho_xx.data());
        pg->transform_wave_to_physical(rho_xy.data());
        pg->transform_wave_to_physical(rho_xz.data());
        pg->transform_wave_to_physical(rho_yy.data());
        pg->transform_wave_to_physical(rho_yz.data());
        pg->transform_wave_to_physical(rho_zz.data());
        pg->transform_wave_to_physical(mx_x.data());   // X momentum
        pg->transform_wave_to_physical(mx_y.data());
        pg->transform_wave_to_physical(mx_z.data());
        pg->transform_wave_to_physical(mx_xx.data());
        pg->transform_wave_to_physical(mx_xy.data());
        pg->transform_wave_to_physical(mx_xz.data());
        pg->transform_wave_to_physical(mx_yy.data());
        pg->transform_wave_to_physical(mx_yz.data());
        pg->transform_wave_to_physical(mx_zz.data());
        pg->transform_wave_to_physical(my_x.data());   // Y momentum
        pg->transform_wave_to_physical(my_y.data());
        pg->transform_wave_to_physical(my_z.data());
        pg->transform_wave_to_physical(my_xx.data());
        pg->transform_wave_to_physical(my_xy.data());
        pg->transform_wave_to_physical(my_xz.data());
        pg->transform_wave_to_physical(my_yy.data());
        pg->transform_wave_to_physical(my_yz.data());
        pg->transform_wave_to_physical(my_zz.data());
        pg->transform_wave_to_physical(mz_x.data());   // Z momentum
        pg->transform_wave_to_physical(mz_y.data());
        pg->transform_wave_to_physical(mz_z.data());
        pg->transform_wave_to_physical(mz_xx.data());
        pg->transform_wave_to_physical(mz_xy.data());
        pg->transform_wave_to_physical(mz_xz.data());
        pg->transform_wave_to_physical(mz_yy.data());
        pg->transform_wave_to_physical(mz_yz.data());
        pg->transform_wave_to_physical(mz_zz.data());
        pg->transform_wave_to_physical(e_x.data());   // Total energy
        pg->transform_wave_to_physical(e_y.data());
        pg->transform_wave_to_physical(e_z.data());
        pg->transform_wave_to_physical(div_grad_e.data());

        // Compute nonlinear operator

        // Retrieve constants and compute derived constants
        const real_t beta             = scenario.beta;
        const real_t gamma            = scenario.gamma;
        const real_t Pr               = scenario.Pr;
        const real_t Re               = scenario.Re;
        const real_t inv_Re           = 1 / Re;
        const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));

        // Temporary storage used within following loop
        Eigen::Vector3d grad_rho;
        Eigen::Matrix3d grad_grad_rho;
        Eigen::Vector3d m;
        Eigen::Matrix3d grad_m;
        Eigen::Vector3d div_grad_m;
        Eigen::Vector3d grad_div_m;
        Eigen::Vector3d grad_e;
        Eigen::Vector3d u;
        Eigen::Matrix3d grad_u;
        Eigen::Vector3d grad_div_u, div_grad_u;
        Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;
        Eigen::Matrix3d tau;
        Eigen::Vector3d div_tau;

        // Walk physical space state storage in linear fashion
        for (// Loop initialization
             real_t *p_rho    = reinterpret_cast<real_t *>(state_rho.origin()),
                    *p_rho_x  = rho_x.physical.begin(),
                    *p_rho_y  = rho_y.physical.begin(),
                    *p_rho_z  = rho_z.physical.begin(),
                    *p_rho_xx = rho_xx.physical.begin(),
                    *p_rho_xy = rho_xy.physical.begin(),
                    *p_rho_xz = rho_xz.physical.begin(),
                    *p_rho_yy = rho_yy.physical.begin(),
                    *p_rho_yz = rho_yz.physical.begin(),
                    *p_rho_zz = rho_zz.physical.begin(),
                    *p_mx     = reinterpret_cast<real_t *>(state_rhou.origin()),
                    *p_mx_x   = mx_x.physical.begin(),
                    *p_mx_y   = mx_y.physical.begin(),
                    *p_mx_z   = mx_z.physical.begin(),
                    *p_mx_xx  = mx_xx.physical.begin(),
                    *p_mx_xy  = mx_xy.physical.begin(),
                    *p_mx_xz  = mx_xz.physical.begin(),
                    *p_mx_yy  = mx_yy.physical.begin(),
                    *p_mx_yz  = mx_yz.physical.begin(),
                    *p_mx_zz  = mx_zz.physical.begin(),
                    *p_my     = reinterpret_cast<real_t *>(state_rhov.origin()),
                    *p_my_x   = my_x.physical.begin(),
                    *p_my_y   = my_y.physical.begin(),
                    *p_my_z   = my_z.physical.begin(),
                    *p_my_xx  = my_xx.physical.begin(),
                    *p_my_xy  = my_xy.physical.begin(),
                    *p_my_xz  = my_xz.physical.begin(),
                    *p_my_yy  = my_yy.physical.begin(),
                    *p_my_yz  = my_yz.physical.begin(),
                    *p_my_zz  = my_zz.physical.begin(),
                    *p_mz     = reinterpret_cast<real_t *>(state_rhow.origin()),
                    *p_mz_x   = mz_x.physical.begin(),
                    *p_mz_y   = mz_y.physical.begin(),
                    *p_mz_z   = mz_z.physical.begin(),
                    *p_mz_xx  = mz_xx.physical.begin(),
                    *p_mz_xy  = mz_xy.physical.begin(),
                    *p_mz_xz  = mz_xz.physical.begin(),
                    *p_mz_yy  = mz_yy.physical.begin(),
                    *p_mz_yz  = mz_yz.physical.begin(),
                    *p_mz_zz  = mz_zz.physical.begin(),
                    *p_e      = reinterpret_cast<real_t *>(state_rhoe.origin()),
                    *p_e_x    = e_x.physical.begin(),
                    *p_e_y    = e_y.physical.begin(),
                    *p_e_z    = e_z.physical.begin(),
                    *p_div_grad_e = div_grad_e.physical.begin();
             // Loop completion test
             p_rho_x != rho_x.physical.end();
             // Loop increment
             ++p_rho,
             ++p_rho_x,
             ++p_rho_y,
             ++p_rho_z,
             ++p_rho_xx,
             ++p_rho_xy,
             ++p_rho_xz,
             ++p_rho_yy,
             ++p_rho_yz,
             ++p_rho_zz,
             ++p_mx,
             ++p_mx_x,
             ++p_mx_y,
             ++p_mx_z,
             ++p_mx_xx,
             ++p_mx_xy,
             ++p_mx_xz,
             ++p_mx_yy,
             ++p_mx_yz,
             ++p_mx_zz,
             ++p_my,
             ++p_my_x,
             ++p_my_y,
             ++p_my_z,
             ++p_my_xx,
             ++p_my_xy,
             ++p_my_xz,
             ++p_my_yy,
             ++p_my_yz,
             ++p_my_zz,
             ++p_mz,
             ++p_mz_x,
             ++p_mz_y,
             ++p_mz_z,
             ++p_mz_xx,
             ++p_mz_xy,
             ++p_mz_xz,
             ++p_mz_yy,
             ++p_mz_yz,
             ++p_mz_zz,
             ++p_e,
             ++p_e_x,
             ++p_e_y,
             ++p_e_z,
             ++p_div_grad_e) {

            // Prepare local density-related quantities
            const real_t rho          = *p_rho;
            grad_rho[0]               = *p_rho_x;
            grad_rho[1]               = *p_rho_y;
            grad_rho[2]               = *p_rho_z;
            const real_t div_grad_rho = *p_rho_xx + *p_rho_yy + *p_rho_zz;
            grad_grad_rho(0,0)        = *p_rho_xx;
            grad_grad_rho(0,1)        = *p_rho_xy;
            grad_grad_rho(0,2)        = *p_rho_xz;
            grad_grad_rho(1,0)        = grad_grad_rho(0,1);
            grad_grad_rho(1,1)        = *p_rho_yy;
            grad_grad_rho(1,2)        = *p_rho_yz;
            grad_grad_rho(2,0)        = grad_grad_rho(0,2);
            grad_grad_rho(2,1)        = grad_grad_rho(1,2);
            grad_grad_rho(2,2)        = *p_rho_zz;

            // Prepare local momentum-related quantities
            m[0]               = *p_mx;
            m[1]               = *p_my;
            m[2]               = *p_mz;
            const real_t div_m = *p_mx_x + *p_my_y + *p_my_z;
            grad_m(0,0)        = *p_mx_x;
            grad_m(0,1)        = *p_mx_y;
            grad_m(0,2)        = *p_mx_z;
            grad_m(1,0)        = *p_my_x;
            grad_m(1,1)        = *p_my_y;
            grad_m(1,2)        = *p_my_z;
            grad_m(2,0)        = *p_mz_x;
            grad_m(2,1)        = *p_mz_y;
            grad_m(2,2)        = *p_mz_z;
            div_grad_m[0]      = *p_mx_xx + *p_mx_yy + *p_mx_zz;
            div_grad_m[1]      = *p_my_xx + *p_my_yy + *p_my_zz;
            div_grad_m[2]      = *p_mz_xx + *p_mz_yy + *p_mz_zz;
            grad_div_m[0]      = *p_mx_xx + *p_mx_xy + *p_mx_xz;
            grad_div_m[1]      = *p_mx_xy + *p_mx_yy + *p_mx_yz;
            grad_div_m[2]      = *p_mx_xz + *p_mx_yz + *p_mx_zz;

            // Prepare local total energy-related quantities
            const real_t e          = *p_e;
            grad_e[0]               = *p_e_x;
            grad_e[1]               = *p_e_y;
            grad_e[2]               = *p_e_z;
            const real_t div_grad_e = *p_div_grad_e;  // FIXME Shadow

            // Prepare quantities derived from local state and its derivatives
            u                  = suzerain::orthonormal::rhome::u(rho, m);
            const real_t div_u = suzerain::orthonormal::rhome::div_u(
                                    rho, grad_rho, m, div_m);
            grad_u             = suzerain::orthonormal::rhome::grad_u(
                                    rho, grad_rho, m, grad_m);
            grad_div_u         = suzerain::orthonormal::rhome::grad_div_u(
                                    rho, grad_rho, grad_grad_rho,
                                    m, div_m, grad_m, grad_div_m);
            div_grad_u         = suzerain::orthonormal::rhome::div_grad_u(
                                    rho, grad_rho, div_grad_rho,
                                    m, grad_m, div_grad_m);
            real_t p, T, mu, lambda;
            suzerain::orthonormal::rhome::p_T_mu_lambda(
                beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
                p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
            const real_t div_grad_p = suzerain::orthonormal::rhome::div_grad_p(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m,
                                        e, grad_e, div_grad_e);
            const real_t div_grad_T = suzerain::orthonormal::rhome::div_grad_T(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        p, grad_p, div_grad_p);
            tau     = suzerain::orthonormal::tau(
                        mu, lambda, div_u, grad_u);
            div_tau = suzerain::orthonormal::div_tau(
                        mu, grad_mu, lambda, grad_lambda,
                        div_u, grad_u, div_grad_u, grad_div_u);

            // Maintain the minimum stable time step
            // TODO Operator knows about the scheme's eigenvalues.  Fix that.
            convective_delta_t = std::min(convective_delta_t,
                suzerain::timestepper::convective_stability_criterion(
                    u.x(), one_over_delta_x,
                    u.y(), one_over_delta_y,
                    u.z(), one_over_delta_z,
                    std::sqrt(3.0),
                    std::sqrt(T))); // nondimensional a = sqrt(T)
            diffusive_delta_t = std::min(diffusive_delta_t,
                suzerain::timestepper::diffusive_stability_criterion(
                    one_over_delta_x,
                    one_over_delta_y,
                    one_over_delta_z,
                    Re, Pr, gamma, 2.512, mu / rho));

            // Continuity equation
            *p_rho = - div_m
                ;

            // Momentum equation
            Eigen::Vector3d momentum =
                - suzerain::orthonormal::div_u_outer_m(m, grad_m, u, div_u)
                - grad_p
                + inv_Re * div_tau
                ;
            *p_mx = momentum[0];
            *p_my = momentum[1];
            *p_mz = momentum[2];

            // Energy equation
            *p_e = - suzerain::orthonormal::div_e_u(
                        e, grad_e, u, div_u
                     )
                   - suzerain::orthonormal::div_p_u(
                        p, grad_p, u, div_u
                     )
                   + inv_Re_Pr_gamma1 * suzerain::orthonormal::div_mu_grad_T(
                        grad_T, div_grad_T, mu, grad_mu
                     )
                   + inv_Re * suzerain::orthonormal::div_tau_u<real_t>(
                        u, grad_u, tau, div_tau
                     )
                   ;
        }

        // Convert collocation point values to wave space
        pg->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rho.origin()));
        pg->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhou.origin()));
        pg->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhov.origin()));
        pg->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhow.origin()));
        pg->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhoe.origin()));

        // Convert collocation point values to Bspline coefficients
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                state_rho.origin(), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                state_rhou.origin(), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                state_rhov.origin(), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                state_rhow.origin(), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                state_rhoe.origin(), 1, state.shape()[1]);

        // ------------------------------------------------------------------
        // BEGIN: Boundary conditions and driving forces
        // See writeup/channel_treatment.tex for full details
        // ------------------------------------------------------------------
        const std::size_t lower_wall = 0;                    // index of wall
        const std::size_t upper_wall = state.shape()[1] - 1; // index of wall

        // Add f_rho to mean density right hand side per writeup step (2)
        if (dkbx == 0 && dkbz == 0) {                // "zero-zero" rank only
            real_t f_rho;
            bspw->integrate(reinterpret_cast<real_t *>(state_rho.origin()),
                            sizeof(complex_t)/sizeof(real_t),
                            &f_rho);
            f_rho /= scenario.Ly;
            for (std::size_t i = 0; i < state.shape()[1]; ++i) {
                state_rho[i][0][0] -= f_rho;
            }
        }

        // Set no-slip condition on walls per writeup step (3)
        // Done as three separate loops to walk memory linearly
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {      // x momentum
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                state_rhou[lower_wall][j][k] = 0;
                state_rhou[upper_wall][j][k] = 0;
            }
        }
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {      // y momentum
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                state_rhov[lower_wall][j][k] = 0;
                state_rhov[upper_wall][j][k] = 0;
            }
        }
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {      // z momentum
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                state_rhow[lower_wall][j][k] = 0;
                state_rhow[upper_wall][j][k] = 0;
            }
        }

        // Set isothermal condition on walls per writeup step (4)
        const real_t inv_gamma_gamma1 = 1 / (gamma * (gamma - 1));
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                state_rhoe[lower_wall][j][k]
                    = inv_gamma_gamma1 * state_rho[lower_wall][j][k];
                state_rhoe[upper_wall][j][k]
                    = inv_gamma_gamma1 * state_rho[upper_wall][j][k];
            }
        }

        // Apply f_{m_x} term to mean x-momentum, mean energy
        if (dkbx == 0 && dkbz == 0) {                // "zero-zero" rank only

            // Compute temporary per writeup implementation step (5)
            real_t alpha;
            bspw->integrate(reinterpret_cast<real_t *>(state_rhou.origin()),
                            sizeof(complex_t)/sizeof(real_t),
                            &alpha);
            alpha /= scenario.Ly;

            // Apply to non-wall mean x-momentum right hand side per step (6)
            for (std::size_t i = lower_wall + 1; i < upper_wall; ++i) {
                state_rhou[i][0][0] -= alpha;
            }

            // Apply to non-wall mean energy right hand side per step (7)
            alpha /= bulk_density;
            for (std::size_t i = lower_wall + 1; i < upper_wall; ++i) {
                state_rhoe[i][0][0] -= alpha * original_state_mx[i];
            }
        }

        // ------------------------------------------------------------------
        // END: Boundary conditions and driving forces
        // See writeup/channel_treatment.tex for full details
        // ------------------------------------------------------------------

        return std::min(convective_delta_t, diffusive_delta_t);
    }

private:
    // Details for suzerain::diffwave::* calls
    const int Ny;
    const int dNx;
    const int dkbx;
    const int dkex;
    const int dNz;
    const int dkbz;
    const int dkez;

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

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

/** Routine to load state from file.  */
static void load_state(esio_handle h, state_type &state)
{
    assert(state.shape()[0] == field_names.static_size);

    // TODO Load state with different Nx
    // TODO Load state with different Nz
    // TODO Load state with different Ny
    int zglobal, xglobal, yglobal, ncomponents;
    esio_field_sizev(
            h, field_names[0], &zglobal, &xglobal, &yglobal, &ncomponents);
    assert(static_cast<unsigned>(zglobal) == grid.Nz);
    assert(static_cast<unsigned>(xglobal) == grid.Nx);
    assert(static_cast<unsigned>(yglobal) == grid.Ny);
    assert(ncomponents == 2);

    esio_field_establish(h, zglobal, state_start[2], state_extent[2],
                            xglobal, state_start[0], state_extent[0],
                            yglobal, state_start[1], state_extent[1]);

    for (size_t i = 0; i < field_names.static_size; ++i) {
        complex_field_read(h, field_names[i], state[i].origin(),
                state.strides()[3], state.strides()[2], state.strides()[1]);
    }
}

/** Routine to output status.  Signature for TimeController use. */
static bool log_status(real_t t, std::size_t nt) {
    INFO("Simulation reached time " << t << " at time step " << nt);
    return true;
}

static std::size_t last_restart_saved_nt = 0;

/** Routine to store a restart file.  Signature for TimeController use. */
static bool save_restart(real_t t, std::size_t nt)
{
    esio_file_clone(esioh, restart.metadata().c_str(),
                    restart.uncommitted().c_str(), 1 /*overwrite*/);

    // Save simulation time information
    store_time(esioh, t);

    // TODO Save only non-dealiased portion of state
    DEBUG("Storing simulation fields at simulation step " << nt);
    esio_field_establish(esioh, grid.Nz, state_start[2], state_extent[2],
                                grid.Nx, state_start[0], state_extent[0],
                                grid.Ny, state_start[1], state_extent[1]);

    for (size_t i = 0; i < field_names.static_size; ++i) {
        complex_field_write(esioh,
                field_names[i],
                (*state_linear)[i].origin(),
                (*state_linear).strides()[3],
                (*state_linear).strides()[2],
                (*state_linear).strides()[1],
                field_descriptions[i]);
    }

    esio_file_close_restart(esioh,
                            restart.desttemplate().c_str(),
                            restart.retain());

    INFO("Successfully wrote restart file at t = " << t << " for nt = " << nt);

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return true; // Continue
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Obtain some basic MPI environment details.
    const int nranks = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Establish MPI-savvy, rank-dependent logging names
    name_logger_within_comm_world();

    INFO("Processing command line arguments and response files");
    {
        suzerain::ProgramOptions options(
                "Suzerain-based explicit compressible channel simulation");
        options.add_definition(
                const_cast<ScenarioDefinition<real_t>&>(scenario));
        options.add_definition(
                const_cast<GridDefinition<real_t>&>(grid));
        options.add_definition(
                const_cast<RestartDefinition<>&>(restart));
        options.add_definition(
                const_cast<TimeDefinition<real_t>&>(timedef));
        options.process(argc, argv);
    }

    // TODO Account for grid differences at load time

    INFO("Loading details from restart file: " << restart.load());
    esio_file_open(esioh, restart.load().c_str(), 0 /* read-only */);
    load(esioh, const_cast<ScenarioDefinition<real_t>&>(scenario));
    load(esioh, const_cast<GridDefinition<real_t>&>(grid));
    load(esioh, bspw, const_cast<GridDefinition<real_t>&>(grid));
    esio_file_close(esioh);

    INFO("Saving metadata template file: " << restart.metadata());
    {
        esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_create(h, restart.metadata().c_str(), 1 /* overwrite */);
        store(h, scenario);
        store(h, grid, scenario.Lx, scenario.Lz);
        store(h, bspw);
        esio_file_close(h);
        esio_handle_finalize(h);
    }

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<suzerain::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    pg = make_shared<suzerain::pencil_grid>(grid.dealiased_extents(),
                                            grid.processor_grid);
    INFO( "Number of MPI ranks:               " << nranks);
    INFO( "Rank grid used for decomposition:  " << pg->processor_grid());
    DEBUG("Local dealiased wave start  (XYZ): " << pg->local_wave_start());
    DEBUG("Local dealiased wave end    (XYZ): " << pg->local_wave_end());
    DEBUG("Local dealiased wave extent (XYZ): " << pg->local_wave_extent());

    // Compute how much non-dealiased XYZ state is local to this rank
    // Additional munging necessary X direction has (Nx/2+1) complex values
    state_start = pg->local_wave_start();
    state_end[0] = std::min<suzerain::pencil_grid::index>(
            grid.global_extents[0]/2+1, pg->local_wave_end()[0]);
    for (int i = 1; i < state_end.static_size; ++i) {
        state_end[i] = std::min<suzerain::pencil_grid::index>(
                grid.global_extents[i], pg->local_wave_end()[i]);
    }
    for (int i = 0; i < state_extent.static_size; ++i) {
        state_extent[i] = std::max<suzerain::pencil_grid::index>(
                state_end[i] - state_start[i], 0);
    }

    DEBUG("Local state wave start  (XYZ): " << state_start);
    DEBUG("Local state wave end    (XYZ): " << state_end);
    DEBUG("Local state wave extent (XYZ): " << state_extent);

    // Create the state storage for the linear and nonlinear operators
    // with appropriate padding to allow nonlinear state to be P3DFFTified
    {
        using suzerain::to_yxz;
        using suzerain::prepend;
        using suzerain::strides_cm;
        state_linear.reset(new state_type(to_yxz(5, state_extent)));
        state_nonlinear.reset(new state_type(
                to_yxz(5, state_extent),
                prepend(pg->local_wave_storage(),
                        strides_cm(to_yxz(pg->local_wave_extent())))));
    }
    if (DEBUG_ENABLED) {
        boost::array<suzerain::pencil_grid::index,4> strides;
        std::copy(state_linear->strides(),
                  state_linear->strides() + 4, strides.begin());
        DEBUG("Linear state strides    (FYXZ): " << strides);
        std::copy(state_nonlinear->strides(),
                  state_nonlinear->strides() + 4, strides.begin());
        DEBUG("Nonlinear state strides (FYXZ): " << strides);
    }

    // Zero out any garbage in state_{non,}linear
    // Use fill rather than state_{non,}linear.scale to wipe any NaNs
    suzerain::multi_array::fill(*state_linear, 0);
    suzerain::multi_array::fill(*state_nonlinear, 0);

    // Load restart information into state_linear, including simulation time
    esio_file_open(esioh, restart.load().c_str(), 0 /* read-only */);
    real_t initial_t;
    load_time(esioh, initial_t);
    load_state(esioh, *state_linear);
    esio_file_close(esioh);

    // Compute bulk density, which we hold constant in time, from the state
    bspw->integrate(reinterpret_cast<real_t *>((*state_linear)[0].origin()),
                    sizeof(complex_t)/sizeof(real_t),
                    &bulk_density);
    DEBUG("Bulk density will be held constant at " << bulk_density);

    // Instantiate the operators and time stepping details
    // See write up section 2.1 (Spatial Discretization) for coefficient origin
    const suzerain::timestepper::lowstorage::SMR91Method<complex_t> smr91;
    MassOperator L(scenario.Lx * scenario.Lz * grid.Nx * grid.Nz);
    NonlinearOperator N;

    // Establish TimeController for use with operators and state storage
    using suzerain::timestepper::TimeController;
    boost::scoped_ptr<TimeController<real_t> > tc(
            make_LowStorageTimeController(
                smr91, L, N, *state_linear, *state_nonlinear, initial_t));

    // Register status callbacks status_{dt,nt}, if requested
    if (timedef.status_dt || timedef.status_nt) {
        tc->add_periodic_callback(
                timedef.status_dt ? timedef.status_dt : tc->forever_t(),
                timedef.status_nt ? timedef.status_nt : tc->forever_nt(),
                &log_status);
    }

    // Register restart-writing callbacks every_{dt,nt}, if requested
    if (restart.every_dt() || restart.every_nt()) {
        tc->add_periodic_callback(
                restart.every_dt() ? restart.every_dt() : tc->forever_t(),
                restart.every_nt() ? restart.every_nt() : tc->forever_nt(),
                &save_restart);
    }

    // Advance time according to advance_dt, advance_nt criteria
    switch ((!!timedef.advance_dt << 1) + !!timedef.advance_nt) {
        case 0:
            INFO("Advancing simulation until forcibly terminated");
            tc->advance();
            break;
        case 1:
            INFO("Advancing simulation " << timedef.advance_nt
                 << " discrete time steps");
            tc->step(timedef.advance_nt);
            break;
        case 2:
            INFO("Advancing simulation by " << timedef.advance_dt
                 << " units of physical time");
            tc->advance(initial_t + timedef.advance_dt);
            break;
        case 3:
            INFO("Advancing simulation by at most " << timedef.advance_dt
                 << " units of physical time");
            INFO("Advancing simulation by at most " << timedef.advance_nt
                 << " discrete time steps");
            tc->advance(initial_t + timedef.advance_dt,
                        timedef.advance_nt);
            break;
        default:
            FATAL("Sanity error in time control");
            return EXIT_FAILURE;
    }

    // Output statistics on time advancement
    INFO("Advanced simulation from t_initial = " << initial_t
         << " to t_final = " << tc->current_t()
         << " in " << tc->current_nt() << " steps");
    INFO("Min/mean/max/standard deviation of delta_t: "
         << tc->taken_min()  << ", "
         << tc->taken_mean() << ", "
         << tc->taken_max()  << ", "
         << tc->taken_stddev());

    // Save a final restart before exit if one was not just saved
    if (last_restart_saved_nt != tc->current_nt()) {
        INFO("Saving final restart file prior to quitting.");
        save_restart(tc->current_t(), tc->current_nt());
    }
}
