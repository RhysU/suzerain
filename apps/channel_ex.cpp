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

#define EIGEN_DEFAULT_IO_FORMAT \
        Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ";", "", "", "[", "]")

#include <suzerain/common.hpp>
#include <esio/esio.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/grid_definition.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/inorder.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
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
using boost::scoped_ptr;
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
static const ScenarioDefinition<real_t> scenario;
static const GridDefinition grid;
static const RestartDefinition<> restart(/* load         */ "",
                                         /* metadata     */ "metadata.h5",
                                         /* uncommitted  */ "uncommitted.h5",
                                         /* desttemplate */ "restart#.h5",
                                         /* retain       */ 1,
                                         /* every_dt     */ 0,
                                         /* every_nt     */ 0);
static const TimeDefinition<real_t> timedef(/* advance_dt           */ 0,
                                            /* advance_nt           */ 0,
                                            /* status_dt            */ 0,
                                            /* status_nt            */ 0,
                                            /* min_dt               */ 0,
                                            /* max_dt               */ 0,
                                            /* evmagfactor per Prem */ 0.72);

// Global grid-details initialized in main()
static shared_ptr<const suzerain::bspline>     bspw;
static shared_ptr<      suzerain::bspline_luz> bspluzw;
static shared_ptr<const suzerain::pencil_grid> dgrid;
static Eigen::ArrayXr                          one_over_delta_y;

// State details specific to this rank initialized in main()
static shared_ptr<state_type> state_linear;
static shared_ptr<state_type> state_nonlinear;

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
        Ny(dgrid->global_wave_extent.y()),
        Nx(grid.N.x()),
        dNx(dgrid->global_wave_extent.x()),
        dkbx(dgrid->local_wave_start.x()),
        dkex(dgrid->local_wave_end.x()),
        Nz(grid.N.z()),
        dNz(dgrid->global_wave_extent.z()),
        dkbz(dgrid->local_wave_start.z()),
        dkez(dgrid->local_wave_end.z()),
        rho_x(*dgrid),  rho_y(*dgrid),  rho_z(*dgrid),
        rho_xx(*dgrid), rho_xy(*dgrid), rho_xz(*dgrid),
                        rho_yy(*dgrid), rho_yz(*dgrid),
                                        rho_zz(*dgrid),
        mx_x(*dgrid),  mx_y(*dgrid),  mx_z(*dgrid),
        mx_xx(*dgrid), mx_xy(*dgrid), mx_xz(*dgrid),
                       mx_yy(*dgrid), mx_yz(*dgrid),
                                      mx_zz(*dgrid),
        my_x(*dgrid),  my_y(*dgrid),  my_z(*dgrid),
        my_xx(*dgrid), my_xy(*dgrid), my_xz(*dgrid),
                       my_yy(*dgrid), my_yz(*dgrid),
                                      my_zz(*dgrid),
        mz_x(*dgrid),  mz_y(*dgrid),  mz_z(*dgrid),
        mz_xx(*dgrid), mz_xy(*dgrid), mz_xz(*dgrid),
                       mz_yy(*dgrid), mz_yz(*dgrid),
                                      mz_zz(*dgrid),
        e_x(*dgrid), e_y(*dgrid), e_z(*dgrid), div_grad_e(*dgrid)
    {
        // NOP
    }

    real_t applyOperator(
        istate_type &istate,
        const real_t evmaxmag_real,
        const real_t evmaxmag_imag,
        const bool delta_t_requested = false)
        const
        throw(std::exception)
    {
        SUZERAIN_UNUSED(delta_t_requested);
        real_t delta_t_candidates[2] = { numeric_limits<real_t>::max(),
                                         numeric_limits<real_t>::max()  };
        real_t &convective_delta_t = delta_t_candidates[0];
        real_t &diffusive_delta_t  = delta_t_candidates[1];
        const real_t one_over_delta_x = scenario.Lx / Nx; // !dNx, dealiasing
        const real_t one_over_delta_z = scenario.Lz / Nz; // !dNz, dealiasing

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
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rho.origin(),
            complex_zero, rho_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rho.origin(),
            complex_zero, rho_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of density at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rho.origin(),
            complex_zero, rho_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rho.origin(),
            complex_zero, rho_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);


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
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhou.origin(),
            complex_zero, mx_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhou.origin(),
            complex_zero, mx_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of X momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhou.origin(),
            complex_zero, mx_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhou.origin(),
            complex_zero, mx_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);


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
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhov.origin(),
            complex_zero, my_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, my_y.wave.begin(),
            complex_zero, my_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhov.origin(),
            complex_zero, my_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of Y momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhov.origin(),
            complex_zero, my_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, my_y.wave.begin(),
            complex_zero, my_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhov.origin(),
            complex_zero, my_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);


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
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhow.origin(),
            complex_zero, mz_xx.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 0, // dx dy
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_xy.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(1, 1, // dx dz
            complex_one, state_rhow.origin(),
            complex_zero, mz_xz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of Z momentum at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhow.origin(),
            complex_zero, mz_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 1, // dy dz
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_yz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhow.origin(),
            complex_zero, mz_zz.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);


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
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(2, 0, // d2x
            complex_one, state_rhoe.origin(),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Compute Z-related derivatives of total energy at collocation points
        suzerain::diffwave::accumulate(0, 1, // dz
            complex_one, state_rhoe.origin(),
            complex_zero, e_z.wave.begin(),
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);
        suzerain::diffwave::accumulate(0, 2, // d2z
            complex_one, state_rhoe.origin(),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            scenario.Lx, scenario.Lz,
            Ny, Nx, dNx, dkbx, dkex, Nz, dNz, dkbz, dkez);

        // Collectively convert state to physical space
        dgrid->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rho.origin()));
        dgrid->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhou.origin()));
        dgrid->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhov.origin()));
        dgrid->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhow.origin()));
        dgrid->transform_wave_to_physical(
                reinterpret_cast<real_t *>(state_rhoe.origin()));

        // Collectively convert state derivatives to physical space
        dgrid->transform_wave_to_physical(rho_x.data());   // density
        dgrid->transform_wave_to_physical(rho_y.data());
        dgrid->transform_wave_to_physical(rho_z.data());
        dgrid->transform_wave_to_physical(rho_xx.data());
        dgrid->transform_wave_to_physical(rho_xy.data());
        dgrid->transform_wave_to_physical(rho_xz.data());
        dgrid->transform_wave_to_physical(rho_yy.data());
        dgrid->transform_wave_to_physical(rho_yz.data());
        dgrid->transform_wave_to_physical(rho_zz.data());
        dgrid->transform_wave_to_physical(mx_x.data());   // X momentum
        dgrid->transform_wave_to_physical(mx_y.data());
        dgrid->transform_wave_to_physical(mx_z.data());
        dgrid->transform_wave_to_physical(mx_xx.data());
        dgrid->transform_wave_to_physical(mx_xy.data());
        dgrid->transform_wave_to_physical(mx_xz.data());
        dgrid->transform_wave_to_physical(mx_yy.data());
        dgrid->transform_wave_to_physical(mx_yz.data());
        dgrid->transform_wave_to_physical(mx_zz.data());
        dgrid->transform_wave_to_physical(my_x.data());   // Y momentum
        dgrid->transform_wave_to_physical(my_y.data());
        dgrid->transform_wave_to_physical(my_z.data());
        dgrid->transform_wave_to_physical(my_xx.data());
        dgrid->transform_wave_to_physical(my_xy.data());
        dgrid->transform_wave_to_physical(my_xz.data());
        dgrid->transform_wave_to_physical(my_yy.data());
        dgrid->transform_wave_to_physical(my_yz.data());
        dgrid->transform_wave_to_physical(my_zz.data());
        dgrid->transform_wave_to_physical(mz_x.data());   // Z momentum
        dgrid->transform_wave_to_physical(mz_y.data());
        dgrid->transform_wave_to_physical(mz_z.data());
        dgrid->transform_wave_to_physical(mz_xx.data());
        dgrid->transform_wave_to_physical(mz_xy.data());
        dgrid->transform_wave_to_physical(mz_xz.data());
        dgrid->transform_wave_to_physical(mz_yy.data());
        dgrid->transform_wave_to_physical(mz_yz.data());
        dgrid->transform_wave_to_physical(mz_zz.data());
        dgrid->transform_wave_to_physical(e_x.data());   // Total energy
        dgrid->transform_wave_to_physical(e_y.data());
        dgrid->transform_wave_to_physical(e_z.data());
        dgrid->transform_wave_to_physical(div_grad_e.data());

        // Compute nonlinear operator

        // Retrieve constants and compute derived constants
        const real_t beta             = scenario.beta;
        const real_t gamma            = scenario.gamma;
        const real_t Pr               = scenario.Pr;
        const real_t Re               = scenario.Re;
        const real_t inv_Re           = 1 / Re;
        const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));

        // Temporary storage used within following loop
        Eigen::Vector3r grad_rho;
        Eigen::Matrix3r grad_grad_rho;
        Eigen::Vector3r m;
        Eigen::Matrix3r grad_m;
        Eigen::Vector3r div_grad_m;
        Eigen::Vector3r grad_div_m;
        Eigen::Vector3r grad_e;
        Eigen::Vector3r u;
        Eigen::Matrix3r grad_u;
        Eigen::Vector3r grad_div_u, div_grad_u;
        Eigen::Vector3r grad_p, grad_T, grad_mu, grad_lambda;
        Eigen::Matrix3r tau;
        Eigen::Vector3r div_tau;

        // Used to track local \frac{1}{\Delta{}y} for time criterion.
        // In physical space storage is X Z Y with Y direction slowest
        // Maintain where we are relative to local rank's starting y offset
        const size_t stride_y = rho_x.physical.size_x * rho_x.physical.size_z;
        size_t ndx            = stride_y * rho_x.physical.start_y;

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
             ++p_div_grad_e,
             ++ndx) {

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
            grad_div_m[0]      = *p_mx_xx + *p_my_xy + *p_mz_xz;
            grad_div_m[1]      = *p_mx_xy + *p_my_yy + *p_mz_yz;
            grad_div_m[2]      = *p_mx_xz + *p_my_yz + *p_mz_zz;

            // Prepare local total energy-related quantities
            const real_t e          = *p_e;
            grad_e[0]               = *p_e_x;
            grad_e[1]               = *p_e_y;
            grad_e[2]               = *p_e_z;

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
                                        e, grad_e, *p_div_grad_e);
            const real_t div_grad_T = suzerain::orthonormal::rhome::div_grad_T(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        p, grad_p, div_grad_p);
            tau     = suzerain::orthonormal::tau(
                        mu, lambda, div_u, grad_u);
            div_tau = suzerain::orthonormal::div_tau(
                        mu, grad_mu, lambda, grad_lambda,
                        div_u, grad_u, div_grad_u, grad_div_u);

            // Continuity equation right hand side
            *p_rho = - div_m
                ;

            // Momentum equation right hand side
            Eigen::Vector3r momentum =
                - suzerain::orthonormal::div_u_outer_m(m, grad_m, u, div_u)
                - grad_p
                + inv_Re * div_tau
                ;
            *p_mx = momentum[0];
            *p_my = momentum[1];
            *p_mz = momentum[2];

            // Energy equation right hand side
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

            // Maintain the minimum observed stable time step
            convective_delta_t = suzerain::math::minnan(
                    convective_delta_t,
                    suzerain::timestepper::convective_stability_criterion(
                            u.x(), one_over_delta_x,
                            u.y(), one_over_delta_y[ndx / stride_y],
                            u.z(), one_over_delta_z,
                            evmaxmag_real,
                            std::sqrt(T)) /* nondimensional a = sqrt(T) */);
            diffusive_delta_t = suzerain::math::minnan(
                    diffusive_delta_t,
                    suzerain::timestepper::diffusive_stability_criterion(
                            one_over_delta_x,
                            one_over_delta_y[ndx / stride_y],
                            one_over_delta_z,
                            Re, Pr, gamma, evmaxmag_imag, mu / rho));
        }

        // Convert collocation point values to wave space
        dgrid->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rho.origin()));
        dgrid->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhou.origin()));
        dgrid->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhov.origin()));
        dgrid->transform_physical_to_wave(
                reinterpret_cast<real_t *>(state_rhow.origin()));
        dgrid->transform_physical_to_wave(
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

        // Perform Allreduce on stable time step sizes
        // Note delta_t_candidates aliases {convective,diffusive}_delta_t
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, delta_t_candidates,
                    sizeof(delta_t_candidates)/sizeof(delta_t_candidates[0]),
                    suzerain::mpi::datatype<real_t>::value,
                    MPI_MIN, MPI_COMM_WORLD));

        // Return minimum of either time step criterion, accounting for NaNs
        return suzerain::math::minnan(convective_delta_t, diffusive_delta_t);
    }

protected:

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

/** See writeup/channel_treatment.tex for full details */
class NonlinearOperatorWithBoundaryConditions : public NonlinearOperator
{
private:
    typedef NonlinearOperator base;

public:

    real_t applyOperator(
        istate_type &istate,
        const real_t evmaxmag_real,
        const real_t evmaxmag_imag,
        const bool delta_t_requested = false)
        const
        throw(std::exception)
    {
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

        // Special handling occurs only on rank holding the "zero-zero" mode
        const bool zero_zero_rank = (dkbx == 0) && (dkbz == 0);

        // Precompute and store quantities necessary for BC implementation
        Eigen::ArrayXr original_state_mx;
        real_t bulk_density = numeric_limits<real_t>::quiet_NaN();
        if (zero_zero_rank) {
            // Save a copy of the constant x-momentum modes
            original_state_mx.resize(state.shape()[1]);
            for (std::size_t i = 0; i < state.shape()[1]; ++i) {
                original_state_mx[i]
                    = suzerain::complex::real(state_rhou[i][0][0]);
            }

            // Compute the bulk density so we can hold it constant in time
            bspw->integrate(reinterpret_cast<real_t *>(state_rho.origin()),
                    sizeof(complex_t)/sizeof(real_t),
                    &bulk_density);
        }

        // Apply an operator that cares nothing about the boundaries
        const real_t delta_t = base::applyOperator(
                istate, evmaxmag_real, evmaxmag_imag, delta_t_requested);

        // Indices that will be useful as shorthand
        const std::size_t lower_wall = 0;                    // index of wall
        const std::size_t upper_wall = state.shape()[1] - 1; // index of wall

        // Add f_rho to mean density per writeup step (2)
        if (zero_zero_rank) {
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
        const real_t inv_gamma_gamma1
            = 1 / (scenario.gamma * (scenario.gamma - 1));
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                state_rhoe[lower_wall][j][k]
                    = inv_gamma_gamma1 * state_rho[lower_wall][j][k];
                state_rhoe[upper_wall][j][k]
                    = inv_gamma_gamma1 * state_rho[upper_wall][j][k];
            }
        }

        // Apply f_{m_x} to mean x-momentum, mean energy
        if (zero_zero_rank) {

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

        // Return the time step found by the BC-agnostic operator
        return delta_t;
    }
};

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

/** <tt>atexit</tt> callback to remove the metadata file. */
static void atexit_metadata(void) {
    if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        if (0 == unlink(restart.metadata().c_str())) {
            DEBUG("Cleaned up temporary file " << restart.metadata());
        } else {
            WARN("Error cleaning up temporary file " << restart.metadata());
        }
    }
}

/** Routine to load state from file.  */
static void load_state(esio_handle h, state_type &state)
{
    // Ensure local state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field_names.size());
    assert(numeric_cast<int>(state.shape()[1]) == dgrid->global_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid->local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid->global_wave_extent.z());

    // Obtain details on the restart field's global sizes
    // TODO Allow Ny != Fy per #1273
    int Fz, Fx, Fy, ncomponents;
    esio_field_sizev(h, field_names[0], &Fz, &Fx, &Fy, &ncomponents);
    assert(static_cast<int>(Fy) == dgrid->global_wave_extent.y());
    assert(ncomponents == 2);

    // Compute wavenumber translation logistics for X direction.
    // Requires turning a C2R FFT complex-valued coefficient count into a
    // real-valued coefficient count.  Further, need to preserve even- or
    // odd-ness of the coefficient count to handle, for example, Fx == 1.
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(2 * (Fx - 1) + (Fx & 1),
                                            grid.dN.x(),
                                            dgrid->local_wave_start.x(),
                                            dgrid->local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Y direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(Fz,
                                            grid.dN.z(),
                                            dgrid->local_wave_start.z(),
                                            dgrid->local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // FIXME Necessary when Fy != Ny
    // Allocate and clear temporary storage able to hold a single scalar field
    // boost::array<suzerain::pencil_grid::index,3> tmp_extent
    //     = {{ global_wave_extent.x(), Fy, global_wave_extent.z() }};
    // state_type tmp_state = suzerain::to_yxz(1, tmp_extent);
    // suzerain::multi_array::fill(tmp_state, 0);
    // Use this storage to transform between wall normal grids

    // Zero state_linear storage
    suzerain::multi_array::fill(*state_linear, 0);

    DEBUG0("Starting to load simulation fields");

    // Load each scalar field in turn...
    for (size_t i = 0; i < field_names.static_size; ++i) {

        // Create a view of the state for just the i-th scalar
        // Not strictly necessary, but aids greatly in debugging
        boost::multi_array_types::index_range all;
        boost::array_view_gen<state_type,3>::type field
            = (*state_linear)[boost::indices[i][all][all][all]];

        // ...which requires two reads per field, once per Z range
        for (int j = 0; j < 2; ++j) {

            // Destination of read is NULL for empty READ operations
            // Required since MultiArray triggers asserts on invalid indices
            complex_t * dest = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                dest = &field[0]
                             [mxb[0] - dgrid->local_wave_start.x()]
                             [mzb[j] - dgrid->local_wave_start.z()];
            }

            // Collectively establish size of read across all ranks
            esio_field_establish(h, Fz, fzb[j], (fze[j] - fzb[j]),
                                    Fx, fxb[0], (fxe[0] - fxb[0]),
                                    Fy,      0, (            Fy));

            // Perform collective read operation into state_linear
            complex_field_read(h, field_names[i], dest,
                               field.strides()[2],
                               field.strides()[1],
                               field.strides()[0]);
        }
    }
}

/** Routine to output status.  Signature for TimeController use. */
static bool log_status(real_t t, std::size_t nt) {
    INFO0("Simulation reached time " << t << " at time step " << nt);
    return true;
}

/** Tracks last time a restart file was written successfully */
static std::size_t last_restart_saved_nt = numeric_limits<std::size_t>::max();

/** Routine to store a restart file.  Signature for TimeController use. */
static bool save_restart(real_t t, std::size_t nt)
{
    esio_file_clone(esioh, restart.metadata().c_str(),
                    restart.uncommitted().c_str(), 1 /*overwrite*/);

    // Save simulation time information
    store_time(esioh, t);

    // Ensure local state storage meets this routine's assumptions
    assert(                  state_linear->shape()[0]  == field_names.size());
    assert(numeric_cast<int>(state_linear->shape()[1]) == dgrid->local_wave_extent.y());
    assert(numeric_cast<int>(state_linear->shape()[2]) == dgrid->local_wave_extent.x());
    assert(numeric_cast<int>(state_linear->shape()[3]) == dgrid->local_wave_extent.z());

    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(grid.N.x(),
                                            grid.dN.x(),
                                            dgrid->local_wave_start.x(),
                                            dgrid->local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(grid.N.z(),
                                            grid.dN.z(),
                                            dgrid->local_wave_start.z(),
                                            dgrid->local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    DEBUG0("Starting to store simulation fields at step " << nt);

    // Save each scalar field in turn...
    for (size_t i = 0; i < field_names.static_size; ++i) {

        // Create a view of the state for just the i-th scalar
        // Not strictly necessary, but aids greatly in debugging
        boost::multi_array_types::index_range all;
        boost::array_view_gen<state_type,3>::type field
            = (*state_linear)[boost::indices[i][all][all][all]];

        // ...which requires two writes per field, once per Z range
        for (int j = 0; j < 2; ++j) {

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &field[0]
                            [mxb[0] - dgrid->local_wave_start.x()]
                            [mzb[j] - dgrid->local_wave_start.z()];
            }

            // Collectively establish size of read across all ranks
            esio_field_establish(esioh,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Perform collective write operation from state_linear
            complex_field_write(esioh, field_names[i], src,
                                field.strides()[2],
                                field.strides()[1],
                                field.strides()[0],
                                field_descriptions[i]);
        }
    }

    esio_file_close_restart(esioh,
                            restart.desttemplate().c_str(),
                            restart.retain());

    INFO0("Successfully wrote restart at t = " << t << " for nt = " << nt);

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return true; // Continue time advancement
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

    DEBUG0("Processing command line arguments and response files");
    bool default_advance_dt;
    bool default_advance_nt;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based explicit compressible channel simulation");
        options.add_definition(
                const_cast<ScenarioDefinition<real_t>&>(scenario));
        options.add_definition(
                const_cast<GridDefinition&>(grid));
        options.add_definition(
                const_cast<RestartDefinition<>&>(restart));
        options.add_definition(
                const_cast<TimeDefinition<real_t>&>(timedef));
        options.process(argc, argv);

        default_advance_dt = options.variables()["advance_dt"].defaulted();
        default_advance_nt = options.variables()["advance_nt"].defaulted();

    }

    if (default_advance_dt && default_advance_nt) {
        FATAL0("Either --advance_dt or --advance_nt is required");
        return EXIT_FAILURE;
    }

    INFO0("Loading details from restart file: " << restart.load());
    esio_file_open(esioh, restart.load().c_str(), 0 /* read-only */);
    load(esioh, const_cast<ScenarioDefinition<real_t>&>(scenario));
    load(esioh, const_cast<GridDefinition&>(grid));
    load(esioh, bspw);
    esio_file_close(esioh);

    // TODO Account for B-spline differences at load time
    assert(bspw->order() == grid.k);

    INFO0("Saving metadata temporary file: " << restart.metadata());
    {
        esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_create(h, restart.metadata().c_str(), 1 /* overwrite */);
        store(h, scenario);
        store(h, grid, scenario.Lx, scenario.Lz);
        store(h, bspw);
        esio_file_close(h);
        esio_handle_finalize(h);
        atexit(&atexit_metadata); // Delete lingering metadata file at exit
    }

    // Initialize array holding \frac{1}{\Delta{}y} grid spacing
    one_over_delta_y.resize(grid.N.y());
    {
        // Determine minimum delta y observable from each collocation point
        real_t a, b, c;
        bspw->collocation_point(0, &a);                  // First point
        bspw->collocation_point(1, &b);
        one_over_delta_y[0] = std::abs(b - a);
        for (int i = 1; i < grid.N.y() - 1; ++i) {          // Intermediates
            bspw->collocation_point(i-1, &a);
            bspw->collocation_point(i,   &b);
            bspw->collocation_point(i+1, &c);
            one_over_delta_y[i] = std::min(std::abs(b-a), std::abs(c-b));
        }
        bspw->collocation_point(grid.N.y() - 2, &a);        // Last point
        bspw->collocation_point(grid.N.y() - 1, &b);
        one_over_delta_y[grid.N.y() - 1] = std::abs(b - a);

        // Invert to find \frac{1}{\Delta{}y}
        for (int i = 0; i < grid.N.y(); ++i) { // Invert
            one_over_delta_y[i] = 1 / one_over_delta_y[i];
        }
    }

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<suzerain::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    // Display global degree of freedom information
    INFO0("Global degrees of freedom  (DOF): " << grid.N.prod());
    INFO0("DOF by direction           (XYZ): " << grid.N);
    INFO0("Dealiased DOF by direction (XYZ): " << grid.dN);
    INFO0("Number of MPI ranks:              " << nranks);

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    dgrid = make_shared<suzerain::pencil_grid>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());
    INFO0("Rank grid used for decomposition: " << dgrid->processor_grid);
    DEBUG("Local wave start      (XYZ): " << dgrid->local_wave_start);
    DEBUG("Local wave end        (XYZ): " << dgrid->local_wave_end);
    DEBUG("Local wave extent     (XYZ): " << dgrid->local_wave_extent);
    DEBUG("Local physical start  (XYZ): " << dgrid->local_physical_start);
    DEBUG("Local physical end    (XYZ): " << dgrid->local_physical_end);
    DEBUG("Local physical extent (XYZ): " << dgrid->local_physical_extent);

    // Create state storage for linear operator
    // TODO Have state_linear only store non-dealiased state
    state_linear = make_shared<state_type>(
            suzerain::to_yxz(5, dgrid->local_wave_extent));

    // Load restart information into state_linear, including simulation time
    esio_file_open(esioh, restart.load().c_str(), 0 /* read-only */);
    real_t initial_t;
    load_time(esioh, initial_t);
    load_state(esioh, *state_linear); // May have large memory overhead!
    esio_file_close(esioh);

    // Create the state storage for nonlinear operator with appropriate padding
    // to allow P3DFFTification.  Must clear to avoid lingering NaN issues.
    {
        using suzerain::to_yxz;
        using suzerain::prepend;
        using suzerain::strides_cm;
        state_nonlinear = make_shared<state_type>(
                to_yxz(5, dgrid->local_wave_extent),
                prepend(dgrid->local_wave_storage(), strides_cm(
                        to_yxz(dgrid->local_wave_extent)))
                );
    }
    suzerain::multi_array::fill(*state_nonlinear, 0);

    // Dump some state shape and stride information for debugging purposes
    DEBUG("Linear state shape      (FYXZ): "
          << suzerain::multi_array::shape_array(*state_linear));
    DEBUG("Nonlinear state shape   (FYXZ): "
          << suzerain::multi_array::shape_array(*state_nonlinear));
    DEBUG("Linear state strides    (FYXZ): "
          << suzerain::multi_array::strides_array(*state_linear));
    DEBUG("Nonlinear state strides (FYXZ): "
          << suzerain::multi_array::strides_array(*state_nonlinear));

    // Instantiate the operators and time stepping details
    // See write up section 2.1 (Spatial Discretization) for coefficient origin
    using suzerain::timestepper::lowstorage::SMR91Method;
    const SMR91Method<complex_t> smr91(timedef.evmagfactor);
    MassOperator L(scenario.Lx * scenario.Lz * grid.N.x() * grid.N.z());
    NonlinearOperatorWithBoundaryConditions N;

    // Establish TimeController for use with operators and state storage
    using suzerain::timestepper::TimeController;
    scoped_ptr<TimeController<real_t> > tc(make_LowStorageTimeController(
                smr91, L, N, *state_linear, *state_nonlinear,
                initial_t, timedef.min_dt, timedef.max_dt));

    // Register status callbacks status_{dt,nt} as requested
    // When either is not provided, default to a reasonable behavior
    {
        TimeController<real_t>::time_type dt;
        if (timedef.status_dt) {
            dt = timedef.status_dt;
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain() / 5;
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (timedef.status_nt) {
            nt = timedef.status_nt;
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain() / 5;
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &log_status);
    }

    // Register restart-writing callbacks every_{dt,nt} as requested
    // When either is not provided, default to a reasonable behavior
    {
        TimeController<real_t>::time_type dt;
        if (restart.every_dt()) {
            dt = restart.every_dt();
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain();
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (restart.every_nt()) {
            nt = restart.every_nt();
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain();
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &save_restart);
    }

    // Advance time according to advance_dt, advance_nt criteria
    bool advance_success = true;
    switch ((!!timedef.advance_dt << 1) + !!timedef.advance_nt) {
        case 3:
            INFO0("Advancing simulation by at most " << timedef.advance_dt
                   << " units of physical time");
            INFO0("Advancing simulation by at most " << timedef.advance_nt
                   << " discrete time steps");
            advance_success = tc->advance(initial_t + timedef.advance_dt,
                                          timedef.advance_nt);
            break;
        case 2:
            INFO0("Advancing simulation by " << timedef.advance_dt
                  << " units of physical time");
            advance_success = tc->advance(initial_t + timedef.advance_dt);
            break;
        case 1:
            INFO0("Advancing simulation " << timedef.advance_nt
                   << " discrete time steps");
            advance_success = tc->step(timedef.advance_nt);
            break;
        case 0:
            if (!default_advance_dt) {
                INFO0("Advancing simulation until forcibly terminated");
                advance_success = tc->advance();
            } else if (!default_advance_nt) {
                WARN0("Simulation will not be advanced");
            } else {
                FATAL0("Sanity error in time control");
                return EXIT_FAILURE;
            }
            break;
        default:
            FATAL0("Sanity error in time control");
            return EXIT_FAILURE;
    }

    // Output statistics on time advancement
    if (!advance_success) {
        WARN0("TimeController stopped advancing time unexpectedly");
    }
    INFO0("Advanced simulation from t_initial = " << initial_t
          << " to t_final = " << tc->current_t()
          << " in " << tc->current_nt() << " steps");
    INFO0("Min/mean/max/stddev of delta_t: "
          << tc->taken_min()  << ", "
          << tc->taken_mean() << ", "
          << tc->taken_max()  << ", "
          << tc->taken_stddev());

    // Save a final restart before exit if one was not just saved
    if (advance_success && last_restart_saved_nt != tc->current_nt()) {
        INFO0("Saving final restart file prior to quitting");
        save_restart(tc->current_t(), tc->current_nt());
    }

    return advance_success ? EXIT_SUCCESS : EXIT_FAILURE;
}
