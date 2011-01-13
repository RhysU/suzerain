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
#include <log4cxx/logger.h>
#include <esio/esio.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline_definition.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/state_impl.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/timestepper.hpp>
#include <suzerain/utility.hpp>

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
namespace sz = ::suzerain;

// Global scenario parameters and Bspline details initialized in main()
static sz::problem::ScenarioDefinition<>  def_scenario(100);
static sz::problem::ChannelDefinition<>   def_grid;
static sz::problem::BsplineDefinition<>   def_bspline;
static sz::problem::RestartDefinition<>   def_restart("",
                                                      "metadata.h5",
                                                      "uncommitted.h5",
                                                      "restart#.h5",
                                                      1);
static boost::shared_ptr<sz::bspline>     bspw;
static boost::shared_ptr<sz::bspline_luz> bspluzw;
static boost::shared_ptr<sz::pencil_grid> pg;

// Explicit timestepping scheme uses only complex double 4D NoninterleavedState
// State indices range over (scalar field, Y, X, Z) in wave space
// Establish some shorthand for some overly-templated operator and state types
typedef std::complex<double> complex_type;
typedef sz::storage::noninterleaved<4> storage_type;
typedef sz::IState<
            storage_type::dimensionality, complex_type, storage_type
        > istate_type;
typedef sz::timestepper::lowstorage::ILinearOperator<
            storage_type::dimensionality, complex_type, storage_type
        > ilinearoperator_type;
typedef sz::timestepper::INonlinearOperator<
            storage_type::dimensionality, complex_type, storage_type
        > inonlinearoperator_type;
typedef sz::NoninterleavedState<
            storage_type::dimensionality, complex_type
        > state_type;

// TODO Incorporate IOperatorLifecycle semantics
// TODO Refactor MassOperator into templated BsplineMassOperator

class MassOperator : public ilinearoperator_type
{
public:
    explicit MassOperator(double scaling = 1.0)
        : opscaling_(scaling), luzw_(*bspw)
    {
        complex_type coefficient;
        sz::complex::assign_complex(coefficient, scaling);
        luzw_.form_general(1, &coefficient, *bspw);
    }

    virtual void applyMassPlusScaledOperator(
            const complex_type scale,
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
            const complex_type scale,
            const istate_type &iinput,
            istate_type &ioutput) const throw (std::exception)
    {
        SUZERAIN_UNUSED(scale);
        const state_type &x = dynamic_cast<const state_type&>(iinput);
        state_type y = dynamic_cast<state_type&>(ioutput);
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
            const complex_type scale,
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
    const double opscaling_;
    sz::bspline_luz luzw_;
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

    double applyOperator(
        istate_type &istate,
        const bool delta_t_requested = false)
        const
        throw(std::exception) {

        SUZERAIN_UNUSED(delta_t_requested);
        double convective_delta_t = std::numeric_limits<double>::max();
        double diffusive_delta_t  = std::numeric_limits<double>::max();
        const double one_over_delta_x = def_grid.Lx() / def_grid.Nx();
        const double one_over_delta_y = def_grid.Ly() / def_grid.Ny();
        const double one_over_delta_z = def_grid.Lz() / def_grid.Nz();

        state_type &state = dynamic_cast<state_type&>(istate);

        const double complex_one[2]  = { 1.0, 0.0 };
        const double complex_zero[2] = { 0.0, 0.0 };

        // All state enters routine as coefficients in X, Y, and Z directions

        // Compute Y derivatives of density at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[0][0][0][0]), 1, state.shape()[1],
            complex_zero, rho_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[0][0][0][0]), 1, state.shape()[1],
            complex_zero, rho_yy.wave.begin(), 1, state.shape()[1]);
        // Compute density at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, &(state[0][0][0][0]), 1, state.shape()[1]);

        // Compute X-related derivatives of density at collocation points
        sz::diffwave::accumulate(1, 0, // dx
            complex_one, &(state[0][0][0][0]),
            complex_zero, rho_x.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(2, 0, // d2x
            complex_one, &(state[0][0][0][0]),
            complex_zero, rho_xx.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 0, // dx dy
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_xy.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 1, // dx dz
            complex_one, &(state[0][0][0][0]),
            complex_zero, rho_xz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Compute Z-related derivatives of density at collocation points
        sz::diffwave::accumulate(0, 1, // dz
            complex_one, &(state[0][0][0][0]),
            complex_zero, rho_z.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 1, // dy dz
            complex_one, rho_y.wave.begin(),
            complex_zero, rho_yz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 2, // d2z
            complex_one, &(state[0][0][0][0]),
            complex_zero, rho_zz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);


        // Compute Y derivatives of X momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[1][0][0][0]), 1, state.shape()[1],
            complex_zero, mx_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[1][0][0][0]), 1, state.shape()[1],
            complex_zero, mx_yy.wave.begin(), 1, state.shape()[1]);
        // Compute X momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, &(state[1][0][0][0]), 1, state.shape()[1]);

        // Compute X-related derivatives of X momentum at collocation points
        sz::diffwave::accumulate(1, 0, // dx
            complex_one, &(state[1][0][0][0]),
            complex_zero, mx_x.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(2, 0, // d2x
            complex_one, &(state[1][0][0][0]),
            complex_zero, mx_xx.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 0, // dx dy
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_xy.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 1, // dx dz
            complex_one, &(state[1][0][0][0]),
            complex_zero, mx_xz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Compute Z-related derivatives of X momentum at collocation points
        sz::diffwave::accumulate(0, 1, // dz
            complex_one, &(state[1][0][0][0]),
            complex_zero, mx_z.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 1, // dy dz
            complex_one, mx_y.wave.begin(),
            complex_zero, mx_yz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 2, // d2z
            complex_one, &(state[1][0][0][0]),
            complex_zero, mx_zz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);


        // Compute Y derivatives of Y momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[2][0][0][0]), 1, state.shape()[1],
            complex_zero, my_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[2][0][0][0]), 1, state.shape()[1],
            complex_zero, my_yy.wave.begin(), 1, state.shape()[1]);
        // Compute Y momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, &(state[2][0][0][0]), 1, state.shape()[1]);

        // Compute X-related derivatives of Y momentum at collocation points
        sz::diffwave::accumulate(1, 0, // dx
            complex_one, &(state[2][0][0][0]),
            complex_zero, my_x.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(2, 0, // d2x
            complex_one, &(state[2][0][0][0]),
            complex_zero, my_xx.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 0, // dx dy
            complex_one, my_y.wave.begin(),
            complex_zero, my_xy.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 1, // dx dz
            complex_one, &(state[2][0][0][0]),
            complex_zero, my_xz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Compute Z-related derivatives of Y momentum at collocation points
        sz::diffwave::accumulate(0, 1, // dz
            complex_one, &(state[2][0][0][0]),
            complex_zero, my_z.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 1, // dy dz
            complex_one, my_y.wave.begin(),
            complex_zero, my_yz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 2, // d2z
            complex_one, &(state[2][0][0][0]),
            complex_zero, my_zz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);


        // Compute Y derivatives of Z momentum at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[3][0][0][0]), 1, state.shape()[1],
            complex_zero, mz_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[3][0][0][0]), 1, state.shape()[1],
            complex_zero, mz_yy.wave.begin(), 1, state.shape()[1]);
        // Compute Y momentum at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, &(state[3][0][0][0]), 1, state.shape()[1]);

        // Compute X-related derivatives of Z momentum at collocation points
        sz::diffwave::accumulate(1, 0, // dx
            complex_one, &(state[3][0][0][0]),
            complex_zero, mz_x.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(2, 0, // d2x
            complex_one, &(state[3][0][0][0]),
            complex_zero, mz_xx.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 0, // dx dy
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_xy.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(1, 1, // dx dz
            complex_one, &(state[3][0][0][0]),
            complex_zero, mz_xz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Compute Z-related derivatives of Z momentum at collocation points
        sz::diffwave::accumulate(0, 1, // dz
            complex_one, &(state[3][0][0][0]),
            complex_zero, mz_z.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 1, // dy dz
            complex_one, mz_y.wave.begin(),
            complex_zero, mz_yz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 2, // d2z
            complex_one, &(state[3][0][0][0]),
            complex_zero, mz_zz.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);


        // Compute Y derivatives of total energy at collocation points
        bspw->accumulate_operator(1, // dy
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[4][0][0][0]), 1, state.shape()[1],
            complex_zero, e_y.wave.begin(), 1, state.shape()[1]);
        bspw->accumulate_operator(2, // d2y
            state.shape()[2]*state.shape()[3],
            complex_one, &(state[4][0][0][0]), 1, state.shape()[1],
            complex_zero, div_grad_e.wave.begin(), 1, state.shape()[1]);
        // Compute total energy at collocation points within state storage
        bspw->apply_operator(0, // eval
            state.shape()[2]*state.shape()[3],
            1.0, &(state[4][0][0][0]), 1, state.shape()[1]);

        // Compute X-related derivatives of total energy at collocation points
        sz::diffwave::accumulate(1, 0, // dx
            complex_one, &(state[4][0][0][0]),
            complex_zero, e_x.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(2, 0, // d2x
            complex_one, &(state[4][0][0][0]),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Compute Z-related derivatives of total energy at collocation points
        sz::diffwave::accumulate(0, 1, // dz
            complex_one, &(state[4][0][0][0]),
            complex_zero, e_z.wave.begin(),
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);
        sz::diffwave::accumulate(0, 2, // d2z
            complex_one, &(state[4][0][0][0]),
            complex_one, div_grad_e.wave.begin(), // sum with contents
            def_grid.Lx(), def_grid.Lz(),
            Ny, def_grid.Nx(), dNx, dkbx, dkex, def_grid.Nz(), dNz, dkbz, dkez);

        // Collectively convert state to physical space
        pg->transform_wave_to_physical(
                reinterpret_cast<double *>(&state[0][0][0][0]));
        pg->transform_wave_to_physical(
                reinterpret_cast<double *>(&state[1][0][0][0]));
        pg->transform_wave_to_physical(
                reinterpret_cast<double *>(&state[2][0][0][0]));
        pg->transform_wave_to_physical(
                reinterpret_cast<double *>(&state[3][0][0][0]));
        pg->transform_wave_to_physical(
                reinterpret_cast<double *>(&state[4][0][0][0]));

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
        const double beta             = def_scenario.beta();
        const double gamma            = def_scenario.gamma();
        const double Pr               = def_scenario.Pr();
        const double Re               = def_scenario.Re();
        const double inv_Re           = 1.0 / Re;
        const double inv_Re_Pr_gamma1 = 1.0 / (Re * Pr * (gamma - 1));

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
             double *p_rho    = reinterpret_cast<double *>(&state[0][0][0][0]),
                    *p_rho_x  = rho_x.physical.begin(),
                    *p_rho_y  = rho_y.physical.begin(),
                    *p_rho_z  = rho_z.physical.begin(),
                    *p_rho_xx = rho_xx.physical.begin(),
                    *p_rho_xy = rho_xy.physical.begin(),
                    *p_rho_xz = rho_xz.physical.begin(),
                    *p_rho_yy = rho_yy.physical.begin(),
                    *p_rho_yz = rho_yz.physical.begin(),
                    *p_rho_zz = rho_zz.physical.begin(),
                    *p_mx     = reinterpret_cast<double *>(&state[1][0][0][0]),
                    *p_mx_x   = mx_x.physical.begin(),
                    *p_mx_y   = mx_y.physical.begin(),
                    *p_mx_z   = mx_z.physical.begin(),
                    *p_mx_xx  = mx_xx.physical.begin(),
                    *p_mx_xy  = mx_xy.physical.begin(),
                    *p_mx_xz  = mx_xz.physical.begin(),
                    *p_mx_yy  = mx_yy.physical.begin(),
                    *p_mx_yz  = mx_yz.physical.begin(),
                    *p_mx_zz  = mx_zz.physical.begin(),
                    *p_my     = reinterpret_cast<double *>(&state[2][0][0][0]),
                    *p_my_x   = my_x.physical.begin(),
                    *p_my_y   = my_y.physical.begin(),
                    *p_my_z   = my_z.physical.begin(),
                    *p_my_xx  = my_xx.physical.begin(),
                    *p_my_xy  = my_xy.physical.begin(),
                    *p_my_xz  = my_xz.physical.begin(),
                    *p_my_yy  = my_yy.physical.begin(),
                    *p_my_yz  = my_yz.physical.begin(),
                    *p_my_zz  = my_zz.physical.begin(),
                    *p_mz     = reinterpret_cast<double *>(&state[3][0][0][0]),
                    *p_mz_x   = mz_x.physical.begin(),
                    *p_mz_y   = mz_y.physical.begin(),
                    *p_mz_z   = mz_z.physical.begin(),
                    *p_mz_xx  = mz_xx.physical.begin(),
                    *p_mz_xy  = mz_xy.physical.begin(),
                    *p_mz_xz  = mz_xz.physical.begin(),
                    *p_mz_yy  = mz_yy.physical.begin(),
                    *p_mz_yz  = mz_yz.physical.begin(),
                    *p_mz_zz  = mz_zz.physical.begin(),
                    *p_e      = reinterpret_cast<double *>(&state[4][0][0][0]),
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
            const double rho          = *p_rho;
            grad_rho[0]               = *p_rho_x;
            grad_rho[1]               = *p_rho_y;
            grad_rho[2]               = *p_rho_z;
            const double div_grad_rho = *p_rho_xx + *p_rho_yy + *p_rho_zz;
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
            const double div_m = *p_mx_x + *p_my_y + *p_my_z;
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
            const double e          = *p_e;
            grad_e[0]               = *p_e_x;
            grad_e[1]               = *p_e_y;
            grad_e[2]               = *p_e_z;
            const double div_grad_e = *p_div_grad_e;  // FIXME Shadow

            // Prepare quantities derived from local state and its derivatives
            u                  = sz::orthonormal::rhome::u(rho, m);
            const double div_u = sz::orthonormal::rhome::div_u(
                                    rho, grad_rho, m, div_m);
            grad_u             = sz::orthonormal::rhome::grad_u(
                                    rho, grad_rho, m, grad_m);
            grad_div_u         = sz::orthonormal::rhome::grad_div_u(
                                    rho, grad_rho, grad_grad_rho,
                                    m, div_m, grad_m, grad_div_m);
            div_grad_u         = sz::orthonormal::rhome::div_grad_u(
                                    rho, grad_rho, div_grad_rho,
                                    m, grad_m, div_grad_m);
            double p, T, mu, lambda;
            sz::orthonormal::rhome::p_T_mu_lambda(
                beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
                p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
            const double div_grad_p = sz::orthonormal::rhome::div_grad_p(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m,
                                        e, grad_e, div_grad_e);
            const double div_grad_T = sz::orthonormal::rhome::div_grad_T(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        p, grad_p, div_grad_p);
            tau     = sz::orthonormal::tau(mu, lambda, div_u, grad_u);
            div_tau = sz::orthonormal::div_tau(
                        mu, grad_mu, lambda, grad_lambda,
                        div_u, grad_u, div_grad_u, grad_div_u);

            // Maintain the minimum stable timestep
            // TODO Operator knows about the scheme's eigenvalues.  Fix that.
            convective_delta_t = std::min(convective_delta_t,
                sz::timestepper::convective_stability_criterion(
                    u.x(), one_over_delta_x,
                    u.y(), one_over_delta_y,
                    u.z(), one_over_delta_z,
                    std::sqrt(3.0),
                    std::sqrt(T))); // nondimensional a = sqrt(T)
            diffusive_delta_t = std::min(diffusive_delta_t,
                sz::timestepper::diffusive_stability_criterion(
                    one_over_delta_x,
                    one_over_delta_y,
                    one_over_delta_z,
                    Re, Pr, gamma, 2.512, mu / rho));

            // Continuity equation
            *p_rho = - div_m
                ;

            // Momentum equation
            Eigen::Vector3d momentum =
                - sz::orthonormal::div_u_outer_m(m, grad_m, u, div_u)
                - grad_p
                + inv_Re * div_tau
                ;
            *p_mx = momentum[0];
            *p_my = momentum[1];
            *p_mz = momentum[2];

            // Energy equation
            *p_e = - sz::orthonormal::div_e_u(e, grad_e, u, div_u)
                   - sz::orthonormal::div_p_u(p, grad_p, u, div_u)
                   + inv_Re_Pr_gamma1 * sz::orthonormal::div_mu_grad_T(
                        grad_T, div_grad_T, mu, grad_mu
                     )
                   + inv_Re * sz::orthonormal::div_tau_u<double>(
                        u, grad_u, tau, div_tau
                     )
                   ;
        }

        // TODO Enforce boundary conditions in wave space
        // Enforce isothermal lower wall boundary condition
        // Enforce isothermal upper wall boundary condition

        // Convert collocation point values to wave space
        pg->transform_physical_to_wave(
                reinterpret_cast<double *>(&state[0][0][0][0]));
        pg->transform_physical_to_wave(
                reinterpret_cast<double *>(&state[1][0][0][0]));
        pg->transform_physical_to_wave(
                reinterpret_cast<double *>(&state[2][0][0][0]));
        pg->transform_physical_to_wave(
                reinterpret_cast<double *>(&state[3][0][0][0]));
        pg->transform_physical_to_wave(
                reinterpret_cast<double *>(&state[4][0][0][0]));

        // Convert collocation point values to Bspline coefficients
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                &(state[0][0][0][0]), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                &(state[1][0][0][0]), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                &(state[2][0][0][0]), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                &(state[3][0][0][0]), 1, state.shape()[1]);
        bspluzw->solve(state.shape()[2]*state.shape()[3],
                &(state[4][0][0][0]), 1, state.shape()[1]);

        return std::min(convective_delta_t, diffusive_delta_t);
    }


private:
    // Details for sz::diffwave::* calls
    const int Ny;
    const int dNx;
    const int dkbx;
    const int dkex;
    const int dNz;
    const int dkbz;
    const int dkez;

    mutable sz::pencil<> rho_x, rho_y, rho_z;
    mutable sz::pencil<> rho_xx, rho_xy, rho_xz, rho_yy, rho_yz, rho_zz;

    mutable sz::pencil<> mx_x, mx_y, mx_z;
    mutable sz::pencil<> mx_xx, mx_xy, mx_xz, mx_yy, mx_yz, mx_zz;

    mutable sz::pencil<> my_x, my_y, my_z;
    mutable sz::pencil<> my_xx, my_xy, my_xz, my_yy, my_yz, my_zz;

    mutable sz::pencil<> mz_x, mz_y, mz_z;
    mutable sz::pencil<> mz_xx, mz_xy, mz_xz, mz_yy, mz_yz, mz_zz;

    mutable sz::pencil<> e_x, e_y, e_z, div_grad_e;
};

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Initialize logger using MPI environment details.
    const int nproc  = sz::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = sz::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr log = log4cxx::Logger::getLogger(
            sz::mpi::comm_rank_identifier(MPI_COMM_WORLD));
    // Log only warnings and above from ranks 1 and higher when not debugging
    if (procid > 0 && !log->isDebugEnabled()) {
        log->setLevel(log4cxx::Level::getWarn());
    }

    // Process command line arguments
    // TODO Rank 0 should read and broadcast these to all other ranks
    // Otherwise we will encounter an IO bottleneck on large job starts
    sz::ProgramOptions options(
            "Suzerain-based explicit compressible channel simulation");
    options.add_definition(def_scenario);
    options.add_definition(def_grid);
    options.add_definition(def_bspline);
    options.add_definition(def_restart);
    if (!procid) {
        options.process(argc, argv, MPI_COMM_WORLD);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv, MPI_COMM_WORLD,
                        nullstream, nullstream, nullstream, nullstream);
    }
    assert(def_grid.DAFy() == 1.0);  // Wall normal dealiasing disallowed

    // Dump relevant global scenario parameters
    LOG4CXX_INFO(log, "Number of MPI ranks:     " << nproc);
    LOG4CXX_INFO(log, "Global extents (XYZ):    " << def_grid.global_extents());
    LOG4CXX_INFO(log, "Dealiased extents (XYZ): " << def_grid.dealiased_extents());
    LOG4CXX_INFO(log, "B-spline order:          " << def_bspline.k());
    LOG4CXX_INFO(log, "Breakpoint stretching:   " << def_bspline.htdelta());
    LOG4CXX_INFO(log, "Reynolds number:         " << def_scenario.Re());
    LOG4CXX_INFO(log, "Prandtl number:          " << def_scenario.Pr());
    LOG4CXX_INFO(log, "gamma = C_p/C_v:         " << def_scenario.gamma());
    LOG4CXX_INFO(log, "beta = ln mu/ln T:       " << def_scenario.beta());

    // Initialize B-spline workspace using [0, Ly] with Ny degrees of freedom
    const int nbreakpoints = def_grid.Ny() - def_bspline.k() + 2;
    assert(nbreakpoints > 1);
    double *breakpoints
        = (double *)sz::blas::malloc(nbreakpoints*sizeof(double));
    assert(breakpoints);
    sz::math::linspace(0.0, 1.0, nbreakpoints, breakpoints); // Uniform [0, 1]
    for (int i = 0; i < nbreakpoints; ++i) {                 // Stretch 'em out
        breakpoints[i] = def_grid.Ly() * suzerain_htstretch2(
                def_bspline.htdelta(), 1.0, breakpoints[i]);
    }
    for (int i = 0; i < nbreakpoints; ++i) {
        LOG4CXX_TRACE(log,
                      "B-spline breakpoint[" << i << "] = " << breakpoints[i]);
    }
    bspw = boost::make_shared<sz::bspline>(
             def_bspline.k(), 2, nbreakpoints, breakpoints);
    assert(static_cast<unsigned>(bspw->ndof()) == def_grid.Ny());
    sz::blas::free(breakpoints);

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = boost::make_shared<sz::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    pg = boost::make_shared<sz::pencil_grid>(def_grid.dealiased_extents(),
                                             def_grid.processor_grid());
    LOG4CXX_INFO(log, "Processor grid used: " << pg->processor_grid());
    LOG4CXX_DEBUG(log, "Local dealiased wave start  (XYZ): "
                       << pg->local_wave_start());
    LOG4CXX_DEBUG(log, "Local dealiased wave end    (XYZ): "
                       << pg->local_wave_end());
    LOG4CXX_DEBUG(log, "Local dealiased wave extent (XYZ): "
                       << pg->local_wave_extent());

    // Compute how much non-dealiased XYZ state is local to this rank
    // Additional munging necessary X direction has (Nx/2+1) complex values
    const boost::array<sz::pencil_grid::index,3> state_start
        = pg->local_wave_start();
    const boost::array<sz::pencil_grid::index,3> state_end = {{
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[0]/2+1,
                                             pg->local_wave_end()[0]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[1],
                                             pg->local_wave_end()[1]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[2],
                                            pg->local_wave_end()[2])
    }};
    const boost::array<sz::pencil_grid::index,3> state_extent = {{
        std::max<sz::pencil_grid::index>(state_end[0] - state_start[0], 0),
        std::max<sz::pencil_grid::index>(state_end[1] - state_start[1], 0),
        std::max<sz::pencil_grid::index>(state_end[2] - state_start[2], 0)
    }};
    LOG4CXX_DEBUG(log, "Local state wave start  (XYZ): " << state_start);
    LOG4CXX_DEBUG(log, "Local state wave end    (XYZ): " << state_end);
    LOG4CXX_DEBUG(log, "Local state wave extent (XYZ): " << state_extent);

    // Create the state storage for the linear and nonlinear operators
    // with appropriate padding to allow nonlinear state to be P3DFFTified
    state_type state_linear(sz::to_yxz(5, state_extent));
    state_type state_nonlinear(
            sz::to_yxz(5, state_extent),
            sz::prepend(pg->local_wave_storage(),
                        sz::strides_cm(sz::to_yxz(pg->local_wave_extent()))));
    if (log->isDebugEnabled()) {
        boost::array<sz::pencil_grid::index,4> strides;
        std::copy(state_linear.strides(),
                  state_linear.strides() + 4, strides.begin());
        LOG4CXX_DEBUG(log, "Linear state strides    (FYXZ): " << strides);
        std::copy(state_nonlinear.strides(),
                  state_nonlinear.strides() + 4, strides.begin());
        LOG4CXX_DEBUG(log, "Nonlinear state strides (FYXZ): " << strides);
    }

    // Instantiate the operators and timestepping details
    // See write up section 2.1 (Spatial Discretization) for coefficient origin
    const sz::timestepper::lowstorage::SMR91Method<complex_type> smr91;
    MassOperator L(def_grid.Lx()*def_grid.Lz()*def_grid.Nx()*def_grid.Nz());
    NonlinearOperator N;

    // Take a timestep
    sz::timestepper::lowstorage::step(smr91, L, N, state_linear, state_nonlinear);
}
