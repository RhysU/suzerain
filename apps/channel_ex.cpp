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
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bspline_definition.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
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
        assert(luzw_.ndof() == state.shape()[1]);
        assert(nrhs == state.strides()[0]);
        bspw->apply_operator(0, nrhs, opscaling_,
                state.memory_begin(), state.strides()[1], state.strides()[2]);
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
            ix < x.index_bases()[0] + x.shape()[0];
            ++ix, ++iy) {

            for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
                lx < x.index_bases()[3] + x.shape()[3];
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
        assert(luzw_.ndof() == state.shape()[1]);
        assert(nrhs == state.strides()[0]);
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

    complex_type applyOperator(
        istate_type &istate,
        const bool delta_t_requested = false)
        const
        throw(std::exception) {

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

        // TODO Compute nonlinear terms
        // TODO Compute stable timestep

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

        return complex_type(0);
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

// TODO This definition of breakpoint locations is lousy.
// Moser mentioned looking at Kolmogorov theory regarding the viscous scales
// near the wall, and then looking at their growth through the log layer.

/**
 * Compute \c n B-spline breakpoints locations across <tt>[xbegin,xend]</tt>
 * using stretching parameter \c alpha to cluster points near the edges.
 *
 * @param[in]  xbegin One edge
 * @param[in]  xend   The opposite edge
 * @param[in]  n      Number of breakpoints
 * @param[in]  alpha  Stretching parameter
 * @param[out] output Starting location in which to save the breakpoints:
 *                    <tt>output[0]</tt> to <tt>output[n-1]</tt> must
 *                    be valid storage locations.
 *
 * @see sz::math::stretchspace for more information on \c alpha
 */
template<typename FPT, typename SizeType>
static void compute_breakpoints(
        FPT xbegin, FPT xend, SizeType n, FPT alpha, FPT *output)
{
    using sz::math::stretchspace;
    const FPT xhalf = (xbegin + xend)/2;

    if (n & 1) {
        stretchspace(xbegin, xhalf, n/2+1, alpha,   output);
        stretchspace(xhalf,  xend,  n/2+1, 1/alpha, output+n/2);
    } else {
        const FPT aitch    = (xend-xbegin)/(n + 1);
        const FPT xhalflen = aitch*(n/2);
        stretchspace(xbegin, xbegin+xhalflen, n/2, alpha,   output);
        stretchspace(xend-xhalflen, xend,     n/2, 1/alpha, output+n/2);
    }

    assert(output[0] == xbegin);
    assert(output[n-1] == xend);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                   // Initialize MPI on startup
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit

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
    sz::ProgramOptions options;
    options.add_definition(def_scenario);
    options.add_definition(def_grid);
    options.add_definition(def_bspline);
    if (!procid) {
        options.process(argc, argv);
    } else {
        boost::onullstream nullstream;
        options.process(argc, argv,
                        nullstream, nullstream, nullstream, nullstream);
    }
    assert(def_grid.DAFy() == 1.0);  // Wall normal dealiasing disallowed

    // Dump relevant global scenario parameters
    LOG4CXX_INFO(log, "Number of MPI ranks:     " << nproc);
    LOG4CXX_INFO(log, "Global extents (XYZ):    " << def_grid.global_extents());
    LOG4CXX_INFO(log, "Dealiased extents (XYZ): " << def_grid.dealiased_extents());
    LOG4CXX_INFO(log, "B-spline order:          " << def_bspline.k());
    LOG4CXX_INFO(log, "B-spline stretching:     " << def_bspline.alpha());
    LOG4CXX_INFO(log, "Reynolds number:         " << def_scenario.Re());
    LOG4CXX_INFO(log, "Prandtl number:          " << def_scenario.Pr());
    LOG4CXX_INFO(log, "gamma = C_p/C_v:         " << def_scenario.gamma());
    LOG4CXX_INFO(log, "beta = ln mu/ln T:       " << def_scenario.beta());

    // Initialize B-spline workspace using [0, Ly] and Ny
    double *breakpoints
        = (double *)sz::blas::malloc(def_grid.Ny()*sizeof(double));
    assert(breakpoints);
    compute_breakpoints(0.0, def_grid.Ly(), def_grid.Ny(),
                        def_bspline.alpha(), breakpoints);
    for (int i = 0; i < def_grid.Ny(); ++i) {
        LOG4CXX_TRACE(log,
                      "B-spline breakpoint[" << i << "] = " << breakpoints[i]);
    }
    bspw = boost::make_shared<sz::bspline>(
             def_bspline.k(), 2, def_grid.Ny(), breakpoints);
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

    // Create the state storage for the linear and nonlinear operators
    // Compute how much non-dealiased XYZ state is local to this rank
    // Additional munging necessary X direction has (Nx/2+1) complex values
    const boost::array<sz::pencil_grid::index,3> state_start
        = pg->local_wave_start();
    const boost::array<sz::pencil_grid::index,3> state_end = {
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[0]/2+1,
                                             pg->local_wave_end()[0]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[1],
                                             pg->local_wave_end()[1]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[2],
                                            pg->local_wave_end()[2])
    };
    const boost::array<sz::pencil_grid::index,3> state_extent = {
        std::max<sz::pencil_grid::index>(state_end[0] - state_start[0], 0),
        std::max<sz::pencil_grid::index>(state_end[1] - state_start[1], 0),
        std::max<sz::pencil_grid::index>(state_end[2] - state_start[2], 0)
    };
    LOG4CXX_DEBUG(log, "Local state wave start  (XYZ): " << state_start);
    LOG4CXX_DEBUG(log, "Local state wave end    (XYZ): " << state_end);
    LOG4CXX_DEBUG(log, "Local state wave extent (XYZ): " << state_extent);
    // Create the state instances with appropriate padding
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

}
