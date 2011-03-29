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
// bspline_operators.hpp: low-storage operators build atop B-splines
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
#ifndef __SUZERAIN_BSPLINE_OPERATORS_HPP
#define __SUZERAIN_BSPLINE_OPERATORS_HPP

#include <suzerain/bspline.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>

/** @file Provides operators built atop B-splines */

namespace suzerain {

/** Forward declaration of BsplineMassOperator */
template<class State, class Enable = void> class BsplineMassOperator;

/**
 * BsplineMassOperator specialization for complex-valued 4D
 *  NoninterleavedState
 */
template<typename FPT>
class BsplineMassOperator<
    NoninterleavedState<4, std::complex<FPT> >,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type
> : public suzerain::timestepper::lowstorage::ILinearOperator<
        NoninterleavedState<4, std::complex<FPT> >
    >
{
public:
    typedef FPT real_t;
    typedef typename std::complex<real_t> complex_t;
    typedef NoninterleavedState<4,complex_t> state_type;

    explicit BsplineMassOperator(
            boost::shared_ptr<const suzerain::bspline> bspw,
            real_t scaling = 1)
        : bspw_(bspw), opscaling_(scaling), luzw_(*bspw)
    {
        complex_t coefficient;
        suzerain::complex::assign_complex(coefficient, scaling);
        luzw_.form_general(1, &coefficient, *bspw);
    }

    virtual void applyMassPlusScaledOperator(
            const complex_t &scale,
            state_type &state) const
    {
        SUZERAIN_UNUSED(scale);

        const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
        assert(1 == state.strides()[1]);
        assert(static_cast<unsigned>(luzw_.ndof()) == state.shape()[1]);
        bspw_->apply_operator(0, nrhs, opscaling_,
                state.memory_begin(), 1, state.strides()[2]);
    }

    virtual void accumulateMassPlusScaledOperator(
            const complex_t &scale,
            const state_type &input,
            state_type &output) const
    {
        SUZERAIN_UNUSED(scale);
        const state_type &x = input;  // Shorthand
        state_type &y       = output; // Shorthand
        assert(x.isIsomorphic(y));

        typedef typename state_type::index index;
        for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
            ix < static_cast<index>(x.index_bases()[0] + x.shape()[0]);
            ++ix, ++iy) {

            for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
                lx < static_cast<index>(x.index_bases()[3] + x.shape()[3]);
                ++lx, ++ly) {

                bspw_->accumulate_operator(0, x.shape()[2], opscaling_,
                        &x[ix][x.index_bases()[1]][x.index_bases()[2]][lx],
                        x.strides()[1], x.strides()[2],
                        1.0, &y[iy][y.index_bases()[1]][y.index_bases()[2]][ly],
                        y.strides()[1], y.strides()[2]);
            }
        }
    }

    virtual void invertMassPlusScaledOperator(
            const complex_t &scale,
            state_type &state) const
    {
        SUZERAIN_UNUSED(scale);

        const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
        assert(1 == state.strides()[1]);
        assert(static_cast<unsigned>(luzw_.ndof()) == state.shape()[1]);
        luzw_.solve(nrhs, state.memory_begin(), 1, state.strides()[2]);
    }

private:
    const boost::shared_ptr<const suzerain::bspline> bspw_;
    const real_t opscaling_;
    suzerain::bspline_luz luzw_;
};

} // end namespace suzerain

#endif /* __SUZERAIN_BSPLINE_OPERATORS_HPP */
