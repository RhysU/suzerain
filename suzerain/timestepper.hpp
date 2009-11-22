/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * timestepper.hpp: low storage Runge-Kutta timestepper interface
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_TIMESTEPPER_HPP
#define __SUZERAIN_TIMESTEPPER_HPP

#include <suzerain/common.hpp>
#include <suzerain/exceptions.hpp>
#include <suzerain/state.hpp>

namespace suzerain
{

namespace timestepper
{

class IOperatorConfig
{
public:
    virtual ~IOperatorConfig() {};
};

class IOperatorSplit
{
public:
    virtual ~IOperatorSplit() {};
};

class IOperatorLifecycle
{
public:
    virtual void init(const IOperatorConfig * const config)
                      throw(suzerain::runtime_error) {};
    virtual void establishSplit(const IOperatorSplit * const split)
                                throw(suzerain::runtime_error) {};
    virtual void destroy() {};
    virtual ~IOperatorLifecycle() {};
};

template< typename FPT >
class INonlinearOperator : public IOperatorLifecycle
{
public:
    virtual void applyOperator(suzerain::IState<FPT> * const state) const
                               throw(suzerain::runtime_error) = 0;
}

template< typename FPT >
class ILinearOperator : public IOperatorLifecycle
{
public:
    virtual void applyIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> * const state) const
                     throw(suzerain::runtime_error) = 0;
    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> * const state) const
                     throw(suzerain::runtime_error) = 0;
}

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
