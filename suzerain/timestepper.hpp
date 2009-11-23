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
                     const suzerain::IState<FPT> * const input,
                           suzerain::IState<FPT> * const output) const
                     throw(suzerain::runtime_error) = 0;
    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> * const state) const
                     throw(suzerain::runtime_error) = 0;
}

template< typename FPT >
class ILowStorageMethod
{
public:
    virtual const char * name() const = 0;
    virtual std::size_t substeps() const = 0;
    virtual FPT alpha(std::size_t substep) const = 0;
    virtual FPT beta(std::size_t substep) const = 0;
    virtual FPT gamma(std::size_t substep) const = 0;
    virtual FPT delta(std::size_t substep) const = 0;
    virtual ~ILowStorageMethod() {};
}

template< typename FPT >
class SMR91Method : ILowStorageMethod<FPT>
{
public:
    virtual const char * name() { return "SMR91"; };
    virtual std::size_t substeps() { return 3; };
    virtual FPT alpha(const std::size_t substep);
    virtual FPT beta(const std::size_t substep);
    virtual FPT gamma(const std::size_t substep);
    virtual FPT zeta(const std::size_t substep);
    virtual ~ILowStorageMethod() {};
}

template< typename FPT >
FPT SMR91Method::alpha(const std::size_t substep)
{
    const FPT coeff[3] = { FPT(29)/FPT(96),  FPT(-3)/FPT(40), FPT(1)/FPT(6) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method::beta(const std::size_t substep)
{
    const FPT coeff[3] = { FPT(37)/FPT(160), FPT(5)/FPT(24), FPT(1)/FPT(6) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method::gamma(const std::size_t substep)
{
    const FPT coeff[3] = { FPT(8)/FPT(15), FPT(5)/FPT(12), FPT(3)/FPT(4) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method::zeta(const std::size_t substep)
{
    const FPT coeff[3] = { FPT(0), FPT(-17)/FPT(60), FPT(-5)/FPT(12) };
    return coeff[substep];
}

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
