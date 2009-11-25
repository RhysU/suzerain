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

class IOperatorLifecycle
{
public:
    virtual void init(const IOperatorConfig * const config)
                      throw(std::runtime_error) {};
    virtual void destroy() {};
    virtual ~IOperatorLifecycle() {};
};

class IOperatorSplit
{
public:
    virtual ~IOperatorSplit() {};
};

class IAdjustableSplitOperator
{
public:
    virtual void establishSplit(const IOperatorSplit * const split)
                                throw(std::exception) {};
    virtual ~IAdjustableSplitOperator() {};
};

template< typename FPT >
class INonlinearOperator : public IOperatorLifecycle
{
public:
    virtual void applyOperator(suzerain::IState<FPT> * const state) const
                               throw(std::exception) = 0;
};

namespace lowstorage
{

template< typename FPT >
class ILinearOperator : public IOperatorLifecycle
{
public:
    virtual void accumulateIdentityPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> * const input,
                           suzerain::IState<FPT> * const output) const
                     throw(std::exception) = 0;
    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> * const state) const
                     throw(std::exception) = 0;
};

template< typename FPT >
class MultiplicativeOperator
    : public ILinearOperator<FPT>, public INonlinearOperator<FPT>
{
private:
    const FPT factor;

public:
    MultiplicativeOperator(const FPT factor) : factor(factor) {};

    virtual void applyOperator(suzerain::IState<FPT> * const state) const
                               throw(std::exception)
    {
        state->scale(factor);
    };

    virtual void accumulateIdentityPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> * const input,
                           suzerain::IState<FPT> * const output) const
                     throw(std::exception)
    {
        output->addScaled(1.0 + scale*factor, input);
    };

    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> * const state) const
                     throw(std::exception)
    {
        state->scale(1/(1+scale*factor));
    };
};

template< typename FPT >
class ILowStorageMethod
{
public:
    virtual const char * name() const = 0;
    virtual std::size_t substeps() const = 0;
    virtual FPT alpha(std::size_t substep) const = 0;
    virtual FPT beta(std::size_t substep) const = 0;
    virtual FPT gamma(std::size_t substep) const = 0;
    virtual FPT zeta(std::size_t substep) const = 0;
    virtual ~ILowStorageMethod() {};
};

template< typename FPT >
class SMR91Method : public ILowStorageMethod<FPT>
{
public:
    SMR91Method() {};
    virtual const char * name() const { return "SMR91"; }
    virtual std::size_t substeps() const { return 3; };
    virtual FPT alpha(const std::size_t substep) const;
    virtual FPT beta(const std::size_t substep) const;
    virtual FPT gamma(const std::size_t substep) const;
    virtual FPT zeta(const std::size_t substep) const;
};

template< typename FPT >
FPT SMR91Method<FPT>::alpha(const std::size_t substep) const
{
    const FPT coeff[3] = { FPT(29)/FPT(96),  FPT(-3)/FPT(40), FPT(1)/FPT(6) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method<FPT>::beta(const std::size_t substep) const
{
    const FPT coeff[3] = { FPT(37)/FPT(160), FPT(5)/FPT(24), FPT(1)/FPT(6) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method<FPT>::gamma(const std::size_t substep) const
{
    const FPT coeff[3] = { FPT(8)/FPT(15), FPT(5)/FPT(12), FPT(3)/FPT(4) };
    return coeff[substep];
}

template< typename FPT >
FPT SMR91Method<FPT>::zeta(const std::size_t substep) const
{
    const FPT coeff[3] = { FPT(0), FPT(-17)/FPT(60), FPT(-5)/FPT(12) };
    return coeff[substep];
}

template< typename FPT >
void substep(const ILowStorageMethod<FPT> * const m,
             const INonlinearOperator<FPT> * const N,
             const ILinearOperator<FPT> * const L,
             const FPT delta_t,
             IState<FPT> * const a,
             IState<FPT> * const b,
             const std::size_t substep_index)
throw(std::exception)
{
    if (substep_index >= m->substeps())
        throw std::invalid_argument("Requested substep too large");

    b->scale(delta_t * m->zeta(substep_index));
    L->accumulateIdentityPlusScaledOperator(
            delta_t * m->alpha(substep_index), a, b);
    N->applyOperator(a);
    b->addScaled(delta_t * m->gamma(substep_index), a);
    L->invertIdentityPlusScaledOperator(
            - delta_t * m->beta(substep_index), b);
}

template< typename FPT >
void substep(const ILowStorageMethod<FPT> * const m,
             const INonlinearOperator<FPT> * const N,
             const FPT delta_t,
             IState<FPT> * const a,
             IState<FPT> * const b,
             const std::size_t substep_index)
throw(std::exception)
{
    MultiplicativeOperator<FPT> zero_operator(FPT(0));
    return substep<FPT>(m, N, &zero_operator, delta_t, a, b, substep_index);
}

} // namespace lowstorage

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
