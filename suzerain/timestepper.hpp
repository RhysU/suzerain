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
    virtual void init(const IOperatorConfig &config)
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
    virtual void establishSplit(const IOperatorSplit &split)
                                throw(std::exception) {};
    virtual ~IAdjustableSplitOperator() {};
};

template< typename FPT >
class INonlinearOperator : public IOperatorLifecycle
{
public:
    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
                              throw(std::exception) = 0;
};

namespace lowstorage
{

template< typename FPT >
class ILinearOperator : public IOperatorLifecycle
{
public:
    virtual void applyIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception) = 0;

    virtual void accumulateIdentityPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> &input,
                           suzerain::IState<FPT> &output) const
                     throw(std::exception) = 0;

    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception) = 0;
};

template< typename FPT >
class MultiplicativeOperator
    : public ILinearOperator<FPT>, public INonlinearOperator<FPT>
{
private:
    const FPT factor;
    const FPT delta_t;

public:
    MultiplicativeOperator(
            const FPT factor,
            const FPT delta_t = std::numeric_limits<FPT>::quiet_NaN())
        : factor(factor), delta_t(delta_t) {};

    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
                              throw(std::exception)
    {
        state.scale(factor);
        return delta_t;
    };

    virtual void applyIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception)
    {
        state.scale(1 + scale*factor);
    };

    virtual void accumulateIdentityPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> &input,
                           suzerain::IState<FPT> &output) const
                     throw(std::exception)
    {
        output.addScaled(1 + scale*factor, input);
    };

    virtual void invertIdentityPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception)
    {
        state.scale(1/(1 + scale*factor));
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

template< typename charT, typename traits, typename FPT >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const ILowStorageMethod<FPT> &m)
{
    return os << m.name();
}

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
void substep(const ILowStorageMethod<FPT> &m,
             const ILinearOperator<FPT> &L,
             const INonlinearOperator<FPT> &N,
             const FPT delta_t,
             IState<FPT> &a,
             IState<FPT> &b,
             const std::size_t substep_index)
throw(std::exception)
{
    if (substep_index >= m.substeps())
        throw std::invalid_argument("Requested substep too large");

    b.scale(delta_t * m.zeta(substep_index));
    L.accumulateIdentityPlusScaledOperator(
            delta_t * m.alpha(substep_index), a, b);
    N.applyOperator(a);
    b.addScaled(delta_t * m.gamma(substep_index), a);
    L.invertIdentityPlusScaledOperator( -delta_t * m.beta(substep_index), b);
}

template< typename FPT >
void step(const ILowStorageMethod<FPT> &m,
          const ILinearOperator<FPT> &L,
          const INonlinearOperator<FPT> &N,
          const FPT delta_t,
          IState<FPT> &a,
          IState<FPT> &b)
throw(std::exception)
{
    IState<FPT> *p_a = &a, *p_b = &b;

    // Even substep counts will the roles of a and b at return.
    // Perform one auxiliary flip for odd substep counts.
    // Possible since b = N(u_{i-1}) is wholly ignored for first substep
    if (m.substeps() & 1) {
        b = a; // TODO Use IState<FPT>::swap(IState<FPT>&) once available
        boost::swap(p_a, p_b);
    }

    for (std::size_t i = 0; i < m.substeps(); ++i) {
        substep(m, L, N, delta_t, *p_a, *p_b, i);
        boost::swap(p_a, p_b);
    }
}

template< typename FPT >
void substep(const ILowStorageMethod<FPT> &m,
             const INonlinearOperator<FPT> &N,
             const FPT delta_t,
             IState<FPT> &a,
             IState<FPT> &b,
             const std::size_t substep_index)
throw(std::exception)
{
    return substep<FPT>(
            m, MultiplicativeOperator<FPT>(0), N,
            delta_t, a, b, substep_index);
}

template< typename FPT >
void step(const ILowStorageMethod<FPT> &m,
          const INonlinearOperator<FPT> &N,
          const FPT delta_t,
          IState<FPT> &a,
          IState<FPT> &b)
throw(std::exception)
{
    return step(m, MultiplicativeOperator<FPT>(0), N, delta_t, a, b);
}

template< typename FPT, typename LinearState, typename NonlinearState >
void step(const ILowStorageMethod<FPT> &m,
          const ILinearOperator<FPT> &L,
          const INonlinearOperator<FPT> &N,
          LinearState &a,
          NonlinearState &b)
throw(std::exception)
{
    // First substep handling is special since we need to determine delta_t
    b = a;
    const FPT delta_t = N.applyOperator(b, true /* we need delta_t */);
    L.applyIdentityPlusScaledOperator(delta_t * m.alpha(0), a);
    a.addScaled(delta_t * m.gamma(0), b);
    L.invertIdentityPlusScaledOperator( -delta_t * m.beta(0), a);

    // Second and subsequent substeps are identical
    for (std::size_t i = 1; i < m.substeps(); ++i) {
        b.scale(delta_t * m.zeta(i));
        L.accumulateIdentityPlusScaledOperator(delta_t * m.alpha(i), a, b);
        b.exchange(a); // Note nonlinear storage controls exchange operation
        N.applyOperator(b, false /* delta_t not needed */);
        a.addScaled(delta_t * m.gamma(i), b);
        L.invertIdentityPlusScaledOperator( -delta_t * m.beta(i), a);
    }
}

template< typename FPT, typename LinearState, typename NonlinearState >
void step(const ILowStorageMethod<FPT> &m,
          const INonlinearOperator<FPT> &N,
          LinearState &a,
          NonlinearState &b)
throw(std::exception)
{
    return step(m, MultiplicativeOperator<FPT>(0), N, a, b);
}

} // namespace lowstorage

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
