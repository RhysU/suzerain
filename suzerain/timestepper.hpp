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

/** @file
 * Provides time integration schemes.
 */

namespace suzerain
{

/**
 * Provides time integration schemes and the associated interfaces that the
 * underlying operators must obey.  Also includes the associated operator
 * configuration and lifecycle interfaces.
 */
namespace timestepper
{

/**
 * A marker interface for configuration data supplied to operator at runtime.
 *
 * @see IOperatorLifecycle for more information.
 */
class IOperatorConfig
{
public:

    /** Virtual destructor to support interface-like behavior. */
    virtual ~IOperatorConfig() {};
};

/**
 * An interface that runtime configurable operators must implement.  The
 * interface follows the J2EE Servlet pattern.  Prior to being applied, the
 * init method should be called.  One or more operator application methods may
 * be called.  Finally, the destroy method will be called.  After destroy, init
 * may again be invoked with a new operator configuration.
 */
class IOperatorLifecycle
{
public:

    /**
     * Initialize the operator so that it may process operator application
     * requests.  The operator should use \c config to determine what resources
     * it needs, and it should use this information to obtain these resources.
     * All expensive or one-time set up cost should be completed during this
     * method.
     *
     * @param config operator configuration to set for subsequent
     *        operator application.
     * @throw std::exception if a non-recoverable error occurs
     *        during initialization.
     */
    virtual void init(const IOperatorConfig &config)
                      throw(std::exception) {};

    /**
     * Destroy any resources associated with the previous \c init invocation.
     */
    virtual void destroy() {};

    /** Virtual destructor to support interface-like behavior. */
    virtual ~IOperatorLifecycle() {};
};

/**
 * A marker interface for operator splitting data provided to the
 * operator at runtime.
 *
 * @see IAdjustableSplitOperator for more information.
 */
class IOperatorSplit
{
public:

    /** Virtual destructor to support interface-like behavior. */
    virtual ~IOperatorSplit() {};
};

/**
 * An interface indicating that an ILinearOperator and INonlinearOperator
 * implementation pair implement to indicate they support adjustable
 * operator splitting.
 */
class IAdjustableSplitOperator
{
public:

    /**
     * Establish an operator split using the given information.
     *
     * @param split operator split information to set for subsequent
     *        operator application.
     * @throw std::exception if a non-recoverable error occurs
     *        during operator splitting.
     */
    virtual void establishSplit(const IOperatorSplit &split)
                                throw(std::exception) {};

    /** Virtual destructor to support interface-like behavior. */
    virtual ~IAdjustableSplitOperator() {};
};

/**
 * Defines the nonlinear operator interface required for timestepping.
 *
 * @see suzerain::IState<FPT> for more information on state vectors.
 * @see suzerain::timestepper::lowstorage for low storage schemes.
 */
template< typename FPT >
class INonlinearOperator : public IOperatorLifecycle
{
public:

    /**
     * Apply the operator in place on a given state location.  That is, compute
     * \f$\mbox{state}\leftarrow{}\mbox{operator}\left(\mbox{state}\right)\f$.
     *
     * @param state The state location to use.  It is expected that certain
     *        implementations will require more specific state types.
     * @param delta_t_requested If true, the operator must compute
     *        and return a stable timestep.
     *
     * @return a stable timestep if \c delta_t_requested is true.  Otherwise
     *        the return value is meaningless.
     */
    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
                              throw(std::exception) = 0;
};


/**
 * Provides low-storage Runge-Kutta time integration schemes and the associated
 * operator interfaces.  A hybrid implicit/explicit scheme is available.  The
 * hybrid integrator advances the state vector \f$ u(t) \f$ to
 * \f$u(t+\Delta{}t)\f$ according to \f$ M u_{t} = Lu + N(u) \f$ where \f$M\f$,
 * \f$L\f$, and \f$N\f$ are linear, linear, and nonlinear operators,
 * respectively.  \f$M\f$ is referred to as the mass matrix.  No operator may
 * depend on time.
 *
 * @see ILowStorageMethod for details on this class of timestepping schemes.
 */
namespace lowstorage
{

/**
 * Defines the linear operator interface required for low storage timestepping.
 *
 * @see suzerain::IState<FPT> for more information on state vectors.
 */
template< typename FPT >
class ILinearOperator : public IOperatorLifecycle
{
public:

    /**
     * Apply \f$M+\phi{}L\f$ in-place for real scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)\mbox{state}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void applyMassPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception) = 0;

    /**
     * Accumulate \f$M+\phi{}L\f$ out-of-place for real scalar \f$\phi\f$.
     * That is, \f$\mbox{output}\leftarrow{}\mbox{output} +
     * \left(M+\phi{}L\right)\mbox{input}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param input State vector on which to apply the operator.
     * @param output State vector into which to accumulate the result.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void accumulateMassPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> &input,
                           suzerain::IState<FPT> &output) const
                     throw(std::exception) = 0;

    /**
     * Invert \f$M+\phi{}L\f$ in-place for real scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)^{-1}\mbox{state}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void invertMassPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception) = 0;
};

/**
 * Implements simple multiplicative operator which scales all state
 * variables by a uniform factor.  The associated mass matrix \f$M\f$
 * is identity.
 */
template< typename FPT >
class MultiplicativeOperator
    : public ILinearOperator<FPT>, public INonlinearOperator<FPT>
{
private:
    const FPT factor;   /**< Uniform scale factor to apply */
    const FPT delta_t;  /**< Uniform stable time step to report */

public:

    /**
     * Construct an instance which scales by \c factor and reports
     * \c delta_t as a stable timestep.
     *
     * @param factor uniform scaling factor to apply.
     * @param delta_t uniform stable timestep to return from ::applyOperator.
     *                If not provided, it defaults to NaN.
     */
    MultiplicativeOperator(
            const FPT factor,
            const FPT delta_t = std::numeric_limits<FPT>::quiet_NaN())
        : factor(factor), delta_t(delta_t) {};

    /**
     * Scale \c state by the factor set at construction time.
     *
     * @param state to scale in place.
     * @param delta_t_requested ignored in this implementation.
     *
     * @return The \c delta_t provided at construction time.
     */
    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
                              throw(std::exception)
    {
        state.scale(factor);
        return delta_t;
    };

    /**
     * Compute \f$\mbox{state}\leftarrow{}
     * \left(I+\phi\times\mbox{factor}\right) \mbox{state}\f$ where \c factor
     * is the scaling factor set at construction time.
     *
     * @param scale Additional scaling \f$\phi\f$ to apply.
     * @param state to modify in place.
     */
    virtual void applyMassPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception)
    {
        state.scale(1 + scale*factor);
    };

    /**
     * Compute \f$\mbox{output}\leftarrow{}\mbox{output}+
     * \left(I+\phi\times\mbox{factor}\right) \mbox{input}\f$ where \c factor
     * is the scaling factor set at construction time.
     *
     * @param scale Additional scaling \f$\phi\f$ to apply.
     * @param input on which to apply the operator.
     * @param output on which to accumulate the result.
     */
    virtual void accumulateMassPlusScaledOperator(
                     const FPT scale,
                     const suzerain::IState<FPT> &input,
                           suzerain::IState<FPT> &output) const
                     throw(std::exception)
    {
        output.addScaled(1 + scale*factor, input);
    };

    /**
     * Compute \f$\mbox{state}\leftarrow{}
     * \left(I+\phi\times\mbox{factor}\right)^{-1} \mbox{state}\f$ where \c
     * factor is the scaling factor set at construction time.
     *
     * @param scale Additional scaling \f$\phi\f$ to apply.
     * @param state to modify in place.
     */
    virtual void invertMassPlusScaledOperator(
                     const FPT scale,
                     suzerain::IState<FPT> &state) const
                     throw(std::exception)
    {
        state.scale(1/(1 + scale*factor));
    };
};

/**
 * Encapsulates a hybrid implicit/explicit low storage Runge-Kutta method to
 * advance \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$.
 * The method consists of one or more substeps governed by coefficients
 * \f$\alpha_i\f$, \f$\beta_i\f$, \f$\gamma_i\f$, and \f$\zeta_i\f$ where
 * \f$i\f$ is less than the number of substeps.  Each substep obeys
 * \f[
 *   \left(M - \Delta{}t\beta_{i}L\right) u^{i+1}
 *   =
 *   \left(M + \Delta{}t\alpha_{i}L\right) u^{i}
 *   + \Delta{}t\gamma_{i}N\left(u^{i}\right)
 *   + \Delta{}t\zeta_{i}N\left(u^{i-1}\right)
 * \f]
 * where \f$\alpha_i+\beta_i=\gamma_i+\zeta_i\f$.  Note that the indexing
 * on the \f$\zeta_i\f$ coefficients differs slightly from other sources.
 *
 * @see ILinearOperator for the interface that \f$L\f$ must implement.
 * @see INonlinearOperator for the interface that \f$N\f$ must implement.
 * @see SMR91Method for an example of a concrete scheme.
 * @see step() or substep() for methods that can advance state variables
 *      according to a timestepping method.
 */
template< typename FPT >
class ILowStorageMethod
{
public:

    /**
     * A human-readable name for the timestepping method.
     * Intended to be used in tracing and logging.
     *
     * @return The time advancement scheme's name.
     */
    virtual const char * name() const = 0;

    /**
     * The number of substeps required to advance \f$u(t)\f$ to
     * \f$u(t+\Delta{}t)\f$.
     *
     * @return The number of substeps per time step.
     */
    virtual std::size_t substeps() const = 0;

    /**
     * Obtain the scheme's \f$\alpha_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual FPT alpha(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\beta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual FPT beta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\gamma_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual FPT gamma(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\zeta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual FPT zeta(std::size_t substep) const = 0;

    /** Virtual destructor to support interface-like behavior. */
    virtual ~ILowStorageMethod() {};
};

/**
 * Output the timestepping scheme <tt>m</tt>'s name on the given
 * output stream.
 *
 * @param os output stream to use.
 * @param m scheme's name to output.
 *
 * @return The output stream.
 */
template< typename charT, typename traits, typename FPT >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const ILowStorageMethod<FPT> &m)
{
    return os << m.name();
}

/**
 * Encapsulates the three stage, third order scheme from Appendix A of Spalart,
 * Moser, and Rogers' 1991 ``Spectral Methods for the Navier-Stokes Equations
 * with One Infinite and Two Periodic Directions'' published in the
 * <em>Journal of Computational Physics</em> volume 96 pages 297-324.
 */
template< typename FPT >
class SMR91Method : public ILowStorageMethod<FPT>
{
public:

    /** Explicit default constructor */
    SMR91Method() {};

    /*! @copydoc ILowStorageMethod::name */
    virtual const char * name() const { return "SMR91"; }

    /*! @copydoc ILowStorageMethod::substeps */
    virtual std::size_t substeps() const { return 3; };

    /*! @copydoc ILowStorageMethod::alpha */
    virtual FPT alpha(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::beta */
    virtual FPT beta(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::gamma */
    virtual FPT gamma(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::zeta */
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

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme to advance the system \f$ M u_t = Lu +
 * N(u) \f$.  Note that the roles of state locations \c a and \c b are flipped
 * after each substep.
 *
 * @param m The low storage scheme to use.  For example, SMR91Method.
 * @param L The linear operator to be treated implicitly.
 * @param N The nonlinear operator to be treated explicitly.
 * @param delta_t The time step \f$\Delta{}t\f$ to take.  The same time step
 *                must be supplied for all substep computations.
 * @param a On entry contains \f$u^{i}\f$ and on exit contains
 *          \f$N\left(u^{i}\right)\f$.
 * @param b On entry contains \f$N\left(u^{i-1}\right)\f$ and on exit contains
 *          \f$u^{i+1}\f$.
 * @param substep_index The substep number to take.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
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
    L.accumulateMassPlusScaledOperator(
            delta_t * m.alpha(substep_index), a, b);
    N.applyOperator(a);
    b.addScaled(delta_t * m.gamma(substep_index), a);
    L.invertMassPlusScaledOperator( -delta_t * m.beta(substep_index), b);
}

/**
 * Using the given method and a linear and nonlinear operator, advance from
 * \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using a hybrid implicit/explicit scheme
 * to advance the system \f$ M u_t = Lu + N(u) \f$.
 *
 * @param m The low storage scheme to use.  For example, SMR91Method.
 * @param L The linear operator to be treated implicitly.
 * @param N The nonlinear operator to be treated explicitly.
 * @param delta_t The time step \f$\Delta{}t\f$ to take.
 * @param a On entry contains \f$u(t)\f$ and on exit contains
 *          \f$u(t+\Delta{}t)\f$.
 * @param b Used as a temporary storage location during the substeps.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
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

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu + N(u)\f$.
 * The time step taken, \f$\Delta{}t\f$ will be computed during the first
 * nonlinear operator application.
 *
 * @param m The low storage scheme to use.  For example, SMR91Method.
 * @param L The linear operator to be treated implicitly.
 * @param N The nonlinear operator to be treated explicitly.
 * @param a On entry contains \f$u^{i}\f$ and on exit contains
 *          \f$N\left(u^{i}\right)\f$.  The linear operator is applied
 *          only to this state storage.
 * @param b On entry contains \f$N\left(u^{i-1}\right)\f$ and on exit contains
 *          \f$u^{i+1}\f$.  The nonlinear operator is applied only
 *          to this state storage.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
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
    L.applyMassPlusScaledOperator(delta_t * m.alpha(0), a);
    a.addScaled(delta_t * m.gamma(0), b);
    L.invertMassPlusScaledOperator( -delta_t * m.beta(0), a);

    // Second and subsequent substeps are identical
    for (std::size_t i = 1; i < m.substeps(); ++i) {
        b.scale(delta_t * m.zeta(i));
        L.accumulateMassPlusScaledOperator(delta_t * m.alpha(i), a, b);
        b.exchange(a); // Note nonlinear storage controls exchange operation
        N.applyOperator(b, false /* delta_t not needed */);
        a.addScaled(delta_t * m.gamma(i), b);
        L.invertMassPlusScaledOperator( -delta_t * m.beta(i), a);
    }
}

} // namespace lowstorage

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
