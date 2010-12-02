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
#include <suzerain/traits.hpp>
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
    virtual void init(const IOperatorConfig& config)
                      throw(std::exception) {
        SUZERAIN_UNUSED(config);
    };

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
    virtual void establishSplit(const IOperatorSplit& split)
                                throw(std::exception) {
        SUZERAIN_UNUSED(split);
    };

    /** Virtual destructor to support interface-like behavior. */
    virtual ~IAdjustableSplitOperator() {};
};

/**
 * Defines the nonlinear operator interface required for timestepping.
 *
 * @see suzerain::IState<NumDims,Element,Storage,CompatibleStorage>
 *      for more information on state vectors.
 * @see suzerain::timestepper::lowstorage for low storage schemes.
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Storage,
    typename CompatibleStorage = Storage
>
class INonlinearOperator : public virtual IOperatorLifecycle
{
public:

    /**
     * Apply the operator in place on a given state location.  That is, compute
     * \f$\mbox{state}\leftarrow{}\mbox{operator}\left(\mbox{state}\right)\f$.
     *
     * @param state The state location to use.  It is expected that certain
     *        implementations will require more specific state types.
     * @param delta_t_requested If true, the operator must compute
     *        and return a stable time step.
     *
     * @return a stable time step if \c delta_t_requested is true.  Otherwise
     *        the return value is meaningless.
     */
    virtual typename suzerain::traits::component<Element>::type applyOperator(
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state,
        const bool delta_t_requested = false)
        const
        throw(std::exception) = 0;
};

/**
 * Compute an approximation to the maximum stable time step according to the
 * convective criterion for the three dimensional Euler equations.  This
 * surrogate is a model for the convective portion of the Navier-Stokes
 * operator.  This criterion appears as equation 2.39 in Wai Y. Kwok's thesis
 * (2002) and equations 4.20 and 4.21 in Stephen Guarini's thesis (1998):
 * \f[
 *     \pi\left(
 *         \frac{\left|u_{x}\right| + a}{\Delta{}x}
 *       + \frac{\left|u_{y}\right| + a}{\Delta{}y}
 *       + \frac{\left|u_{z}\right| + a}{\Delta{}z}
 *     \right) \Delta{}t \leq \left|\lambda_{I}\Delta{}t\right|_{\mbox{max}}
 * \f]
 * where \f$a\f$ is the local acoustic velocity, \f$u_{x}\f$ denotes the
 * velocity in the X direction, \f$\Delta{}x\f$ represents the grid size in the
 * X direction, etc.  The maximum pure imaginary eigenvalue magnitude,
 * \f$\left|\lambda_{I}\Delta{}t\right|_{\mbox{max}}\f$, is a feature of the
 * chosen timestepping method.  For example, it is \f$\sqrt{3}\f$ for the SMR91
 * scheme.
 *
 * @note Using a hybrid implicit/explicit %timestepper with acoustic terms
 * computed implicitly effectively sets the sound speed to be zero for this CFL
 * calculation.
 *
 * @param u_x              Velocity in the X direction \f$u_{x}\f$
 * @param one_over_delta_x Inverse local X grid spacing \f$1/\Delta{}x\f$
 * @param u_y              Velocity in the Y direction \f$u_{y}\f$
 * @param one_over_delta_y Inverse local Y grid spacing \f$1/\Delta{}y\f$
 * @param u_z              Velocity in the Z direction \f$u_{z}\f$
 * @param one_over_delta_z Inverse local Z grid spacing \f$1/\Delta{}z\f$
 * @param evmaxmag_imag    The maximum pure imaginary eigenvalue magnitude
 *        for some Runge-Kutta scheme, denoted
 *        \f$\left|\lambda_{I}\Delta_{}t\right|_{\mbox{max}}\f$
 *        in Guarini's thesis.
 * @param a                The local sound speed \f$a\f$
 *
 * @return The maximum stable time step \f$\Delta{}t\f$ according to
 *         the convective CFL criterion.
 */
template<typename FPT>
FPT convective_stability_criterion(
        const FPT u_x,
        const FPT one_over_delta_x,
        const FPT u_y,
        const FPT one_over_delta_y,
        const FPT u_z,
        const FPT one_over_delta_z,
        const FPT evmaxmag_imag,
        const FPT a = 0)
{
    // Precision for a 128-bit IEEE quad found via Sage's N(1/pi,digits=34)
    const FPT one_over_pi = (FPT) 0.3183098861837906715377675267450287L;
    return (evmaxmag_imag * one_over_pi)
        /  (   (std::abs(u_x) + a)*one_over_delta_x
             + (std::abs(u_y) + a)*one_over_delta_y
             + (std::abs(u_z) + a)*one_over_delta_z);
}

/**
 * Compute an approximation to the maximum stable time step according to the
 * viscous stability criterion for a three dimensional model diffusion
 * equation.  This surrogate is a model for the diffusive part of the
 * Navier-Stokes operator.  This criterion appears as equation 2.40 in Wai Y.
 * Kwok's thesis (2002) and equations 4.29 and 4.30 in Stephen Guarini's thesis
 * (1998):
 * \f[
 *   \mbox{max}\!\left(
 *     \left|\frac{\gamma\left(\nu-\nu_{0}\right)}{\mbox{Re}\mbox{Pr}}\right|,
 *     \left|\frac{\nu-\nu_{0}}{\mbox{Re}}\right|
 *   \right)
 *   \pi^{2}
 *   \left(
 *       \frac{1}{\Delta{}x^{2}}
 *     + \frac{1}{\Delta{}y^{2}}
 *     + \frac{1}{\Delta{}z^{2}}
 *   \right)
 *   \Delta{}t \leq \left|\lambda_{R}\Delta_{}t\right|_{\mbox{max}}
 * \f]
 * The maximum pure real eigenvalue magnitude,
 * \f$\left|\lambda_{R}\Delta{}t\right|_{\mbox{max}}\f$, is a feature of the
 * chosen timestepping method.  For example, it is 2.512 for the SMR91 scheme.
 * The absolute values within the maximum operation account for the possibility
 * that \f$\nu<\nu_{0}\f$.
 *
 * @note Using a hybrid implicit/explicit %timestepper with viscous terms
 * computed implicitly sets \f$\nu_0\f$ to be the reference kinematic viscosity
 * about which the viscous terms were linearized.
 *
 * @param one_over_delta_x Inverse local X grid spacing \f$1/\Delta{}x\f$
 * @param one_over_delta_y Inverse local Y grid spacing \f$1/\Delta{}y\f$
 * @param one_over_delta_z Inverse local Z grid spacing \f$1/\Delta{}z\f$
 * @param Re               The Reynolds number \f$\mbox{Re}\f$
 * @param Pr               The Prandtl number \f$\mbox{Pr}\f$
 * @param gamma            The ratio of specific heats \f$\gamma\f$
 * @param evmaxmag_real    The maximum pure real eigenvalue magnitude
 *        for some Runge-Kutta scheme, denoted
 *        \f$\left|\lambda_{R}\Delta_{}t\right|_{\mbox{max}}\f$
 *        in Guarini's thesis.
 * @param nu                The local kinematic viscosity \f$\nu\f$
 * @param nu0               The kinematic viscosity reference value \f$\nu_{0}\f$
 *
 * @return The maximum stable time step \f$\Delta{}t\f$ according to
 *         the diffusive criterion.
 */
template<typename FPT>
FPT diffusive_stability_criterion(
        const FPT one_over_delta_x,
        const FPT one_over_delta_y,
        const FPT one_over_delta_z,
        const FPT Re,
        const FPT Pr,
        const FPT gamma,
        const FPT evmaxmag_real,
        const FPT nu,
        const FPT nu0 = 0)
{
    // Precision for a 128-bit quad found via Sage's N(1/(pi*pi),digits=34)
    const FPT one_over_pi_squared = (FPT)0.1013211836423377714438794632097276L;
    const FPT nu_less_nu0 = std::abs(nu - nu0);
    const FPT maxcoeff = std::max((gamma*nu_less_nu0)/(Re*Pr), nu_less_nu0/Re);
    return (evmaxmag_real * one_over_pi_squared)
        /  (maxcoeff * (   one_over_delta_x*one_over_delta_x
                         + one_over_delta_y*one_over_delta_y
                         + one_over_delta_z*one_over_delta_z));
}

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
 * @see suzerain::IState<NumDims,Element,Storage,CompatibleStorage>
 *      for more information on state vectors.
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Storage,
    typename CompatibleStorage = Storage
>
class ILinearOperator : public virtual IOperatorLifecycle
{
public:

    /**
     * Apply \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)\mbox{state}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void applyMassPlusScaledOperator(
        const Element scale,
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state)
        const
        throw(std::exception) = 0;

    /**
     * Accumulate \f$M+\phi{}L\f$ out-of-place for scalar \f$\phi\f$.
     * That is, \f$\mbox{output}\leftarrow{}\mbox{output} +
     * \left(M+\phi{}L\right)\mbox{input}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param input State vector on which to apply the operator.
     * @param output State vector into which to accumulate the result.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void accumulateMassPlusScaledOperator(
        const Element scale,
        const suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& input,
        suzerain::IState<NumDims,Element,CompatibleStorage,Storage>& output)
        const
        throw(std::exception) = 0;

    /**
     * Invert \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)^{-1}\mbox{state}\f$.
     *
     * @param scale Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @throw std::exception if an unrecoverable error occurs.
     */
    virtual void invertMassPlusScaledOperator(
        const Element scale,
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state)
        const
        throw(std::exception) = 0;
};

/**
 * Implements simple multiplicative operator which scales all state
 * variables by a uniform factor.  The associated mass matrix \f$M\f$
 * is the identity matrix.
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Storage,
    typename CompatibleStorage = Storage
>
class MultiplicativeOperator
    : public ILinearOperator<NumDims,Element,Storage,CompatibleStorage>,
      public INonlinearOperator<NumDims,Element,Storage,CompatibleStorage>
{
private:
    /** Uniform scale factor to apply as part of operator */
    const Element factor;

    /** Uniform stable time step to report */
    const typename suzerain::traits::component<Element>::type delta_t;

public:

    /**
     * Construct an instance which scales by \c factor and reports
     * \c delta_t as a stable time step.
     *
     * @param factor uniform scaling factor to apply.
     * @param delta_t uniform, presumably stable time step to
     *        always return from ::applyOperator.
     */
    template< typename FactorType, typename DeltaTType >
    MultiplicativeOperator(
            const FactorType& factor,
            const DeltaTType& delta_t)
        : factor(factor), delta_t(delta_t) {};

    /**
     * Construct an instance which scales by \c factor and reports the
     * maximum representable floating point value as a stable time step.
     * Useful in testing contexts or during operator composition.
     *
     * @param factor uniform scaling factor to apply.
     */
    template< typename FactorType >
    MultiplicativeOperator(
            const FactorType& factor)
        : factor(factor),
          delta_t(std::numeric_limits<
                    typename suzerain::traits::component<Element>::type
                >::quiet_NaN()) {};

    /**
     * Scale \c state by the factor set at construction time.
     *
     * @param state to scale in place.
     * @param delta_t_requested ignored in this implementation.
     *
     * @return The \c delta_t provided at construction time.
     */
    virtual typename suzerain::traits::component<Element>::type applyOperator(
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state,
        const bool delta_t_requested = false)
        const
        throw(std::exception)
    {
        SUZERAIN_UNUSED(delta_t_requested);
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
        const Element scale,
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state)
        const
        throw(std::exception)
    {
        // Assumes Element has operator* and operator+
        state.scale(scale*factor + Element(1));
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
        const Element scale,
        const suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& input,
        suzerain::IState<NumDims,Element,CompatibleStorage,Storage>& output)
        const
        throw(std::exception)
    {
        // Assumes Element has operator* and operator+
        output.addScaled(scale*factor + Element(1), input);
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
        const Element scale,
        suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& state)
        const
        throw(std::exception)
    {
        // Assumes Element has a single argument constructor
        state.scale((Element(1))/(scale*factor + Element(1)));
    };
};

/**
 * Create a multiplicative operator matching particular state types.
 *
 * @see MultiplicativeOperator<NumDims,Element,Storage,CompatibleStorage> for
 *      more details.
 */
template<
    typename FactorType,
    typename DeltaTType,
    std::size_t NumDims,
    typename Element,
    typename Storage,
    typename CompatibleStorage
>
MultiplicativeOperator<NumDims,Element,Storage,CompatibleStorage>
make_multiplicator_operator(
    const FactorType& factor,
    const DeltaTType& delta_t,
    const suzerain::IState<NumDims,Element,Storage,CompatibleStorage>& input,
    const suzerain::IState<NumDims,Element,CompatibleStorage,Storage>& output)
{
    MultiplicativeOperator<NumDims,Element,Storage,CompatibleStorage>
        retval(factor,delta_t);
    return retval;
}

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
template< typename Element >
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
    virtual Element alpha(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\beta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual Element beta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\gamma_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual Element gamma(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\zeta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual Element zeta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's maximum pure real eigenvalue magnitude.
     *
     * @return The scheme's maximum pure real eigenvalue magnitude.
     * @see diffusive_stability_criterion() for one use of this magnitude.
     */
    virtual typename suzerain::traits::component<Element>::type
        evmaxmag_real() const = 0;

    /**
     * Obtain the scheme's maximum pure imaginary eigenvalue magnitude.
     *
     * @return The scheme's maximum pure imaginary eigenvalue magnitude.
     * @see convective_stability_criterion() for one use of this magnitude.
     */
    virtual typename suzerain::traits::component<Element>::type
        evmaxmag_imag() const = 0;

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
template< typename charT, typename traits, typename Element >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits>& os,
        const ILowStorageMethod<Element>& m)
{
    return os << m.name();
}

/**
 * Encapsulates the three stage, third order scheme from Appendix A of Spalart,
 * Moser, and Rogers' 1991 ``Spectral Methods for the Navier-Stokes Equations
 * with One Infinite and Two Periodic Directions'' published in the
 * <em>Journal of Computational Physics</em> volume 96 pages 297-324.
 */
template< typename Element >
class SMR91Method : public ILowStorageMethod<Element>
{
public:

    /** Explicit default constructor */
    SMR91Method() {};

    /*! @copydoc ILowStorageMethod::name */
    virtual const char * name() const { return "SMR91"; }

    /*! @copydoc ILowStorageMethod::substeps */
    virtual std::size_t substeps() const { return 3; };

    /*! @copydoc ILowStorageMethod::alpha */
    virtual Element alpha(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::beta */
    virtual Element beta(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::gamma */
    virtual Element gamma(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::zeta */
    virtual Element zeta(const std::size_t substep) const;

    /*! @copydoc ILowStorageMethod::evmaxmag_real */
    virtual typename suzerain::traits::component<Element>::type
        evmaxmag_real() const { return 2.512; }

    /*! @copydoc ILowStorageMethod::evmaxmag_imag */
    virtual typename suzerain::traits::component<Element>::type
        evmaxmag_imag() const { return std::sqrt(3.0); }
};

template< typename Element >
Element SMR91Method<Element>::alpha(const std::size_t substep) const
{
    static const Element coeff[3] = { Element( 29)/Element(96),
                                      Element(- 3)/Element(40),
                                      Element(  1)/Element( 6)  };
    return coeff[substep];
}

template< typename Element >
Element SMR91Method<Element>::beta(const std::size_t substep) const
{
    static const Element coeff[3] = { Element(37)/Element(160),
                                      Element( 5)/Element( 24),
                                      Element( 1)/Element(  6)  };
    return coeff[substep];
}

template< typename Element >
Element SMR91Method<Element>::gamma(const std::size_t substep) const
{
    static const Element coeff[3] = { Element(8)/Element(15),
                                      Element(5)/Element(12),
                                      Element(3)/Element( 4)  };
    return coeff[substep];
}

template< typename Element >
Element SMR91Method<Element>::zeta(const std::size_t substep) const
{
    static const Element coeff[3] = { Element(  0),
                                      Element(-17)/Element(60),
                                      Element(- 5)/Element(12)  };
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
 * @param a On entry contains \f$u^{i}\f$ and on exit contains
 *          \f$N\left(u^{i}\right)\f$.
 * @param b On entry contains \f$N\left(u^{i-1}\right)\f$ and on exit contains
 *          \f$u^{i+1}\f$.
 * @param substep_index The substep number to take.
 * @param delta_t The time step \f$\Delta{}t\f$ to take.  The same time step
 *                must be supplied for all substep computations.
 * @return The time step \f$\Delta{}t\f$ taken.
 *         It will always equal \c delta_t.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 * @see The method step() provides more convenient ways to perform multiple
 *      substeps, including dynamic step size computation.
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Storage
>
const typename suzerain::traits::component<Element>::type substep(
    const ILowStorageMethod<Element>& m,
    const ILinearOperator<NumDims,Element,Storage>& L,
    const INonlinearOperator<NumDims,Element,Storage>& N,
    IState<NumDims,Element,Storage>& a,
    IState<NumDims,Element,Storage>& b,
    const typename suzerain::traits::component<Element>::type delta_t,
    const std::size_t substep_index)
throw(std::exception)
{
    if (SUZERAIN_UNLIKELY(substep_index >= m.substeps()))
        throw std::invalid_argument("Requested substep too large");

    b.scale(delta_t * m.zeta(substep_index));
    L.accumulateMassPlusScaledOperator(
            delta_t * m.alpha(substep_index), a, b);
    N.applyOperator(a);
    b.addScaled(delta_t * m.gamma(substep_index), a);
    L.invertMassPlusScaledOperator( -delta_t * m.beta(substep_index), b);

    return delta_t;
}

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu + N(u)\f$.
 * The time step taken, \f$\Delta{}t\f$ will be based on a stable value
 * computed during the first nonlinear operator application as well as an
 * optional fixed maximum step size.
 *
 * @param m The low storage scheme to use.  For example, SMR91Method.
 * @param L The linear operator to be treated implicitly.
 * @param N The nonlinear operator to be treated explicitly.
 * @param a On entry contains \f$u(t)\f$ and on exit contains
 *          \f$u(t+\Delta{}t)\f$.  The linear operator is applied
 *          only to this state storage.
 * @param b Used as a temporary storage location during substeps.
 *          The nonlinear operator is applied only to this state storage.
 * @param max_delta_t An optional maximum time step size.
 * @return The time step \f$\Delta{}t\f$ taken.
 *         It may be less than \c max_delta_t.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
template<
    std::size_t NumDims,
    typename Element,
    typename StorageA,
    typename StorageB
>
const typename suzerain::traits::component<Element>::type step(
    const ILowStorageMethod<Element>& m,
    const ILinearOperator<NumDims,Element,StorageA,StorageB>& L,
    const INonlinearOperator<NumDims,Element,StorageB,StorageA>& N,
    IState<NumDims,Element,StorageA,StorageB>& a,
    IState<NumDims,Element,StorageB,StorageA>& b,
    const typename suzerain::traits::component<Element>::type max_delta_t = 0)
throw(std::exception)
{
    // First substep handling is special since we need to determine delta_t
    b.assign(a);
    typename suzerain::traits::component<Element>::type delta_t
        = N.applyOperator(b, true /* need delta_t */);
    if (delta_t > 0) {
        if (max_delta_t > 0) delta_t = std::min(delta_t, max_delta_t);
    } else if (max_delta_t > 0) {
        delta_t = max_delta_t;
    } else {
        // Both delta_t, max_delta_t are either non-positive or NaN
        throw std::logic_error("No usable delta_t available for time step");
    }
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

    return delta_t;
}

} // namespace lowstorage

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
