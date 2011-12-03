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
#include <suzerain/math.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timecontroller.hpp>
#include <suzerain/traits.hpp>

/** @file
 * Provides time integration-related operator details and advancement schemes.
 */

namespace suzerain
{

/**
 * Provides time integration schemes and the associated interfaces that the
 * underlying operators must obey.
 */
namespace timestepper
{

/**
 * Defines the nonlinear operator interface required for timestepping.
 *
 * @see tparam State A state type which must provide an \c element \c typedef
 *                   containing its real- or complex-valued scalar type.
 * @see suzerain::timestepper::lowstorage for low storage schemes.
 */
template<typename State>
class INonlinearOperator
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename State::element element;

    /** The real-valued scalar corresponding to \c element */
    typedef typename suzerain::traits::component<element>::type component;

    /**
     * Apply the operator in place on a given state location.  That is, compute
     * \f$\mbox{state}\leftarrow{}\mbox{operator}\left(\mbox{state}\right)\f$.
     *
     * @param time The time at which to apply the operator.
     * @param state The state location to use.  It is expected that certain
     *        implementations will require more specific state types.
     * @param evmaxmag_real The timestepping scheme's maximum pure real
     *        eigenvalue magnitude.  When <tt>substep_index == 0</tt>
     *        the operator should use this information to compute a
     *        stable time step per its convective criteria.
     * @param evmaxmag_imag The timestepping scheme's maximum pure imaginary
     *        eigenvalue magnitude.  When <tt>substep_index == 0</tt>
     *        the operator should use this information to compute a
     *        stable time step per its diffusive criteria.
     * @param substep_index The (zero-indexed) current time step substep.
     *        If zero the operator should compute and return stable time
     *        steps according to one or more criteria.
     *
     * @return stable time step sizes according to one or more criteria
     *         if <tt>substep_index == 0</tt>.  The overall stable time
     *         step is the minimum of all elements in the returned vector.
     *         Returning results from multiple criteria allows monitoring
     *         the relative restrictness of each criterion.
     */
    virtual std::vector<component> applyOperator(
            const component time,
            State& state,
            const component evmaxmag_real,
            const component evmaxmag_imag,
            const std::size_t substep_index) const = 0;

    /** Virtual destructor for peace of mind. */
    virtual ~INonlinearOperator() {}
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
 * For formulations in which an explicit Mach number
 * \f$\mbox{Ma}=\frac{u_0}{a_0}\f$ appears, one \em must provide the velocities
 * and the sound speed \em both nondimensionalized using \f$u_0\f$.  That
 * expressions like \f$\left|u\right| + \frac{a}{\mbox{Ma}}\f$ are appropriate
 * in that context can be seen by finding the eigenvalues of the Euler
 * equations in such a nondimensionalization.
 *
 * Using a hybrid implicit/explicit %timestepper with acoustic terms computed
 * implicitly effectively sets the sound speed to be zero for this CFL
 * calculation.
 *
 * @note Prem Venugopal's 2003 thesis used a nearly identical conservative
 * stability criterion (equation 3.10).  Venugopal found the constraint to be
 * overly conservative in the wall-normal direction because the derivation
 * assumed periodicity.  In section 3.2 he presents a linearized analysis
 * taking into account the inhomogeneous nature of the wall-normal direction.
 * He determined that the wall-normal imaginary eigenvalue magnitude dropped by
 * nearly an order of magnitude after taking into account the inhomogeneity.
 * He concluded that using an effective \f$1/\Delta{}y\f$ <i>four<i> times
 * smaller than the nominal value was feasible (equation 3.29).  His approach
 * can be accomplished by specifying <tt>one_over_delta_y / 4</tt> when
 * invoking this method.
 *
 * @param u_x              Velocity in the X direction \f$u_{x}\f$
 * @param one_over_delta_x Inverse local X grid spacing \f$1/\Delta{}x\f$
 * @param u_y              Velocity in the Y direction \f$u_{y}\f$
 * @param one_over_delta_y Inverse local Y grid spacing \f$1/\Delta{}y\f$
 * @param u_z              Velocity in the Z direction \f$u_{z}\f$
 * @param one_over_delta_z Inverse local Z grid spacing \f$1/\Delta{}z\f$
 * @param evmaxmag_imag    The maximum pure imaginary eigenvalue magnitude
 *                         for some Runge-Kutta scheme, denoted \f$\left|
 *                         \lambda_{I}\Delta_{}t \right|_{\mbox{max}}\f$
 *                         in Guarini's thesis.
 * @param a                The local sound speed \f$a\f$.
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
    static const FPT one_over_pi
        = (FPT) 0.3183098861837906715377675267450287L;

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
 * (1998)
 * \f[
 *   \mbox{max}\!\left(
 *     \left|\frac{\gamma\left(\nu-\nu_{0}\right)}{\mbox{Re}\mbox{Pr}}\right|,
 *     \left|\frac{\nu-\nu_{0}}{\mbox{Re}}\right|,
 *     \left|\frac{\nu_{B}-\nu_{B0}}{\mbox{Re}}\right|
 *   \right)
 *   \pi^{2}
 *   \left(
 *       \frac{1}{\Delta{}x^{2}}
 *     + \frac{1}{\Delta{}y^{2}}
 *     + \frac{1}{\Delta{}z^{2}}
 *   \right)
 *   \Delta{}t \leq \left|\lambda_{R}\Delta_{}t\right|_{\mbox{max}}
 * \f]
 * where a bulk kinematic viscosity \f$\nu_{B}\f$ has
 * been added.  The maximum pure real eigenvalue magnitude,
 * \f$\left|\lambda_{R}\Delta{}t\right|_{\mbox{max}}\f$, is a feature
 * of the chosen timestepping method.  For example, it is 2.512 for
 * the SMR91 scheme.  The absolute values within the maximum operation
 * account for the possibility that \f$\nu<\nu_{0}\f$.
 *
 * Using a hybrid implicit/explicit %timestepper with viscous terms computed
 * implicitly sets \f$\nu_0\f$ to be the reference kinematic viscosity about
 * which the viscous terms were linearized.
 *
 * @param one_over_delta_x Inverse local X grid spacing \f$1/\Delta{}x\f$
 * @param one_over_delta_y Inverse local Y grid spacing \f$1/\Delta{}y\f$
 * @param one_over_delta_z Inverse local Z grid spacing \f$1/\Delta{}z\f$
 * @param Re               The Reynolds number \f$\mbox{Re}\f$
 * @param Pr               The Prandtl number \f$\mbox{Pr}\f$
 * @param gamma            The ratio of specific heats \f$\gamma\f$
 * @param evmaxmag_real    The maximum pure real eigenvalue magnitude
 *                         for some Runge-Kutta scheme, denoted \f$\left|
 *                         \lambda_{R}\Delta_{}t \right|_{\mbox{max}}\f$
 *                         in Guarini's thesis.
 * @param nu               The local kinematic viscosity \f$\nu\f$
 * @param nu0              The kinematic viscosity reference value
 *                         \f$\nu_{0}\f$
 * @param nuB              The local bulk kinematic viscosity \f$\nu\f$
 * @param nuB0             The bulk kinematic viscosity reference value
 *                         \f$\nu_{0}\f$
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
        const FPT nu0 = 0,
        const FPT nuB = 0,
        const FPT nuB0 = 0)
{
    // Precision for a 128-bit quad found via Sage's N(1/(pi*pi),digits=34)
    static const FPT one_over_pi_squared
        = (FPT) 0.1013211836423377714438794632097276L;

    // Kinematic viscosity and bulk kinematic viscosity enter identically
    const FPT nu_less_nu0 = std::max(std::abs(nu - nu0), std::abs(nuB - nuB0));

    const FPT maxcoeff = std::max((gamma*nu_less_nu0)/(Re*Pr), nu_less_nu0/Re);
    return (evmaxmag_real * one_over_pi_squared)
        /  (maxcoeff * (   one_over_delta_x*one_over_delta_x
                         + one_over_delta_y*one_over_delta_y
                         + one_over_delta_z*one_over_delta_z));
}

namespace { // anonymous

/** A comparator which considers NaNs to be less than all other values */
struct NanIsMinimumComparator {
    template< typename T >
    bool operator()(const T &x, const T &y) const {
        return SUZERAIN_UNLIKELY((boost::math::isnan)(x)) || x < y;
    }
};

} // end namespace anonymous

/**
 * A functor that finds the minimum stable time step given a vector of
 * candidates.  Any NaN appearing in the candidates causes a NaN to be
 * returned.  The candidates may not be empty.
 */
struct DeltaTReducer {
    template< typename T>
    T operator()(const std::vector<T> & candidates) {
        assert(candidates.size() > 0);
        return *std::min_element(candidates.begin(), candidates.end(),
                                 NanIsMinimumComparator());
    }
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
 * @see tparam StateA A state type which must provide an \c element \c typedef
 *                    containing its real- or complex-valued scalar type.
 * @see tparam StateB A state type which must provide an \c element \c typedef
 *                    containing its real- or complex-valued scalar type.
 */
template<typename StateA, typename StateB = StateA>
class ILinearOperator
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename StateA::element element;
    BOOST_STATIC_ASSERT(
            (boost::is_same<element,typename StateB::element>::value));

    /** The real-valued scalar corresponding to \c element */
    typedef typename suzerain::traits::component<element>::type component;

    /**
     * Apply \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)\mbox{state}\f$.
     *
     * @param phi Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @param delta_t The size of the currently active time step.
     * @param substep_index The (zero-indexed) time stepper substep index.
     */
    virtual void applyMassPlusScaledOperator(
            const element& phi,
            StateA& state,
            const component delta_t,
            const std::size_t substep_index) const = 0;

    /**
     * Accumulate \f$M+\phi{}L\f$ out-of-place for scalar \f$\phi\f$.
     * That is, \f$\mbox{output}\leftarrow{}
     * \left(M+\phi{}L\right)\mbox{input}+\beta\mbox{output}\f$.
     *
     * @param phi Scale factor \f$\phi\f$ to use.
     * @param input State vector on which to apply the operator.
     * @param beta  Scale factor for output vector during accumulation.
     * @param output State vector into which to accumulate the result.
     * @param delta_t The size of the currently active time step.
     * @param substep_index The (zero-indexed) time stepper substep index.
     */
    virtual void accumulateMassPlusScaledOperator(
            const element& phi,
            const StateA& input,
            const element& beta,
            StateB& output,
            const component delta_t,
            const std::size_t substep_index) const = 0;

    /**
     * Invert \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.
     * That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)^{-1}\mbox{state}\f$.
     *
     * @param phi Scale factor \f$\phi\f$ to use.
     * @param state State vector on which to apply the operator.
     * @param delta_t The size of the currently active time step.
     * @param substep_index The (zero-indexed) time stepper substep index.
     * @param iota The \f$iota_i\f$ value appropriate for \c substep_index.
     *             See \ref ILowStorageMethod for details on how to use \c iota.
     */
    virtual void invertMassPlusScaledOperator(
            const element& phi,
            StateA& state,
            const component delta_t,
            const std::size_t substep_index,
            const component iota) const = 0;

    /** Virtual destructor for peace of mind. */
    virtual ~ILinearOperator() {}
};

/**
 * Implements simple multiplicative operator which scales all state
 * variables by a uniform factor.  The associated mass matrix \f$M\f$
 * is the identity matrix.
 *
 * @tparam StateA A state type which should descend from suzerain::StateBase.
 * @tparam StateB A state type which should descend from suzerain::StateBase.
 */
template<typename StateA,typename StateB = StateA>
class MultiplicativeOperator
    : public ILinearOperator<StateA,StateB>,
      public INonlinearOperator<StateB>
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename StateA::element element;
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename StateB::element>::value));

    /** The real-valued scalar corresponding to \c element */
    typedef typename suzerain::traits::component<element>::type component;

    /**
     * Construct an instance which scales by \c factor and reports
     * \c delta_t as a stable time step.
     *
     * @param factor uniform scaling factor to apply.
     * @param delta_t uniform, presumably stable time step to
     *        always return from ::applyOperator.
     */
    template< typename FactorType, typename DeltaTType >
    MultiplicativeOperator(const FactorType& factor, const DeltaTType& delta_t)
        : factor(factor), delta_t(delta_t) {}

    /**
     * Construct an instance which scales by \c factor and reports the
     * maximum representable floating point value as a stable time step.
     * Useful in testing contexts or during operator composition.
     *
     * @param factor uniform scaling factor to apply.
     */
    template< typename FactorType >
    MultiplicativeOperator(const FactorType& factor)
        : factor(factor),
          delta_t(std::numeric_limits<component>::infinity()) {}

    /**
     * Scale \c state by the factor set at construction time.
     *
     * @param time The time at which to apply the operator (ignored).
     * @param state to scale in place.
     * @param evmaxmag_real Ignored in this implementation.
     * @param evmaxmag_imag Ignored in this implementation.
     * @param substep_index Ignored in this implementation.
     *
     * @return The \c delta_t provided at construction time.
     */
    virtual std::vector<component> applyOperator(
            const component time,
            StateB& state,
            const component evmaxmag_real,
            const component evmaxmag_imag,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(evmaxmag_real);
        SUZERAIN_UNUSED(evmaxmag_imag);
        SUZERAIN_UNUSED(substep_index);
        state.scale(factor);
        return std::vector<component>(1, delta_t);
    }

    /**
     * Compute \f$\mbox{state}\leftarrow{}
     * \left(I+\phi\times\mbox{factor}\right) \mbox{state}\f$ where \c factor
     * is the scaling factor set at construction time.
     *
     * @param phi Additional scaling \f$\phi\f$ to apply.
     * @param state to modify in place.
     * @param delta_t Ignored.
     * @param substep_index Ignored.
     */
    virtual void applyMassPlusScaledOperator(
            const element& phi,
            StateA& state,
            const component delta_t = 0,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(delta_t);
        SUZERAIN_UNUSED(substep_index);
        state.scale(phi*factor + element(1));
    }

    /**
     * Compute \f$\mbox{output}\leftarrow{}\mbox{output}+
     * \left(I+\phi\times\mbox{factor}\right) \mbox{input}\f$ where \c factor
     * is the scaling factor set at construction time.
     *
     * @param phi Additional scaling \f$\phi\f$ to apply.
     * @param input on which to apply the operator.
     * @param beta  Scale factor for output vector during accumulation.
     * @param output on which to accumulate the result.
     * @param delta_t Ignored.
     * @param substep_index Ignored.
     */
    virtual void accumulateMassPlusScaledOperator(
            const element& phi,
            const StateA& input,
            const element& beta,
            StateB& output,
            const component delta_t = 0,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(delta_t);
        SUZERAIN_UNUSED(substep_index);
        output.scale(beta);
        output.addScaled(phi*factor + element(1), input);
    }

    /**
     * Compute \f$\mbox{state}\leftarrow{}
     * \left(I+\phi\times\mbox{factor}\right)^{-1} \mbox{state}\f$ where \c
     * factor is the scaling factor set at construction time.
     *
     * @param phi Additional scaling \f$\phi\f$ to apply.
     * @param state to modify in place.
     * @param delta_t Ignored.
     * @param substep_index Ignored.
     */
    virtual void invertMassPlusScaledOperator(
            const element& phi,
            StateA& state,
            const component delta_t = 0,
            const std::size_t substep_index = 0,
            const component iota = 0) const
    {
        SUZERAIN_UNUSED(delta_t);
        SUZERAIN_UNUSED(substep_index);
        SUZERAIN_UNUSED(iota);
        state.scale((element(1))/(phi*factor + element(1)));
    }

private:
    /** Uniform scale factor to apply as part of operator */
    const element factor;

    /** Uniform stable time step to report */
    const component delta_t;
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
    typename A,
    typename B
>
MultiplicativeOperator<A,B>
make_multiplicator_operator(
    const FactorType& factor,
    const DeltaTType& delta_t,
    const StateBase<A>& input,
    const StateBase<B>& output)
{
    MultiplicativeOperator<A,B> retval(factor,delta_t);
    return retval;
}

// TODO ILowStorageMethod could employ the CRTP to avoid virtual overhead

/**
 * Encapsulates a hybrid implicit/explicit low storage Runge-Kutta method to
 * advance \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$.  The method consists of one or
 * more substeps governed by coefficients \f$\alpha_i\f$, \f$\beta_i\f$,
 * \f$\gamma_i\f$, \f$\zeta_i\f$, and \f$\eta_i\f$ where \f$i\f$ is less than
 * the number of substeps.  Each substep obeys
 * \f[
 *   \left(M - \Delta{}t\beta_{i}L\right) u^{i+1}
 *   =
 *   \left(M + \Delta{}t\alpha_{i}L\right) u^{i}
 *   + \Delta{}t\gamma_{i}
 *     N\left(u^{i},t + \eta_{i}\Delta{}t\right)
 *   + \Delta{}t\zeta_{i}
 *     N\left(u^{i-1}, t + \eta_{i}\Delta{}t\right)
 * \f]
 * where \f$\alpha_i+\beta_i=\gamma_i+\zeta_i\f$ and \f$\eta_{i} =
 * \sum_{i=0}^{i-1} \alpha_{i} + \beta{i} \f$. Note that the indexing on the
 * \f$\zeta_i\f$ coefficients differs slightly from other sources.  Note also
 * that not all literature sources include the possibility of a time-dependent
 * \f$N\f$ operator.
 *
 * The <tt>i</tt>th substep has a "duration" \f$d_i = \left(\eta_{i+1} -
 * \eta_i\right)\Delta{}t\f$ where \f$i+1\f$ greater than the number of
 * substeps is treated as 1.  A running mean quantity across \f$N\f$ substeps,
 * denoted \f$\bar{q}_{N-1}\f$, may be accumulated via
 * \f[
 *   \bar{q}_{N-1} = \frac{\sum_{i=0}^{N-1} d_{i} q_{i}}
 *                        {\sum_{i=0}^{N-1} d_{i}      }.
 * \f]
 * Rearranging the result in terms of means from earlier substeps,
 * \f[
 *   \bar{q}_{N-1} = \frac{\sum_{i=0}^{N-2}}
 *                        {d_{N-1} + \sum_{i=0}^{N-2} d_i} \bar{q}_{N-2}
 *                 + \frac{d_{N-1}}
 *                        {d_{N-1} + \sum_{i=0}^{N-2} d_i} q_{N-1}
 * \f]
 * motivates defining
 * \f[
 *   \iota_{N-1} = \frac{d_i}{\sum_{i=0}^{N-1} d_i}
 *               = \frac{\eta_{i+1} - \eta_i}{\eta_{i+1} - \eta_{0}}
 * \f]
 * since the factor \f$\Delta{}t\f$ may be omitted without changing
 * \f$\bar{q}_{N-1}\f$.  Then maintaining the running mean
 * \f[
 *   \bar{q}_{N-1} = \left(1 - \iota_{N-1}\right) \bar{q}_{N-2}
 *                 + \iota_{N-1} q_{N-1}
 * \f]
 * may be done via a simple update operation <tt>mean += iota_i * (sample -
 * mean)</tt>.  One use case for \f$\iota_i\f$ is obtaining running means of
 * implicitly computed quantities using minimal space and time overhead.
 *
 * @see ILinearOperator for the interface that \f$L\f$ must implement.
 * @see INonlinearOperator for the interface that \f$N\f$ must implement.
 * @see SMR91Method and Yang11Method for examples of concrete schemes.
 * @see step() or substep() for methods that can advance state variables
 *      according to a timestepping method.
 */
template<typename Element>
class ILowStorageMethod
{
public:

    /** The real-valued scalar corresponding to \c Element */
    typedef typename suzerain::traits::component<Element>::type component;

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
    virtual component alpha(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\beta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component beta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\gamma_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component gamma(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\zeta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component zeta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\eta_i\f$ coefficient, which is
     * derived from other scheme coefficients.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component eta(std::size_t substep) const = 0;

    /**
     * Compute the scheme's derived \f$\iota_i\f$ coefficient, which is used to
     * accumulate a running time-averaged value across substeps.
     *
     * @param substep A substep number \f$i\f$ in the range [0,::substeps).
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component iota(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's maximum pure real eigenvalue magnitude.
     *
     * @return The scheme's maximum pure real eigenvalue magnitude.
     * @see diffusive_stability_criterion() for one use of this magnitude.
     */
    virtual component evmaxmag_real() const = 0;

    /**
     * Obtain the scheme's maximum pure imaginary eigenvalue magnitude.
     *
     * @return The scheme's maximum pure imaginary eigenvalue magnitude.
     * @see convective_stability_criterion() for one use of this magnitude.
     */
    virtual component evmaxmag_imag() const = 0;

    /** Virtual destructor to support interface-like behavior. */
    virtual ~ILowStorageMethod() {}
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

    /** The real-valued scalar corresponding to \c Element */
    typedef typename ILowStorageMethod<Element>::component component;

    /**
     * Explicit constructor.
     *
     * @param evmagfactor The multiplicative factor to use when reporting
     *                    maximum pure real and pure imaginary eigenvalue
     *                    magnitudes in evmaxmag_real() and evmaxmag_imag(),
     *                    respectively.
     */
    explicit SMR91Method(component evmagfactor = 1)
        : evmaxmag_real_(evmagfactor * component(2.51274532661832862402373L)),
          evmaxmag_imag_(evmagfactor * std::sqrt(component(3)))
    {
        assert(evmagfactor > 0);
    }

    /** @copydoc ILowStorageMethod::name */
    virtual const char * name() const
    {
        return "SMR91";
    }

    /** @copydoc ILowStorageMethod::substeps */
    virtual std::size_t substeps() const
    {
        return 3;
    }

    /** @copydoc ILowStorageMethod::alpha */
    virtual component alpha(const std::size_t substep) const
    {
        static const component coeff[3] = { component( 29)/component(96),
                                            component(- 3)/component(40),
                                            component(  1)/component( 6)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::beta */
    virtual component beta(const std::size_t substep) const
    {
        static const component coeff[3] = { component(37)/component(160),
                                            component( 5)/component( 24),
                                            component( 1)/component(  6)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::gamma */
    virtual component gamma(const std::size_t substep) const
    {
        static const component coeff[3] = { component(8)/component(15),
                                            component(5)/component(12),
                                            component(3)/component( 4)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::zeta */
    virtual component zeta(const std::size_t substep) const
    {
        static const component coeff[3] = { component(  0),
                                            component(-17)/component(60),
                                            component(- 5)/component(12)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::eta */
    virtual component eta(const std::size_t substep) const
    {
        static const component coeff[3] = { component(0),
                                            component(8)/component(15),
                                            component(2)/component( 3)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::iota */
    virtual component iota(const std::size_t substep) const
    {
        static const component coeff[3] = { component(1),
                                            component(1)/component(5),
                                            component(1)/component(3)  };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::evmaxmag_real */
    virtual component evmaxmag_real() const
    {
        return evmaxmag_real_;
    }

    /** @copydoc ILowStorageMethod::evmaxmag_imag */
    virtual component evmaxmag_imag() const {
        return evmaxmag_imag_;
    }

private:

    /** Value to report from evmaxmag_real(). */
    component evmaxmag_real_;

    /** Value to report from evmaxmag_imag(). */
    component evmaxmag_imag_;
};

/**
 * Encapsulates the three stage, second-order, adjoint-consistent scheme from
 * Shan Yang's 2011 thesis ``A shape Hessian based analysis of roughness
 * effects on fluid flows''.
 */
template< typename Element >
class Yang11Method : public ILowStorageMethod<Element>
{
public:

    /** The real-valued scalar corresponding to \c Element */
    typedef typename ILowStorageMethod<Element>::component component;

    /**
     * Explicit constructor.
     *
     * @param evmagfactor The multiplicative factor to use when reporting
     *                    maximum pure real and pure imaginary eigenvalue
     *                    magnitudes in evmaxmag_real() and evmaxmag_imag(),
     *                    respectively.
     */
    explicit Yang11Method(component evmagfactor = 1)
        : evmaxmag_real_(evmagfactor * component(2.51274532661832862402373L)),
          evmaxmag_imag_(evmagfactor * std::sqrt(component(3)))
    {
        assert(evmagfactor > 0);
    }

    /*! @copydoc ILowStorageMethod::name */
    virtual const char * name() const
    {
        return "Yang11";
    }

    /*! @copydoc ILowStorageMethod::substeps */
    virtual std::size_t substeps() const
    {
        return 3;
    }

    /*! @copydoc ILowStorageMethod::alpha */
    virtual component alpha(const std::size_t substep) const
    {
        static const component coeff[3] = { component( 1)/component(3),
                                            component(-1)/component(2),
                                            component( 1)/component(3)  };
        return coeff[substep];
    }

    /*! @copydoc ILowStorageMethod::beta */
    virtual component beta(const std::size_t substep) const
    {
        static const component coeff[3] = { component(1)/component(6),
                                            component(2)/component(3),
                                            component(0)               };
        return coeff[substep];
    }

    /*! @copydoc ILowStorageMethod::gamma */
    virtual component gamma(const std::size_t substep) const
    {
        static const component coeff[3] = { component(1)/component(2),
                                            component(1)/component(3),
                                            component(1)               };
        return coeff[substep];
    }

    /*! @copydoc ILowStorageMethod::zeta */
    virtual component zeta(const std::size_t substep) const
    {
        static const component coeff[3] = { component( 0),
                                            component(-1)/component(6),
                                            component(-2)/component(3)  };
        return coeff[substep];
    }

    /*! @copydoc ILowStorageMethod::eta */
    virtual component eta(const std::size_t substep) const
    {
        static const component coeff[3] = { component(0),
                                            component(1)/component(2),
                                            component(2)/component(3)   };
        return coeff[substep];
    }

    /** @copydoc ILowStorageMethod::iota */
    virtual component iota(const std::size_t substep) const
    {
        static const component coeff[3] = { component(1),
                                            component(1)/component(4),
                                            component(1)/component(3)  };
        return coeff[substep];
    }

    /*! @copydoc ILowStorageMethod::evmaxmag_real */
    virtual component evmaxmag_real() const
    {
        return evmaxmag_real_;
    }

    /*! @copydoc ILowStorageMethod::evmaxmag_imag */
    virtual component evmaxmag_imag() const {
        return evmaxmag_imag_;
    }

private:

    /**< Value to report from evmaxmag_real(). */
    component evmaxmag_real_;

    /**< Value to report from evmaxmag_imag(). */
    component evmaxmag_imag_;
};


/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme to advance the system \f$ M u_t = Lu +
 * \chi N(u,t) \f$.  Note that the roles of state locations \c a and \c b are
 * flipped after each substep.
 *
 * @param m The low storage scheme to use.  For example, SMR91Method.
 * @param L The linear operator to be treated implicitly.
 * @param chi The factor \f$\chi\f$ used to scale the nonlinear operator.
 * @param N The nonlinear operator to be treated explicitly.
 * @param time The simulation time \f$t\f$ during which the substep occurs.
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
template< typename Element, typename A, typename B >
const typename suzerain::traits::component<Element>::type substep(
    const ILowStorageMethod<Element>& m,
    const ILinearOperator<A,B>& L,
    const typename suzerain::traits::component<Element>::type chi,
    const INonlinearOperator<B>& N,
    const typename suzerain::traits::component<Element>::type time,
    A& a,  // FIXME: Would love a StateBase subclass restriction here
    B& b,  // FIXME: Would love a StateBase subclass restriction here
    const typename suzerain::traits::component<Element>::type delta_t,
    const std::size_t substep_index)
{
    BOOST_STATIC_ASSERT((boost::is_same<Element,typename A::element>::value));
    BOOST_STATIC_ASSERT((boost::is_same<Element,typename B::element>::value));

    if (SUZERAIN_UNLIKELY(substep_index >= m.substeps()))
        throw std::invalid_argument("Requested substep too large");

    L.accumulateMassPlusScaledOperator(
                  delta_t * m.alpha(substep_index), a,
            chi * delta_t * m.zeta(substep_index),  b,
            delta_t, substep_index);
    N.applyOperator(time + delta_t * m.eta(substep_index), a,
                    m.evmaxmag_real(), m.evmaxmag_imag(), substep_index);
    b.addScaled(chi * delta_t * m.gamma(substep_index), a);
    L.invertMassPlusScaledOperator(-delta_t * m.beta(substep_index), b,
                                   delta_t, substep_index,
                                   m.iota(substep_index));

    return delta_t;
}

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu +
 * N(u,t)\f$.  The time step taken, \f$\Delta{}t\f$, will be based on a stable
 * value computed during the first nonlinear operator application as well as an
 * optional fixed maximum step size.
 *
 * @param m       The low storage scheme to use.  For example, SMR91Method.
 * @param reducer A stateful functor taking a vector of stable time step
 *                candidates down to a single stable time step.  Users may
 *                employ a custom functor compatible with DeltaTReducer to add
 *                logic for monitoring or manipulating stable step criteria.
 * @param L       The linear operator to be treated implicitly.
 * @param chi     The factor \f$\chi\f$ used to scale the nonlinear operator.
 * @param N       The nonlinear operator to be treated explicitly.
 * @param time    The simulation time \f$t\f$ at which to take the step.
 * @param a       On entry contains \f$u(t)\f$ and on exit contains
 *                \f$u(t+\Delta{}t)\f$.  The linear operator is applied
 *                only to this state storage.
 * @param b       Used as a temporary storage location during substeps.
 *                The nonlinear operator is applied only to this state storage.
 * @param max_delta_t An optional maximum time step size.
 *
 * @return The time step \f$\Delta{}t\f$ taken.
 *         It may be less than \c max_delta_t.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
template< typename Element, typename A, typename B, typename Reducer >
const typename suzerain::traits::component<Element>::type step(
    const ILowStorageMethod<Element>& m,
    Reducer& reducer,
    const ILinearOperator<A,B>& L,
    const typename suzerain::traits::component<Element>::type chi,
    const INonlinearOperator<B>& N,
    const typename suzerain::traits::component<Element>::type time,
    A& a,  // FIXME: Would love a StateBase subclass restriction here
    B& b,  // FIXME: Would love a StateBase subclass restriction here
    const typename suzerain::traits::component<Element>::type max_delta_t = 0)
{
    BOOST_STATIC_ASSERT((boost::is_same<Element,typename A::element>::value));
    BOOST_STATIC_ASSERT((boost::is_same<Element,typename B::element>::value));
    typedef typename suzerain::traits::component<Element>::type component_type;

    // First substep handling is special since we need to determine delta_t
    b.assign(a);
    const std::vector<component_type> delta_t_candidates
        = N.applyOperator(time, b, m.evmaxmag_real(), m.evmaxmag_imag(), 0);
    component_type delta_t = reducer(delta_t_candidates);

    if (max_delta_t > 0) {
        delta_t = suzerain::math::minnan(delta_t, max_delta_t);
    }
    L.applyMassPlusScaledOperator(delta_t * m.alpha(0), a, delta_t, 0);
    a.addScaled(chi * delta_t * m.gamma(0), b);
    L.invertMassPlusScaledOperator(-delta_t * m.beta(0), a,
                                   delta_t, 0, m.iota(0));

    // Second and subsequent substeps are identical
    for (std::size_t i = 1; i < m.substeps(); ++i) {
        L.accumulateMassPlusScaledOperator(
                delta_t * m.alpha(i),      a,
                chi * delta_t * m.zeta(i), b,
                delta_t, i);
        b.exchange(a); // Note nonlinear storage controls exchange operation
        N.applyOperator(time + delta_t * m.eta(i), b,
                        m.evmaxmag_real(), m.evmaxmag_imag(), i);
        a.addScaled(chi * delta_t * m.gamma(i), b);
        L.invertMassPlusScaledOperator(-delta_t * m.beta(i), a,
                                       delta_t, i, m.iota(i));
    }

    return delta_t;
}

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu +
 * N(u,t)\f$.  The time step taken, \f$\Delta{}t\f$, will be based on a stable
 * value computed during the first nonlinear operator application as well as an
 * optional fixed maximum step size.
 *
 * @param m       The low storage scheme to use.  For example, SMR91Method.
 * @param L       The linear operator to be treated implicitly.
 * @param chi     The factor \f$\chi\f$ used to scale the nonlinear operator.
 * @param N       The nonlinear operator to be treated explicitly.
 * @param time    The simulation time \f$t\f$ at which to take the step.
 * @param a       On entry contains \f$u(t)\f$ and on exit contains
 *                \f$u(t+\Delta{}t)\f$.  The linear operator is applied
 *                only to this state storage.
 * @param b       Used as a temporary storage location during substeps.
 *                The nonlinear operator is applied only to this state storage.
 * @param max_delta_t An optional maximum time step size.
 *
 * @return The time step \f$\Delta{}t\f$ taken.
 *         It may be less than \c max_delta_t.
 *
 * @see ILowStorageMethod for the equation governing time advancement.
 */
template< typename Element, typename A, typename B >
const typename suzerain::traits::component<Element>::type step(
    const ILowStorageMethod<Element>& m,
    const ILinearOperator<A,B>& L,
    const typename suzerain::traits::component<Element>::type chi,
    const INonlinearOperator<B>& N,
    const typename suzerain::traits::component<Element>::type time,
    A& a,  // FIXME: Would love a StateBase subclass restriction here
    B& b,  // FIXME: Would love a StateBase subclass restriction here
    const typename suzerain::traits::component<Element>::type max_delta_t = 0)
{
    DeltaTReducer reducer;
    return step<Element, A, B, DeltaTReducer>(m, reducer, L, chi, N,
                                              time, a, b, max_delta_t);
}

/**
 * Provides higher-level control mechanisms built atop low storage time
 * integration schemes.  This class is a thin wrapper combining TimeController
 * with step().
 *
 * @see TimeController for details on the time controller logic.
 * @see make_LowStorageTimeController for an easy way to create
 *      an instance with the appropriate type signature.
 */
template< typename A, typename B, typename Reducer >
class LowStorageTimeController
    : public TimeController< typename suzerain::traits::component<
            typename A::element
      >::type >
{
protected:

    /** Shorthand for the superclass */
    typedef TimeController< typename suzerain::traits::component<
                typename A::element
            >::type > super;

public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename A::element element;
    BOOST_STATIC_ASSERT((boost::is_same<element, typename B::element>::value));

    /**
     * Construct an instance that will advance a simulation built atop the
     * given operators and storage.
     *
     * @param m         The low storage scheme to use.
     *                  For example, SMR91Method.
     * @param reducer   A stateful functor taking a vector of stable time step
     *                  candidates down to a single stable time step.  Users
     *                  may employ a custom functor compatible with
     *                  DeltaTReducer to add logic for monitoring or
     *                  manipulating stable step criteria.
     * @param L         The linear operator to be treated implicitly.
     * @param chi       The factor \f$\chi\f$ used to scale the nonlinear
     *                  operator.
     * @param N         The nonlinear operator to be treated explicitly.
     * @param a         On entry contains \f$u(t)\f$ and on exit contains
     *                  \f$u(t+\Delta{}t)\f$.  The linear operator is applied
     *                  only to this state storage.
     * @param b         Used as a temporary storage location during substeps.
     *                  The nonlinear operator is applied only to this state
     *                  storage.
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::epsilon()</tt>.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::max()</tt>.
     *                  See max_dt() for the associated semantics.
     *
     * @see The method step() for more details on
     *      \c m, \c reducer, \c L, \c N, \c a, and \c b.
     */
    LowStorageTimeController(
            const ILowStorageMethod<element>& m,
            Reducer& reducer,
            const ILinearOperator<A,B>& L,
            const typename suzerain::traits::component<element>::type chi,
            const INonlinearOperator<B>& N,
            A& a,  // FIXME: Would love a StateBase subclass restriction here
            B& b,  // FIXME: Would love a StateBase subclass restriction here
            typename super::time_type initial_t = 0,
            typename super::time_type min_dt = 0,
            typename super::time_type max_dt = 0)
        : super(boost::bind(&LowStorageTimeController::stepper, this, _1),
                initial_t,
                min_dt,
                max_dt),
          m(m), reducer(reducer), L(L), chi(chi), N(N), a(a), b(b) {}

private:

    const ILowStorageMethod<element>& m;
    Reducer &reducer;
    const ILinearOperator<A,B>& L;
    const typename suzerain::traits::component<element>::type chi;
    const INonlinearOperator<B>& N;
    A& a;
    B& b;

    typename super::time_type stepper(typename super::time_type max_dt)
    {
        return suzerain::timestepper::lowstorage::step(
                m, reducer, L, chi, N, super::current_t(), a, b, max_dt);
    }

};

/**
 * A partial specialization of the LowStorageTimeController template for the
 * case when default DeltaTReducer behavior is desired.
 */
template< typename A, typename B >
class LowStorageTimeController<A, B, void>
    : private DeltaTReducer,
      public LowStorageTimeController<A,B,DeltaTReducer>
{

protected:

    typedef typename LowStorageTimeController<A,B,DeltaTReducer>::super super;

public:

    typedef typename LowStorageTimeController<A,B,DeltaTReducer>::element element;

    /**
     * Construct an instance that will advance a simulation built atop the
     * given operators and storage.
     *
     * @param m         The low storage scheme to use.
     *                  For example, SMR91Method.
     * @param L         The linear operator to be treated implicitly.
     * @param chi       The factor \f$\chi\f$ used to scale the nonlinear
     *                  operator.
     * @param N         The nonlinear operator to be treated explicitly.
     * @param a         On entry contains \f$u(t)\f$ and on exit contains
     *                  \f$u(t+\Delta{}t)\f$.  The linear operator is applied
     *                  only to this state storage.
     * @param b         Used as a temporary storage location during substeps.
     *                  The nonlinear operator is applied only to this state
     *                  storage.
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::epsilon()</tt>.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::max()</tt>.
     *                  See max_dt() for the associated semantics.
     *
     * @see The method step() for more details on
     *      \c m, \c reducer, \c L, \c N, \c a, and \c b.
     */
    LowStorageTimeController(
            const ILowStorageMethod<element>& m,
            const ILinearOperator<A,B>& L,
            const typename suzerain::traits::component<element>::type chi,
            const INonlinearOperator<B>& N,
            A& a,  // FIXME: Would love a StateBase subclass restriction here
            B& b,  // FIXME: Would love a StateBase subclass restriction here
            typename super::time_type initial_t = 0,
            typename super::time_type min_dt = 0,
            typename super::time_type max_dt = 0)
        : DeltaTReducer(),
          LowStorageTimeController<A,B,DeltaTReducer>(
                  m,
                  *reinterpret_cast<DeltaTReducer*>(this),
                  L, chi, N, a, b, initial_t, min_dt, max_dt)
    {}

};

/**
 * A helper method so the compiler can deduce the appropriate template
 * types for a LowStorageTimeController.
 *
 * \copydoc #LowStorageTimeController
 */
template< typename A, typename B, typename ChiType, typename Reducer >
LowStorageTimeController<A,B,Reducer>*
make_LowStorageTimeController(
        const ILowStorageMethod<typename StateBase<A>::element>& m,
        Reducer &reducer,
        const ILinearOperator<A,B>& L,
        const ChiType chi,
        const INonlinearOperator<B>& N,
        A& a,  // FIXME: Would love a StateBase subclass restriction here
        B& b,  // FIXME: Would love a StateBase subclass restriction here
        typename LowStorageTimeController<A,B,Reducer>::time_type initial_t = 0,
        typename LowStorageTimeController<A,B,Reducer>::time_type min_dt = 0,
        typename LowStorageTimeController<A,B,Reducer>::time_type max_dt = 0)
{
    return new LowStorageTimeController<A,B,Reducer>(
            m, reducer, L, chi, N, a, b, initial_t, min_dt, max_dt);
}

/**
 * A helper method so the compiler can deduce the appropriate template
 * types for a LowStorageTimeController.
 *
 * \copydoc #LowStorageTimeController
 */
template< typename A, typename B, typename ChiType >
LowStorageTimeController<A,B,void>*
make_LowStorageTimeController(
        const ILowStorageMethod<typename StateBase<A>::element>& m,
        const ILinearOperator<A,B>& L,
        const ChiType chi,
        const INonlinearOperator<B>& N,
        A& a,  // FIXME: Would love a StateBase subclass restriction here
        B& b,  // FIXME: Would love a StateBase subclass restriction here
        typename LowStorageTimeController<A,B,void>::time_type initial_t = 0,
        typename LowStorageTimeController<A,B,void>::time_type min_dt = 0,
        typename LowStorageTimeController<A,B,void>::time_type max_dt = 0)
{
    return new LowStorageTimeController<A,B,void>(
            m, L, chi, N, a, b, initial_t, min_dt, max_dt);
}


} // namespace lowstorage

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMESTEPPER_HPP
