//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
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

#ifndef SUZERAIN_LOWSTORAGE_HPP
#define SUZERAIN_LOWSTORAGE_HPP

/** @file
 * Provides time integration-related operator details and advancement schemes.
 */

#include <suzerain/common.hpp>
#include <suzerain/math.hpp>
#include <suzerain/timecontroller.hpp>
#include <suzerain/timers.h>
#include <suzerain/traits.hpp>

namespace suzerain {

/**
 * Provides low-storage Runge-Kutta time integration schemes and the associated
 * operator interfaces.
 *
 * @see method_interface for details on this class of timestepping schemes.
 */
namespace lowstorage {

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
 *               = \frac{\eta_{i+1} - \eta_i}{\eta_{i+1}}
 * \f]
 * since the factor \f$\Delta{}t\f$ may be omitted without changing
 * \f$\bar{q}_{N-1}\f$ and \f$\eta_{0} = 0\f$.  Then maintaining the running
 * mean
 * \f[
 *   \bar{q}_{N-1} = \left(1 - \iota_{N-1}\right) \bar{q}_{N-2}
 *                 + \iota_{N-1} q_{N-1}
 * \f]
 * may be done via a simple update operation <tt>mean += iota_i * (sample -
 * mean)</tt>.  One use case for \f$\iota_i\f$ is obtaining running means of
 * quantities varying at each substep using minimal space and time overhead.
 *
 * Similarly to \f$\iota_i\f$, for quantities in effect only during linear
 * operator application, \f$\Delta{}t \alpha_i L\f$, one can define
 * \f[
 *   \iota_{\alpha,N-1} = \frac{\alpha_i}{\eta_{i+1}}.
 * \f]
 * Analogously,
 * \f[
 *   \iota_{\beta,N-1} = \frac{\beta_i}{\eta_{i+1}}
 * \f]
 * for is defined for quantities in effect only during linear operator
 * "inversion", \f$\left(M - \Delta{}t \beta_i L\right)^{-1}\f$.  The former is
 * useful for obtaining averages implicitly-computed forcing terms which are
 * applied for duration \f$\alpha_i \Delta{}t\f$ but computed during the
 * operator inversion stage.  In both cases, running weighted means are
 * computed in-place as <tt>mean += iota_alpha_i * (sample - mean)</tt> or
 * <tt>mean += iota_beta_i * (sample - mean)</tt>.
 *
 * @see linear_operator for the interface that \f$L\f$ must implement.
 * @see nonlinear_operator for the interface that \f$N\f$ must implement.
 * @see \ref smr91 and \ref yang11 for examples of essential information for a
 *      concrete scheme.
 * @see step() or substep() for methods that can advance state variables
 *      according to a timestepping method.
 */
template<typename Element>
class method_interface
{
public:

    /** The real-valued scalar corresponding to \c Element */
    typedef typename traits::component<Element>::type component;

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
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component alpha(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\beta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component beta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\gamma_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component gamma(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\zeta_i\f$ coefficient.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component zeta(std::size_t substep) const = 0;

    /**
     * Obtain the scheme's \f$\eta_i\f$ coefficient, which is
     * derived from other scheme coefficients.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps()]</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component eta(std::size_t substep) const = 0;

    /**
     * Compute the scheme's derived \f$\iota_i\f$ coefficient, which is used to
     * accumulate a running time-averaged value across substeps.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component iota(std::size_t substep) const = 0;

    /**
     * Compute the scheme's derived \f$\iota_{\alpha,i}\f$ coefficient, which is
     * used to accumulate a running time-averaged value across substeps.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component iota_alpha(std::size_t substep) const = 0;

    /**
     * Compute the scheme's derived \f$\iota_{\beta,i}\f$ coefficient, which is
     * used to accumulate a running time-averaged value across substeps.
     *
     * @param substep A substep number \f$i\f$ within <tt>[0,substeps())</tt>.
     *
     * @return The coefficient associated with the requested substep.
     */
    virtual component iota_beta(std::size_t substep) const = 0;

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
    virtual ~method_interface() {}
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
        const method_interface<Element>& m)
{
    return os << m.name();
}

/**
 * Defines the nonlinear operator interface required for timestepping.
 *
 * @see tparam State A state type which must provide an \c element \c typedef
 *                   containing its real- or complex-valued scalar type.
 */
template<typename State>
class nonlinear_operator
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename State::element element;

    /** The real-valued scalar corresponding to \c element */
    typedef typename traits::component<element>::type component;

    /**
     * Apply the operator in place on a given state location.  That is, compute
     * \f$\mbox{state}\leftarrow{}\mbox{operator}\left(\mbox{state}\right)\f$.
     *
     * @param time The time at which to apply the operator.
     * @param state The state location to use.  It is expected that certain
     *        implementations will require more specific state types.
     * @param method The low storage scheme being used.
     *        When <tt>substep_index == 0</tt> the operator should use this
     *        information to compute a stable time step.
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
    virtual std::vector<component> apply_operator(
            const component time,
            State& state,
            const method_interface<element>& method,
            const std::size_t substep_index) const = 0;

    /** Virtual destructor for peace of mind. */
    virtual ~nonlinear_operator() {}
};

/**
 * Compute an approximation to the maximum stable time step according to the
 * convective criterion for the three dimensional Euler equations.  This
 * surrogate is a model for the convective portion of the Navier-Stokes
 * operator.  This criterion appears as equation 2.39 in Wai Y. Kwok's thesis
 * (2002) and equations 4.20 and 4.21 in Stephen Guarini's thesis (1998):
 * \f[
 *     \Delta{}t \leq
 *     \frac{
 *       \left|\lambda_{I}\Delta{}t\right|_{\max}
 *     }{
 *         \left(\left|u_{x}\right| + a\right) \lambda^{(1)}_x
 *       + \left(\left|u_{y}\right| + a\right) \lambda^{(1)}_y
 *       + \left(\left|u_{z}\right| + a\right) \lambda^{(1)}_z
 *     }
 * \f]
 * where \f$a\f$ is the local acoustic velocity, \f$u_{x}\f$ denotes the
 * velocity in the X direction, \f$\lambda^{(1)}_x\f$ represents the maximum
 * imaginary eigenvalue magnitude of the first derivative operator in the $x$
 * direction, etc.  The maximum pure imaginary eigenvalue magnitude,
 * \f$\left|\lambda_{I}\Delta{}t\right|_{\mbox{max}}\f$, is a feature of the
 * chosen timestepping method.  For example, it is \f$\sqrt{3}\f$ for the \ref
 * smr91 scheme.  See the Suzerain model document for details on estimating
 * \f$\lambda^{(1)}_x\f$, \f$\lambda^{(1)}_y\f$, and \f$\lambda^{(1)}_z\f$.
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
 * calculation.  Three different acoustic velocities \f$a_x\f$, \f$a_y\f$ and
 * \f$a_z\f$ may be specified to accommodate when only some directions have
 * implicitly treated acoustics.
 *
 * @param u_x              Velocity in the X direction \f$u_{x}\f$
 * @param a_x              The local sound speed \f$a\f$ for the X direction.
 * @param lambda1_x        The maximum pure imaginary eigenvalue magnitude
 *                         for the first derivative in the X direction.
 * @param u_y              Velocity in the Y direction \f$u_{y}\f$
 * @param a_y              The local sound speed \f$a\f$ for the Y direction.
 * @param lambda1_y        The maximum pure imaginary eigenvalue magnitude
 *                         for the first derivative in the Y direction.
 * @param u_z              Velocity in the Z direction \f$u_{z}\f$
 * @param a_z              The local sound speed \f$a\f$ for the Z direction.
 * @param lambda1_z        The maximum pure imaginary eigenvalue magnitude
 *                         for the first derivative in the Z direction.
 * @param evmaxmag_imag    The maximum pure imaginary eigenvalue magnitude
 *                         for some Runge-Kutta scheme, denoted \f$\left|
 *                         \lambda_{I}\Delta_{}t \right|_{\mbox{max}}\f$
 *                         in Guarini's thesis.
 *
 * @return The maximum stable time step \f$\Delta{}t\f$ according to
 *         the convective CFL criterion.
 */
template<typename FPT>
FPT convective_stability_criterion(
        const FPT u_x,
        const FPT a_x,
        const FPT lambda1_x,
        const FPT u_y,
        const FPT a_y,
        const FPT lambda1_y,
        const FPT u_z,
        const FPT a_z,
        const FPT lambda1_z,
        const FPT evmaxmag_imag)
{
    using ::std::abs;
    return  evmaxmag_imag
        /  (   (abs(u_x) + a_x)*lambda1_x
             + (abs(u_y) + a_y)*lambda1_y
             + (abs(u_z) + a_z)*lambda1_z);
}

/**
 * Compute an approximation to the maximum stable time step according to the
 * convective criterion for the one dimensional Euler equations.  This
 * surrogate is a model for the convective portion of the one-dimensional
 * Navier-Stokes operator:
 * \f[
 *     \Delta{}t \leq
 *     \frac{
 *       \left|\lambda_{I}\Delta{}t\right|_{\max}
 *     }{
 *         \left(\left|u\right| + a\right) \lambda^{(1)}
 *     }
 * \f]
 * where \f$a\f$ is the local acoustic velocity, \f$u\f$ denotes the velocity,
 * and represents the maximum imaginary eigenvalue magnitude of the first
 * derivative operator.
 *
 * @param u                Velocity \f$u\f$
 * @param a                The local sound speed \f$a\f$.
 * @param lambda1          The maximum pure imaginary eigenvalue magnitude
 *                         for the first derivative.
 * @param evmaxmag_imag    The maximum pure imaginary eigenvalue magnitude
 *                         for some Runge-Kutta scheme, denoted \f$\left|
 *                         \lambda_{I}\Delta_{}t \right|_{\mbox{max}}\f$
 *                         in Guarini's thesis.
 *
 * @see convective_stability_criterion(FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT)
 *      for more details.
 * @return The maximum stable time step \f$\Delta{}t\f$ according to
 *         the convective CFL criterion.
 */
template<typename FPT>
FPT convective_stability_criterion(
        const FPT u,
        const FPT a,
        const FPT lambda1,
        const FPT evmaxmag_imag)
{
    using ::std::abs;
    return evmaxmag_imag / ((abs(u) + a)*lambda1);
}

/**
 * Compute an approximation to the maximum stable time step according to the
 * viscous stability criterion for a three-dimensional model diffusion
 * equation.  This surrogate is a model for the diffusive part of the
 * Navier-Stokes operator.  This criterion appears as equation 2.40 in Wai Y.
 * Kwok's thesis (2002) and equations 4.29 and 4.30 in Stephen Guarini's thesis
 * (1998):
 * \f[
 *  \Delta{}t \leq
 *  \frac{
 *      \left|\lambda_{R}\Delta_{}t\right|_{\max}
 *  }{
 *    \max\left(
 *      \left|\frac{\gamma\left(\nu-\nu_{0}\right)}{\mbox{Re}\mbox{Pr}}\right|,
 *      \left|\frac{\nu-\nu_{0}}{\mbox{Re}}\right|,
 *      \left|\frac{\nu_{B}-\nu_{B0}}{\mbox{Re}}\right|
 *    \right)
 *    \left(
 *        \lambda^{(2)}_x
 *      + \lambda^{(2)}_y
 *      + \lambda^{(2)}_z
 *    \right)
 *  }
 * \f]
 * where a bulk kinematic viscosity \f$\nu_{B}\f$ has been added.  See the
 * Suzerain model document for details on estimating \f$\lambda^{(2)}_x\f$,
 * \f$\lambda^{(2)}_y\f$, and \f$\lambda^{(2)}_z\f$.  The maximum pure real
 * eigenvalue magnitude, \f$\left|\lambda_{R}\Delta{}t\right|_{\mbox{max}}\f$,
 * is a feature of the chosen timestepping method.  For example, it is 2.512
 * for the \ref smr91 scheme.  The absolute values within the maximum operation
 * account for the possibility that \f$\nu<\nu_{0}\f$.
 *
 * Using a hybrid implicit/explicit %timestepper with viscous terms computed
 * implicitly sets \f$\nu_0\f$ to be the reference kinematic viscosity about
 * which the viscous terms were linearized.
 *
 * @param lambda2_x        The maximum pure real eigenvalue magnitude
 *                         for the second derivative in the X direction.
 * @param lambda2_y        The maximum pure real eigenvalue magnitude
 *                         for the second derivative in the Y direction.
 * @param lambda2_z        The maximum pure real eigenvalue magnitude
 *                         for the second derivative in the Z direction.
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
        const FPT lambda2_x,
        const FPT lambda2_y,
        const FPT lambda2_z,
        const FPT Re,
        const FPT Pr,
        const FPT gamma,
        const FPT evmaxmag_real,
        const FPT nu,
        const FPT nu0,
        const FPT nuB,
        const FPT nuB0)
{
    // Find maximum diffusive coefficient amongst all possible criteria
    using std::abs;
    using std::max;
    const FPT maxcoeff = max(max(gamma/Pr,FPT(1))*abs(nu  - nu0 )/Re,
                                                  abs(nuB - nuB0)/Re);

    return  evmaxmag_real
        /  (maxcoeff * (lambda2_x + lambda2_y + lambda2_z));
}

/**
 * Compute an approximation to the maximum stable time step according to the
 * viscous stability criterion for a one-dimensional model diffusion
 * equation.  This surrogate is a model for the diffusive part of the
 * Navier-Stokes operator:
 * \f[
 *  \Delta{}t \leq
 *  \frac{
 *      \left|\lambda_{R}\Delta_{}t\right|_{\max}
 *  }{
 *    \max\left(
 *      \left|\frac{\gamma\left(\nu-\nu_{0}\right)}{\mbox{Re}\mbox{Pr}}\right|,
 *      \left|\frac{\nu-\nu_{0}}{\mbox{Re}}\right|,
 *      \left|\frac{\nu_{B}-\nu_{B0}}{\mbox{Re}}\right|
 *    \right)
 *    \left(
 *        \lambda^{(2)}
 *    \right)
 *  }
 * \f]
 * where a bulk kinematic viscosity \f$\nu_{B}\f$ has been added.
 *
 * @param lambda2          The maximum pure real eigenvalue magnitude
 *                         for the second derivative.
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
 * @see diffusive_stability_criterion(FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT,FPT)
 *      for more details.
 * @return The maximum stable time step \f$\Delta{}t\f$ according to
 *         the diffusive criterion.
 */
template<typename FPT>
FPT diffusive_stability_criterion(
        const FPT lambda2,
        const FPT Re,
        const FPT Pr,
        const FPT gamma,
        const FPT evmaxmag_real,
        const FPT nu,
        const FPT nu0,
        const FPT nuB,
        const FPT nuB0)
{

    // Find maximum diffusive coefficient amongst all possible criteria
    using std::abs;
    using std::max;
    const FPT maxcoeff = max(max(gamma/Pr,FPT(1))*abs(nu  - nu0 )/Re,
                                                  abs(nuB - nuB0)/Re);

    return  evmaxmag_real / (maxcoeff * lambda2);
}

namespace { // anonymous

/** A comparator which considers NaNs to be less than all other values */
struct NaN_is_minimum_comparator {
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
struct delta_t_reducer {
    template< typename T>
    T operator()(const std::vector<T> & candidates) {
        assert(candidates.size() > 0);
        return *std::min_element(candidates.begin(), candidates.end(),
                                 NaN_is_minimum_comparator());
    }
};

/**
 * Defines the linear operator interface required for low storage timestepping.
 *
 * @see tparam StateA A state type which must provide an \c element \c typedef
 *                    containing its real- or complex-valued scalar type.
 * @see tparam StateB A state type which must provide an \c element \c typedef
 *                    containing its real- or complex-valued scalar type.
 */
template<typename StateA, typename StateB = StateA>
class linear_operator
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename StateA::element element;
    BOOST_STATIC_ASSERT(
            (boost::is_same<element,typename StateB::element>::value));

    /** The real-valued scalar corresponding to \c element */
    typedef typename traits::component<element>::type component;

    /**
     * Apply \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.  That is,
     * \f$\mbox{state}\leftarrow{} \left(M+\phi{}L\right)\mbox{state}\f$.
     *
     * @param phi           Scale factor \f$\phi\f$ to use.
     * @param state         State vector on which to apply the operator.
     * @param substep_index The (zero-indexed) time stepper substep index.
     */
    virtual void apply_mass_plus_scaled_operator(
            const element& phi,
            StateA& state,
            const std::size_t substep_index) const = 0;

    /**
     * Accumulate \f$M+\phi{}L\f$ out-of-place for scalar \f$\phi\f$.
     * That is, \f$\mbox{output}\leftarrow{}
     * \left(M+\phi{}L\right) \mbox{input}+\beta\mbox{output}\f$.
     *
     * @param phi           Scale factor \f$\phi\f$ to use.
     * @param input         State vector on which to apply the operator.
     * @param beta          Scale factor for output vector during accumulation.
     * @param output        State vector into which to accumulate the result.
     * @param substep_index The (zero-indexed) time stepper substep index.
     */
    virtual void accumulate_mass_plus_scaled_operator(
            const element& phi,
            const StateA& input,
            const element& beta,
            StateB& output,
            const std::size_t substep_index) const = 0;

    /**
     * Invert \f$M+\phi{}L\f$ in-place for scalar \f$\phi\f$.  That is,
     * \f$\mbox{state}\leftarrow{}\left(M+\phi{}L\right)^{-1}\mbox{state}\f$.
     *
     * Low storage timestepping schemes do not make make use of \c ic0.  It is
     * intended as a problem-specific hook for implementing integral
     * constraints.  When \c ic0 is non-NULL, some portion of the inverse
     * operator also must be applied to some portion of the data in \c ic0.
     * That is, \f$\mbox{ic0} \leftarrow{} \left(M+\phi{}L\right)^{-1}
     * \mbox{ic0}\f$.  In a mixed Fourier discretiation, a mean integral
     * constraint is applied via applying the inverse ``zero-zero'' operator
     * and hence the name.  The exact behavior is subclass dependent.
     *
     * @param phi           Scale factor \f$\phi\f$ to use.
     * @param state         State vector on which to apply the operator.
     * @param method        The low storage scheme being used
     * @param delta_t       The size of the currently active time step.
     * @param substep_index The (zero-indexed) time stepper substep index.
     * @param ic0           Additional state, often used for imposing
     *                      integral constraints, which is ignored when NULL.
     */
    virtual void invert_mass_plus_scaled_operator(
            const element& phi,
            StateA& state,
            const method_interface<element>& method,
            const component delta_t,
            const std::size_t substep_index,
            StateA *ic0 = NULL) const = 0;

    /** Virtual destructor for peace of mind. */
    virtual ~linear_operator() {}
};

/**
 * Implements simple multiplicative operator which scales all state
 * variables by a uniform factor.  The associated mass matrix \f$M\f$
 * is the identity matrix.
 *
 * @tparam StateA A type which likely descends from state_base.
 * @tparam StateB A type which likely descends from state_base.
 */
template< typename StateA, typename StateB = StateA >
class multiplicative_operator
    : public linear_operator<StateA,StateB>,
      public nonlinear_operator<StateB>
{
public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename StateA::element element;
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename StateB::element>::value));

    /** The real-valued scalar corresponding to \c element */
    typedef typename traits::component<element>::type component;

    /**
     * Construct an instance which scales by \c factor and reports
     * \c delta_t as a stable time step.
     *
     * @param factor uniform scaling factor to apply.
     * @param delta_t uniform, presumably stable time step to
     *        always return from apply_operator().
     */
    template< typename FactorType, typename DeltaTType >
    multiplicative_operator(const FactorType& factor, const DeltaTType& delta_t)
        : factor(factor)
        , delta_t(delta_t)
    {}

    /**
     * Construct an instance which scales by \c factor and reports the
     * maximum representable floating point value as a stable time step.
     * Useful in testing contexts or during operator composition.
     *
     * @param factor uniform scaling factor to apply.
     */
    template< typename FactorType >
    multiplicative_operator(const FactorType& factor)
        : factor(factor)
        , delta_t(std::numeric_limits<component>::infinity())
    {}

    /**
     * Scale \c state by the factor set at construction time.
     *
     * @param time The time at which to apply the operator (ignored).
     * @param state to scale in place.
     * @param method Ignored in this implementation.
     * @param substep_index Ignored in this implementation.
     *
     * @return The \c delta_t provided at construction time.
     */
    virtual std::vector<component> apply_operator(
            const component time,
            StateB& state,
            const method_interface<element>& method,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(method);
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
     * @param substep_index Ignored.
     */
    virtual void apply_mass_plus_scaled_operator(
            const element& phi,
            StateA& state,
            const std::size_t substep_index = 0) const
    {
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
     * @param substep_index Ignored.
     */
    virtual void accumulate_mass_plus_scaled_operator(
            const element& phi,
            const StateA& input,
            const element& beta,
            StateB& output,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(substep_index);
        output.scale(beta);
        output.add_scaled(phi*factor + element(1), input);
    }

    /**
     * Compute \f$\mbox{state}\leftarrow{}
     * \left(I+\phi\times\mbox{factor}\right)^{-1} \mbox{state}\f$ where \c
     * factor is the scaling factor set at construction time.
     *
     * @param phi           Additional scaling \f$\phi\f$ to apply.
     * @param state         State to modify in place.
     * @param delta_t       Ignored.
     * @param substep_index Ignored.
     * @param ic0           When non-null, modified as \c state.
     */
    virtual void invert_mass_plus_scaled_operator(
            const element& phi,
            StateA& state,
            const method_interface<element>& method,
            const component delta_t = 0,
            const std::size_t substep_index = 0,
            StateA *ic0 = NULL) const
    {
        SUZERAIN_UNUSED(method);
        SUZERAIN_UNUSED(delta_t);
        SUZERAIN_UNUSED(substep_index);
        state.scale((element(1))/(phi*factor + element(1)));

        if (ic0) ic0->scale((element(1))/(phi*factor + element(1)));
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
 * @see multiplicative_operator<NumDims,Element,Storage,CompatibleStorage> for
 *      more details.
 */
template<
    typename FactorType,
    typename DeltaTType,
    typename StateA,
    typename StateB
>
multiplicative_operator<StateA,StateB>
make_multiplicative_operator(
    const FactorType& factor,
    const DeltaTType& delta_t,
    const StateA& input,
    const StateB& output)
{
    multiplicative_operator<StateA,StateB> retval(factor, delta_t);
    return retval;
}

/**
 * Compute redundant information for an method_interface given essential,
 * scheme-specific constants.
 *
 * The \c Scheme parameter must be a templated, static-only struct
 * encapsulating a small number of independent pieces of information:
 * <dl>
 * <dt>name</dt>
 * <dd>A human-readable name for the scheme</dd>
 * <dt>substeps</dt>
 * <dd>Number of substeps within the scheme</dd>
 * <dt>evmaxmag_real</dt>
 * <dd>The maximum purely real stable eigenvalue magnitude</dd>
 * <dt>evmaxmag_imag</dt>
 * <dd>The maximum purely imaginary stable eigenvalue magnitude</dd>
 * <dt>alpha_numerator</dt>
 * <dd>The numerators for the \f$\alpha_i\f$ written using \c denominator</dd>
 * <dt>beta_numerator</dt><dd></dd>
 * <dd>The numerators for the \f$\beta_i\f$ written using \c denominator</dd>
 * <dt>gamma_numerator</dt><dd></dd>
 * <dd>The numerators for the \f$\gamma_i\f$ written using \c denominator</dd>
 * <dt>denominator</dt><dd></dd>
 * <dd>The least common multiple of all denominators appearing within
 *     \f$\alpha\f$, \f$\beta\f$, or \f$\gamma\f$.
 * </dl>
 * Numerators and denominators are tracked separately to permit performing
 * integer operations wherever possible.  This permits maintaining full
 * precision constants via very simple rational operations.  All computed
 * results may be accessed as if contained in static arrays.  For example,
 * <code>eta[1]</code> or <code>iota_beta[2]</code>.
 *
 * @tparam Scheme A class encapsulating the minimum constants
 *         required to describe a low storage scheme.
 * @tparam Component Real-valued type for which constants are returned.
 * @tparam Integer Signed integer type used for indexing and computation.
 *
 * @see method_interface for a full definition of the constants involved.
 * @see \ref smr91 for an example implementation.
 */
template <template <typename,typename> class Scheme,
          typename Component,
          typename Integer = std::ptrdiff_t>
class constants : private Scheme<Component,Integer>
{

private:

    /** Encapsulates only the constants necessary to define the scheme. */
    typedef Scheme<Component,Integer> scheme;

public:

    using scheme::name;
    using scheme::substeps;
    using scheme::evmaxmag_real;
    using scheme::evmaxmag_imag;

private:

    /** Helper for implementing #alpha */
    struct alpha_type {

        /** Computes \f$\alpha_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(scheme::alpha_numerator[i]) / scheme::denominator;
        }

    };

    /** Helper for implementing #beta */
    struct beta_type {

        /** Computes \f$\beta_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(scheme::beta_numerator[i]) / scheme::denominator;
        }

    };

    /** Helper for implementing #gamma */
    struct gamma_type {

        /** Computes \f$\gamma_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(scheme::gamma_numerator[i]) / scheme::denominator;
        }

    };

    /** Helper for implementing #zeta */
    struct zeta_type {

        /** Computes numerator of \f$\zeta_i\f$ given \c Scheme */
        Integer numerator(const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return scheme::alpha_numerator[i]
                 + scheme::beta_numerator [i]
                 - scheme::gamma_numerator[i];
        }

        /** Computes \f$\zeta_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(numerator(i)) / scheme::denominator;
        }

    };

    /** Helper for implementing #eta */
    struct eta_type {

        /** Computes numerator of \f$\eta_i\f$ given \c Scheme */
        Integer numerator(const Integer i) const
        {
            assert(0 <= i && i <= substeps); // i == substeps OK
            Integer v = 0;
            for (Integer j = i; j --> 0 ;) {
                v += scheme::alpha_numerator[j];
                v += scheme::beta_numerator [j];
            }
            return v;
        }

        /** Computes \f$\eta_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i <= substeps); // i == substeps OK
            return Component(numerator(i)) / scheme::denominator;
        }

    };

    /** Helper for implementing #iota */
    struct iota_type {

        /** Computes \f$\iota_i\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(eta.numerator(i+1) - eta.numerator(i))
                 / eta.numerator(i+1);
        }

    };

    /** Helper for implementing #iota_alpha */
    struct iota_alpha_type {

        /** Computes \f$\iota_{\alpha,i}\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(scheme::alpha_numerator[i]) / eta.numerator(i+1);
        }

    };

    /** Helper for implementing #iota_beta */
    struct iota_beta_type {

        /** Computes \f$\iota_{\beta,i}\f$ given \c Scheme */
        Component operator[](const Integer i) const
        {
            assert(0 <= i && i < substeps);
            return Component(scheme::beta_numerator[i]) / eta.numerator(i+1);
        }

    };

public:

    /**
     * Computes \f$\alpha_i\f$ given \c Scheme as <tt>alpha[i]</tt>.
     * @see method_interface::alpha
     */
    static const alpha_type alpha;

    /**
     * Computes \f$\beta_i\f$ given \c Scheme as <tt>beta[i]</tt>.
     * @see method_interface::beta
     */
    static const beta_type beta;

    /**
     * Computes \f$\gamma_i\f$ given \c Scheme as <tt>gamma[i]</tt>.
     * @see method_interface::gamma
     */
    static const gamma_type gamma;

    /**
     * Computes \f$\zeta_i\f$ given \c Scheme as <tt>zeta[i]</tt>.
     * @see method_interface::zeta
     */
    static const zeta_type zeta;

    /**
     * Computes \f$\eta_i\f$ given \c Scheme as <tt>eta[i]</tt>.
     * @see method_interface::eta
     */
    static const eta_type eta;

    /**
     * Computes \f$\iota_i\f$ given \c Scheme as <tt>iota[i]</tt>.
     * @see method_interface::iota
     */
    static const iota_type iota;

    /**
     * Computes \f$\iota_{\alpha,i}\f$ given
     * \c Scheme as <tt>iota_alpha[i]</tt>.
     * @see method_interface::iota_alpha
     */
    static const iota_alpha_type iota_alpha;

    /**
     * Computes \f$\iota_{\beta,i}\f$ given
     * \c Scheme as <tt>iota_beta[i]</tt>.
     * @see method_interface::iota_beta
     */
    static const iota_beta_type iota_beta;

};

// *************************************************************
// BEGIN hideousness for static constant structs within Contants
// *************************************************************

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::alpha_type
constants<Scheme,Component,Integer>::alpha = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::beta_type
constants<Scheme,Component,Integer>::beta = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::gamma_type
constants<Scheme,Component,Integer>::gamma = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::zeta_type
constants<Scheme,Component,Integer>::zeta = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::eta_type
constants<Scheme,Component,Integer>::eta = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::iota_type
constants<Scheme,Component,Integer>::iota = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::iota_alpha_type
constants<Scheme,Component,Integer>::iota_alpha = {};

template <template <typename,typename> class Scheme,
          typename Component, typename Integer>
const typename constants<Scheme,Component,Integer>::iota_beta_type
constants<Scheme,Component,Integer>::iota_beta = {};

// ***********************************************************
// END hideousness for static constant structs within Contants
// ***********************************************************

/**
 * Given a \c Scheme and \c Element type, encapsulates a low storage
 * method behind the method_interface interface using constants.
 *
 * @tparam Scheme An essential set of scheme-specific constants usable
 *         within the constants class.
 * @tparam A real- or complex-valued scalar type to be used.
 *         constants are returned as the corresponding real type,
 *         called \c component.
 * @see \ref smr91 or \ref yang11 for examples of valid Schemes to supply.
 */
template <template <typename,typename> class Scheme, typename Element>
class method : public method_interface<Element>
{

public:

    /** The real-valued scalar corresponding to \c Element */
    typedef typename method_interface<Element>::component component;

    /** The concrete scheme type providing method-specific constants  */
    typedef constants<Scheme,component> scheme;

    /**
     * Explicit constructor.
     *
     * @param evmagfactor The multiplicative factor to use when reporting
     *                    maximum pure real and pure imaginary eigenvalue
     *                    magnitudes in evmaxmag_real() and evmaxmag_imag(),
     *                    respectively.
     */
    explicit method(component evmagfactor = 1)
        : evmaxmag_real_(evmagfactor * scheme::evmaxmag_real()),
          evmaxmag_imag_(evmagfactor * scheme::evmaxmag_imag())
        { assert(evmagfactor > 0); }

    /** @copydoc method_interface::name */
    virtual const char * name() const
    { return scheme::name; }

    /** @copydoc method_interface::substeps */
    virtual std::size_t substeps() const
    { return scheme::substeps; }

    /** @copydoc method_interface::alpha */
    virtual component alpha(const std::size_t substep) const
    { return scheme::alpha[substep]; }

    /** @copydoc method_interface::beta */
    virtual component beta(const std::size_t substep) const
    { return scheme::beta[substep]; }

    /** @copydoc method_interface::gamma */
    virtual component gamma(const std::size_t substep) const
    { return scheme::gamma[substep]; }

    /** @copydoc method_interface::zeta */
    virtual component zeta(const std::size_t substep) const
    { return scheme::zeta[substep]; }

    /** @copydoc method_interface::eta */
    virtual component eta(const std::size_t substep) const
    { return scheme::eta[substep]; }

    /** @copydoc method_interface::iota */
    virtual component iota(const std::size_t substep) const
    { return scheme::iota[substep]; }

    /** @copydoc method_interface::iota_alpha */
    virtual component iota_alpha(const std::size_t substep) const
    { return scheme::iota_alpha[substep]; }

    /** @copydoc method_interface::iota_beta */
    virtual component iota_beta(const std::size_t substep) const
    { return scheme::iota_beta[substep]; }

    /** @copydoc method_interface::evmaxmag_real */
    virtual component evmaxmag_real() const
    { return evmaxmag_real_; }

    /** @copydoc method_interface::evmaxmag_imag */
    virtual component evmaxmag_imag() const
    { return evmaxmag_imag_; }

private:

    /** Value to report from evmaxmag_real(). */
    component evmaxmag_real_;

    /** Value to report from evmaxmag_imag(). */
    component evmaxmag_imag_;
};

/**
 * Encapsulates essential constants for the three stage, third order scheme
 * from Appendix A of Spalart, Moser, and Rogers' 1991 ``Spectral Methods for
 * the Navier-Stokes Equations with One Infinite and Two Periodic Directions''
 * published in the <em>Journal of Computational Physics</em> volume 96 pages
 * 297-324.
 *
 * @see Designed to be used with the constants template.
 */
template <typename Component, typename Integer>
struct smr91
{
    /** A human-readable name for the scheme */
    static const char *name;

    /** Number of substeps within the scheme */
    static const Integer substeps = 3;

    /**
     * The maximum purely real eigenvalue magnitude for the scheme, denoted
     * \f$\left| \lambda_{R}\Delta_{}t \right|_{\mbox{max}}\f$ in Guarini's
     * thesis.
     */
    static Component evmaxmag_real()
    { return Component(2.51274532661832862402373L); }

    /**
     * The maximum purely imaginary eigenvalue magnitude for the scheme, denoted
     * \f$\left| \lambda_{I}\Delta_{}t \right|_{\mbox{max}}\f$ in Guarini's
     * thesis.
     */
    static Component evmaxmag_imag()
    { using std::sqrt; return sqrt(3); }

    /** The numerators for the \f$\alpha_i\f$ written using \c denominator */
    static const Integer alpha_numerator[substeps];

    /** The numerators for the \f$\beta_i\f$ written using \c denominator */
    static const Integer beta_numerator [substeps];

    /** The numerators for the \f$\gamma_i\f$ written using \c denominator */
    static const Integer gamma_numerator[substeps];

    /**
     * The least common multiple of all denominators appearing within
     * \f$\alpha\f$, \f$\beta\f$, or \f$\gamma\f$.
     */
    static const Integer denominator = 480;
};

template <typename Component, typename Integer>
const char * smr91<Component,Integer>::name = "smr91";

template <typename Component, typename Integer>
const Integer smr91<Component,Integer>::alpha_numerator[substeps] = {
         29 * (denominator / 96),
        - 3 * (denominator / 40),
          1 * (denominator /  6)
};

template <typename Component, typename Integer>
const Integer smr91<Component,Integer>::beta_numerator[substeps] = {
         37 * (denominator / 160),
          5 * (denominator /  24),
          1 * (denominator /   6)
};

template <typename Component, typename Integer>
const Integer smr91<Component,Integer>::gamma_numerator[substeps] = {
          8 * (denominator /  15),
          5 * (denominator /  12),
          3 * (denominator /   4)
};

/**
 * Encapsulates the three stage, second-order, adjoint-consistent scheme from
 * Shan Yang's 2011 thesis ``A shape Hessian based analysis of roughness
 * effects on fluid flows''.
 *
 * @see Designed to be used with the constants template.
 */
template <typename Component, typename Integer>
struct yang11
{
    /** A human-readable name for the scheme */
    static const char *name;

    /** Number of substeps within the scheme */
    static const Integer substeps = 3;

    /** The maximum purely real eigenvalue magnitude for the scheme. */
    static Component evmaxmag_real()
    { return Component(2.51274532661832862402373L); }

    /** The maximum purely imaginary eigenvalue magnitude for the scheme. */
    static Component evmaxmag_imag()
    { using std::sqrt; return sqrt(3); }

    /** The numerators for the \f$\alpha_i\f$ written using \c denominator */
    static const Integer alpha_numerator[substeps];

    /** The numerators for the \f$\beta_i\f$ written using \c denominator */
    static const Integer beta_numerator [substeps];

    /** The numerators for the \f$\gamma_i\f$ written using \c denominator */
    static const Integer gamma_numerator[substeps];

    /**
     * The least common multiple of all denominators appearing within
     * \f$\alpha\f$, \f$\beta\f$, or \f$\gamma\f$.
     */
    static const Integer denominator = 6;
};

template <typename Component, typename Integer>
const char * yang11<Component,Integer>::name = "Yang11";

template <typename Component, typename Integer>
const Integer yang11<Component,Integer>::alpha_numerator[substeps] = {
         1 * (denominator / 3),
        -1 * (denominator / 2),
         1 * (denominator / 3)
};

template <typename Component, typename Integer>
const Integer yang11<Component,Integer>::beta_numerator[substeps] = {
         1 * (denominator / 6),
         2 * (denominator / 3),
         0 * (denominator / 1)
};

template <typename Component, typename Integer>
const Integer yang11<Component,Integer>::gamma_numerator[substeps] = {
         1 * (denominator / 2),
         1 * (denominator / 3),
         1 * (denominator / 1)
};

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme to advance the system \f$ M u_t = Lu +
 * \chi N(u,t) \f$.  Notice the presence of scaling factor \f$\chi\f$ within
 * the substep equation. Because the roles of state locations \c a and \c b are
 * swapped after each substep, and in contrast to step(), only a single state
 * type may be supplied to this method.
 *
 * @param m             The low storage scheme to use.
 *                      For example, \ref method in conjunction with \ref smr91.
 * @param L             The linear operator to be treated implicitly.
 * @param chi           The factor \f$\chi\f$ scaling the nonlinear operator.
 * @param N             The nonlinear operator to be treated explicitly.
 * @param time          The simulation time \f$t\f$ at which the substep occurs.
 * @param a             On entry contains \f$u^{i}\f$ and on exit contains
 *                      \f$N\left(u^{i}\right)\f$.
 * @param b             On entry contains \f$N\left(u^{i-1}\right)\f$ and
 *                      on exit contains \f$u^{i+1}\f$.
 * @param substep_index The substep number to take.
 * @param delta_t       The time step \f$\Delta{}t\f$ to take.  The same time
 *                      step should be supplied for all substep computations
 *                      throughout a given step.
 * @return The time step \f$\Delta{}t\f$ taken. It will always equal \c delta_t.
 *
 * @see method_interface for the equation governing time advancement.
 * @see The method step() provides more convenient ways to perform multiple
 *      substeps, including dynamic step size computation.
 */
template< typename Element,
          typename LinearState,
          typename NonlinearState,
          typename State >
const typename traits::component<Element>::type substep(
    const method_interface<Element>& m,
    const linear_operator<LinearState>& L,
    const typename traits::component<Element>::type chi,
    const nonlinear_operator<NonlinearState>& N,
    const typename traits::component<Element>::type time,
    State& a,
    State& b,
    const typename traits::component<Element>::type delta_t,
    const std::size_t substep_index)
{
    BOOST_STATIC_ASSERT(
            (boost::is_same<Element,typename    LinearState::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<Element,typename NonlinearState::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<Element,typename          State::element>::value));

    if (SUZERAIN_UNLIKELY(substep_index >= m.substeps()))
        throw std::invalid_argument("Requested substep too large");

    L.accumulate_mass_plus_scaled_operator(
                  delta_t * m.alpha(substep_index), a,
            chi * delta_t * m.zeta(substep_index),  b,
            substep_index);
    N.apply_operator(time + delta_t * m.eta(substep_index), a,
                     m, substep_index);
    b.add_scaled(chi * delta_t * m.gamma(substep_index), a);
    L.invert_mass_plus_scaled_operator(-delta_t * m.beta(substep_index), b,
                                       m, delta_t, substep_index);

    return delta_t;
}

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu + \chi
 * N(u,t)\f$.  Notice the presence of scaling factor \f$\chi\f$ within the
 * substep equation.  The time step taken, \f$\Delta{}t\f$, will be based on a
 * stable value computed during the first nonlinear operator application as
 * well as an optional fixed maximum step size.
 *
 * @param m       The low storage scheme to use.
 *                For example, \ref method in conjunction with \ref smr91.
 * @param reducer A stateful functor taking a vector of stable time step
 *                candidates down to a single stable time step.  Users may
 *                employ a custom functor compatible with delta_t_reducer to
 *                add logic for monitoring or manipulating stable step criteria.
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
 * @see method_interface for the equation governing time advancement.
 */
template< typename Element, typename Reducer,
          typename LinearA, typename LinearB, typename NonlinearB,
          typename StateA, typename StateB >
const typename traits::component<Element>::type step(
    const method_interface<Element>& m,
    Reducer& reducer,
    const linear_operator<LinearA,LinearB>& L,
    const typename traits::component<Element>::type chi,
    const nonlinear_operator<NonlinearB>& N,
    const typename traits::component<Element>::type time,
    StateA& a,
    StateB& b,
    const typename traits::component<Element>::type max_delta_t = 0)
{
    SUZERAIN_TIMER_SCOPED("lowstorage::step");

    using boost::is_same;
    BOOST_STATIC_ASSERT((is_same<Element,typename    LinearA::element>::value));
    BOOST_STATIC_ASSERT((is_same<Element,typename    LinearB::element>::value));
    BOOST_STATIC_ASSERT((is_same<Element,typename NonlinearB::element>::value));
    BOOST_STATIC_ASSERT((is_same<Element,typename     StateA::element>::value));
    BOOST_STATIC_ASSERT((is_same<Element,typename     StateB::element>::value));
    typedef typename traits::component<Element>::type component_type;

    // First substep handling is special since we need to determine delta_t
    b.assign_from(a);
    const std::vector<component_type> delta_t_candidates
        = N.apply_operator(time, b, m, 0);
    component_type delta_t = reducer(delta_t_candidates);

    if (max_delta_t > 0) {
        delta_t = math::minnan(delta_t, max_delta_t);
    }
    L.apply_mass_plus_scaled_operator(delta_t * m.alpha(0), a, 0);
    a.add_scaled(chi * delta_t * m.gamma(0), b);
    L.invert_mass_plus_scaled_operator(-delta_t * m.beta(0), a, m, delta_t, 0);

    // Second and subsequent substeps are identical
    for (std::size_t i = 1; i < m.substeps(); ++i) {
        L.accumulate_mass_plus_scaled_operator(
                delta_t * m.alpha(i),      a,
                chi * delta_t * m.zeta(i), b,
                i);
        b.exchange(a); // Note nonlinear storage controls exchange operation
        N.apply_operator(time + delta_t * m.eta(i), b,
                         m, i);
        a.add_scaled(chi * delta_t * m.gamma(i), b);
        L.invert_mass_plus_scaled_operator(-delta_t * m.beta(i),
                                           a, m, delta_t, i);
    }

    return delta_t;
}

/**
 * Using the given method and a linear and nonlinear operator, take substep \c
 * substep_index while advancing from \f$u(t)\f$ to \f$u(t+\Delta{}t)\f$ using
 * a hybrid implicit/explicit scheme using the system \f$ M u_t = Lu + \chi
 * N(u,t)\f$.  Notice the presence of the factor \f$chi\f$ within the substep
 * equation.  The time step taken, \f$\Delta{}t\f$, will be based on a stable
 * value computed during the first nonlinear operator application as well as an
 * optional fixed maximum step size.
 *
 * @param m       The low storage scheme to use.
 *                For example, \ref method in conjunction with \ref smr91.
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
 * @see method_interface for the equation governing time advancement.
 */
template< typename Element,
          typename LinearA, typename LinearB, typename NonlinearB,
          typename StateA, typename StateB >
const typename traits::component<Element>::type step(
    const method_interface<Element>& m,
    const linear_operator<LinearA,LinearB>& L,
    const typename traits::component<Element>::type chi,
    const nonlinear_operator<NonlinearB>& N,
    const typename traits::component<Element>::type time,
    StateA& a,
    StateB& b,
    const typename traits::component<Element>::type max_delta_t = 0)
{
    delta_t_reducer reducer;
    return step<
            Element, delta_t_reducer,
            LinearA, LinearB, NonlinearB,
            StateA, StateB
        >(m, reducer, L, chi, N, time, a, b, max_delta_t);
}

/**
 * Provides higher-level control mechanisms built atop low storage time
 * integration schemes.  This class is a thin wrapper combining
 * \ref timecontroller with step().
 *
 * @see \ref timecontroller for details on the time controller logic.
 * @see make_controller for an easy way to create
 *      an instance with the appropriate type signature.
 */
template< typename StateA,
          typename StateB,
          typename Reducer,
          typename LinearA    = StateA,
          typename LinearB    = StateB,
          typename NonlinearB = StateB >
class controller
    : public timecontroller< typename traits::component<
            typename StateA::element
      >::type >
{
protected:

    /** Shorthand for the superclass */
    typedef timecontroller< typename traits::component<
                typename StateA::element
            >::type > super;

public:

    /** The real- or complex-valued scalar type the operator understands */
    typedef typename StateA::element element;
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename    LinearA::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename    LinearB::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename NonlinearB::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename     StateA::element>::value));
    BOOST_STATIC_ASSERT(
            (boost::is_same<element, typename     StateB::element>::value));

    /**
     * Construct an instance that will advance a simulation built atop the
     * given operators and storage.
     *
     * @param m         The low storage scheme to use.
     *                  For example, \ref method in conjunction with \ref smr91.
     * @param reducer   A stateful functor taking a vector of stable time step
     *                  candidates down to a single stable time step.  Users
     *                  may employ a custom functor compatible with
     *                  delta_t_reducer to add logic for monitoring or
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
    controller(
            const method_interface<element>& m,
            Reducer& reducer,
            const linear_operator<LinearA,LinearB>& L,
            const typename traits::component<element>::type chi,
            const nonlinear_operator<NonlinearB>& N,
            StateA& a,
            StateB& b,
            typename super::time_type initial_t = 0,
            typename super::time_type min_dt = 0,
            typename super::time_type max_dt = 0)
        : super(boost::bind(&controller::stepper, this, _1),
                initial_t,
                min_dt,
                max_dt),
          m(m), reducer(reducer), L(L), chi(chi), N(N), a(a), b(b) {}

private:

    const method_interface<element>& m;
    Reducer &reducer;
    const linear_operator<LinearA,LinearB>& L;
    const typename traits::component<element>::type chi;
    const nonlinear_operator<NonlinearB>& N;
    StateA& a;
    StateB& b;

    typename super::time_type stepper(typename super::time_type max_dt)
    {
        return lowstorage::step(
                m, reducer, L, chi, N, super::current_t(), a, b, max_dt);
    }

};

/**
 * A partial specialization of the controller template for the
 * case when default delta_t_reducer behavior is desired.  Empty base class
 * optimization eliminates the delta_t_reducer instance overhead.
 */
template< typename StateA,
          typename StateB,
          typename LinearA,
          typename LinearB,
          typename NonlinearB >
class controller<StateA,StateB,void,LinearA,LinearB,NonlinearB>
    : private delta_t_reducer,
      public controller<
            StateA,StateB,delta_t_reducer,LinearA,LinearB,NonlinearB
        >
{

protected:

    typedef typename controller<
            StateA,StateB,delta_t_reducer,LinearA,LinearB,NonlinearB
        >::super super;

public:

    typedef typename controller<
            StateA,StateB,delta_t_reducer,LinearA,LinearB,NonlinearB
        >::element element;

    /**
     * Construct an instance that will advance a simulation built atop the
     * given operators and storage.
     *
     * @param m         The low storage scheme to use.
     *                  For example, \ref method in conjunction with \ref smr91.
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
    controller(
            const method_interface<element>& m,
            const linear_operator<LinearA,LinearB>& L,
            const typename traits::component<element>::type chi,
            const nonlinear_operator<NonlinearB>& N,
            StateA& a,
            StateB& b,
            typename super::time_type initial_t = 0,
            typename super::time_type min_dt = 0,
            typename super::time_type max_dt = 0)
        : delta_t_reducer(),
          controller<
                StateA,StateB,delta_t_reducer,LinearA,LinearB,NonlinearB
            >(m, *reinterpret_cast<delta_t_reducer*>(this),
              L, chi, N, a, b, initial_t, min_dt, max_dt)
    {}

};

/**
 * A helper method so the compiler can deduce the appropriate template
 * types for a controller employing a custom Reducer.
 *
 * \copydoc controller
 */
template< typename StateA,
          typename StateB,
          typename Reducer,
          typename LinearA,
          typename LinearB,
          typename NonlinearB,
          typename ChiType >
controller<StateA,StateB,Reducer,LinearA,LinearB,NonlinearB>*
make_controller(
        const method_interface<typename StateA::element>& m,
        Reducer &reducer,
        const linear_operator<LinearA,LinearB>& L,
        const ChiType chi,
        const nonlinear_operator<NonlinearB>& N,
        StateA& a,
        StateB& b,
        typename controller<
                StateA,StateB,Reducer,LinearA,LinearB,NonlinearB
            >::time_type initial_t = 0,
        typename controller<
                StateA,StateB,Reducer,LinearA,LinearB,NonlinearB
            >::time_type min_dt = 0,
        typename controller<
                StateA,StateB,Reducer,LinearA,LinearB,NonlinearB
            >::time_type max_dt = 0)
{
    return new controller<
                StateA,StateB,Reducer,LinearA,LinearB,NonlinearB
        >(m, reducer, L, chi, N, a, b, initial_t, min_dt, max_dt);
}

/**
 * A helper method so the compiler can deduce the appropriate template
 * types for a controller employing delta_t_reducer.
 *
 * \copydoc controller
 */
template< typename StateA,
          typename StateB,
          typename LinearA,
          typename LinearB,
          typename NonlinearB,
          typename ChiType >
controller<StateA,StateB,void,LinearA,LinearB,NonlinearB>*
make_controller(
        const method_interface<typename StateA::element>& m,
        const linear_operator<LinearA,LinearB>& L,
        const ChiType chi,
        const nonlinear_operator<NonlinearB>& N,
        StateA& a,
        StateB& b,
        typename controller<
                StateA,StateB,void,LinearA,LinearB,NonlinearB
            >::time_type initial_t = 0,
        typename controller<
                StateA,StateB,void,LinearA,LinearB,NonlinearB
            >::time_type min_dt = 0,
        typename controller<
                StateA,StateB,void,LinearA,LinearB,NonlinearB
            >::time_type max_dt = 0)
{
    return new controller<
            StateA,StateB,void,LinearA,LinearB,NonlinearB
        >(m, L, chi, N, a, b, initial_t, min_dt, max_dt);
}

} // namespace lowstorage

} // namespace suzerain

#endif // SUZERAIN_LOWSTORAGE_HPP
