// Copyright (C) 2012, 2013 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef AR_HPP
#define AR_HPP

/**
 * @mainpage
 *
 * \ref ar implements \ref ar "modeling tools" for autoregressive processes in
 * header-only C++.
 *
 * See the current <a
 * href="https://github.com/RhysU/ar/blob/master/README.rst"> README</a> for a
 * more detailed overview and http://github.com/RhysU/ar for project
 * information.
 */

/** @file
 * Autoregressive process modeling tools in header-only C++.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * Autoregressive process modeling tools in header-only C++.
 *
 * All routines estimate and/or evaluate autoregressive models of the form
 * \f{align}{
 *     x_n + a_1 x_{n - 1} + \dots + a_p x_{n - p} &= \epsilon_n
 *     &
 *     \epsilon_n &\sim{} N\left(0, \sigma^2_\epsilon\right)
 *     \\
 *     \sigma^2_x \left(
 *        \rho_0 + a_1 \rho_{1} + \dots + a_p \rho_{p}
 *     \right) &= \sigma^2_\epsilon
 *     &
 *     \rho_0 &= 1
 *     \\
 *     \rho_k + a_1 \rho_{k-1} + \dots + a_p \rho_{k-p} &= 0
 *     &
 *     k &\geq{} p
 * \f}
 * where \f$x_i\f$ are the process values, \f$a_i\f$ are the model parameters,
 * and \f$\rho_i\f$ are the lag \f$i\f$ autocorrelations.  The white noise
 * input process \f$\epsilon_n\f$ has variance \f$\sigma^2_\epsilon\f$.  The
 * model has output variance \f$\sigma^2_x\f$ and therefore a gain equal to
 * \f$\sigma^2_x / \sigma^2_\epsilon\f$.
 */
namespace ar
{

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Preprocessor macros to simplify implementation.
 *
 * @{
 */

/** Helper for defining GCC version checks. */
#define AR_GCC_VERSION \
    (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

/** Helper macro for implementing \ref AR_STRINGIFY. */
#define AR_STRINGIFY_HELPER(x) #x

/** Expand and stringify the provided argument. */
#define AR_STRINGIFY(x) AR_STRINGIFY_HELPER(x)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown with
 * message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define AR_ENSURE_MSGEXCEPT(expr, msg, except) \
    if (!(expr)) throw except(msg)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown
 * with message \c msg.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define AR_ENSURE_MSG(expr, msg) \
    AR_ENSURE_MSGEXCEPT(expr, msg, std::logic_error)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then a <tt>std::logic_error</tt> is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define AR_ENSURE(expr) \
    AR_ENSURE_MSG(expr, AR_STRINGIFY(expr)" false")

/**
 * Ensure that the argument-related \c expr evaluates to boolean \c true at
 * runtime.  If \c expr evaluates to boolean \c false, then a
 * <tt>std::invalid_argument</tt> is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define AR_ENSURE_ARG(expr) \
    AR_ENSURE_MSGEXCEPT(expr, AR_STRINGIFY(expr)" false", std::invalid_argument)

/**
 * Ensure that \c expr evaluates to boolean \c true at runtime.  If \c expr
 * evaluates to boolean \c false, then an exception \c except is thrown.
 *
 * This macro is intended for <tt>assert</tt>-like checks which should always
 * be performed regardless of whether or not \c NDEBUG is <tt>#define</tt>d.
 */
#define AR_ENSURE_EXCEPT(expr, except) \
    AR_ENSURE_MSGEXCEPT(expr, AR_STRINGIFY(expr)" false", except)

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Stable, one-pass algorithms for computing variances and covariances.
 *
 * @{
 */

/**
 * Compute the mean and the number of samples, N, times the population variance
 * using Welford's algorithm.  The latter quantity is effectively the centered
 * sum of squares. The algorithm is found in Knuth's TAOCP volume 2 section
 * 4.2.2.A on page 232.  The implementation follows
 * http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance.
 *
 * @param[in]  first Beginning of the input data range.
 * @param[in]  last  Exclusive end of the input data range.
 * @param[out] mean  Mean of the data in <tt>[first, last)</tt>.
 * @param[out] nvar  N times the variance of the data.
 *
 * @returns the number data values processed within <tt>[first, last)</tt>.
 */
template <typename InputIterator,
          typename OutputType1,
          typename OutputType2>
std::size_t welford_nvariance(InputIterator first,
                              InputIterator last,
                              OutputType1&   mean,
                              OutputType2&   nvar)
{
    using std::iterator_traits;
    using std::size_t;
    typedef typename iterator_traits<InputIterator>::value_type value;

    size_t N  = 1;  // Running next sample number
    value  m  = 0;  // Running mean of data thus far
    value  nv = 0;  // Running variance times the number of samples

    while (first != last)
    {
        value x  = *first++;
        value d  = x - m;
        m       += d / N++;
        nv      += d*(x - m);
    }

    mean = m;
    nvar = nv;
    return N-1;
}

/**
 * Compute the mean and population variance using Welford's algorithm.
 *
 * @param[in]  first Beginning of the input data range.
 * @param[in]  last  Exclusive end of the input data range.
 * @param[out] mean  Mean of the data in <tt>[first, last)</tt>.
 * @param[out] var   The population variance of the data.
 *
 * @returns N, the number data values processed within <tt>[first, last)</tt>.
 */
template <typename InputIterator,
          typename OutputType1,
          typename OutputType2>
std::size_t welford_variance_population(InputIterator first,
                                        InputIterator last,
                                        OutputType1&  mean,
                                        OutputType2&  var)
{
    using std::size_t;
    size_t N = welford_nvariance(first, last, mean, var);
    var /= N;
    return N;
}

/**
 * Compute the mean and sample variance using Welford's algorithm.
 *
 * @param[in]  first Beginning of the input data range.
 * @param[in]  last  Exclusive end of the input data range.
 * @param[out] mean  Mean of the data in <tt>[first, last)</tt>.
 * @param[out] var   The sample variance of the data.
 *
 * @returns N, the number data values processed within <tt>[first, last)</tt>.
 */
template <typename InputIterator,
          typename OutputType1,
          typename OutputType2>
std::size_t welford_variance_sample(InputIterator first,
                                    InputIterator last,
                                    OutputType1&  mean,
                                    OutputType2&  var)
{
    using std::size_t;
    size_t N = welford_nvariance(first, last, mean, var);
    var /= (N - 1);
    return N;
}

/**
 * Compute means and the number of samples, N, times the population covariance
 * using Welford's algorithm.  The implementation follows the covariance
 * section of http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance.
 *
 * @param[in]  first1 Beginning of the first input range.
 * @param[in]  last1  Exclusive end of first input range.
 * @param[in]  first2 Beginning of the second input range.
 * @param[out] mean1  Mean of the first data set.
 * @param[out] mean2  Mean of the second data set.
 * @param[out] ncovar N times the covariance of the two sets.
 *
 * @returns the number data values processed within <tt>[first1, last1)</tt>.
 */
template <typename InputIterator1,
          typename InputIterator2,
          typename OutputType1,
          typename OutputType2,
          typename OutputType3>
std::size_t welford_ncovariance(InputIterator1 first1,
                                InputIterator1 last1,
                                InputIterator2 first2,
                                OutputType1&   mean1,
                                OutputType2&   mean2,
                                OutputType3&   ncovar)
{
    using std::iterator_traits;
    using std::size_t;
    typedef typename iterator_traits<InputIterator1>::value_type value1;
    typedef typename iterator_traits<InputIterator2>::value_type value2;

    size_t      N  = 1;  // Running next sample number
    value1      m1 = 0;  // Running mean of first data set thus far
    value2      m2 = 0;  // Running mean of second data set thus far
    OutputType3 nc = 0;  // Running covariance times the number of samples

    while (first1 != last1)
    {
        value1 x1  = *first1++;
        value1 d1  = x1 - m1;
        m1        += d1 / N;

        value2 x2  = *first2++;
        value2 d2  = x2 - m2;
        m2        += d2 / N;

        nc += d1*(x2 - m2);

        ++N;
    }

    mean1  = m1;
    mean2  = m2;
    ncovar = nc;
    return N-1;
}

/**
 * Compute means and the population covariance using Welford's algorithm.
 *
 * @param[in]  first1 Beginning of the first input range.
 * @param[in]  last1  Exclusive end of first input range.
 * @param[in]  first2 Beginning of the second input range.
 * @param[out] mean1  Mean of the first data set.
 * @param[out] mean2  Mean of the second data set.
 * @param[out] covar  The covariance of the two sets.
 *
 * @returns the number data values processed within <tt>[first1, last1)</tt>.
 */
template <typename InputIterator1,
          typename InputIterator2,
          typename OutputType1,
          typename OutputType2,
          typename OutputType3>
std::size_t welford_covariance_population(InputIterator1 first1,
                                          InputIterator1 last1,
                                          InputIterator2 first2,
                                          OutputType1&   mean1,
                                          OutputType2&   mean2,
                                          OutputType3&   covar)
{
    using std::size_t;
    size_t N = welford_ncovariance(first1, last1, first2, mean1, mean2, covar);
    covar /= N;
    return N;
}

/**
 * Compute means and the sample covariance using Welford's algorithm.
 *
 * @param[in]  first1 Beginning of the first input range.
 * @param[in]  last1  Exclusive end of first input range.
 * @param[in]  first2 Beginning of the second input range.
 * @param[out] mean1  Mean of the first data set.
 * @param[out] mean2  Mean of the second data set.
 * @param[out] covar  The covariance of the two sets.
 *
 * @returns the number data values processed within <tt>[first1, last1)</tt>.
 */
template <typename InputIterator1,
          typename InputIterator2,
          typename OutputType1,
          typename OutputType2,
          typename OutputType3>
std::size_t welford_covariance_sample(InputIterator1 first1,
                                      InputIterator1 last1,
                                      InputIterator2 first2,
                                      OutputType1&   mean1,
                                      OutputType2&   mean2,
                                      OutputType3&   covar)
{
    using std::size_t;
    size_t N = welford_ncovariance(first1, last1, first2, mean1, mean2, covar);
    covar /= (N - 1);
    return N;
}

/**
 * Compute the inner product of <tt>[first, last)</tt> with itself using \ref
 * welford_nvariance.  Welford's algorithm is combined with the linearity of
 * the expectation operator to compute a more expensive but also more
 * numerically stable result than can be had using <tt>std::inner_product</tt>.
 *
 * @param[in] first Beginning of the input data range.
 * @param[in] last  Exclusive end of the input data range.
 * @param[in] init  Initial value, often zero, establishing the result type.
 *
 * @returns The inner product of <tt>[first, last)</tt> with itself.
 */
template <typename InputIterator,
          typename ValueType>
ValueType welford_inner_product(InputIterator first,
                                InputIterator last,
                                ValueType     init)
{
    typename std::iterator_traits<InputIterator>::value_type mean;
    ValueType nvar;
    const std::size_t N = welford_nvariance(first, last, mean, nvar);
    return init + (nvar + N*(mean*mean));
}

/**
 * Compute the inner product of <tt>[first1, last1)</tt> against <tt>[first2,
 * ...)</tt> using \ref welford_ncovariance.  Welford's algorithm is combined
 * with the linearity of the expectation operator to compute a more expensive
 * but also numerically stable result than can be had using
 * <tt>std::inner_product</tt>.
 *
 * @param[in] first1 Beginning of the first input range.
 * @param[in] last1  Exclusive end of first input range.
 * @param[in] first2 Beginning of the second input range.
 * @param[in] init   Initial value, often zero, establishing the result type.
 *
 * @returns The inner product of <tt>[first1, last1)</tt> against
 *          <tt>[first2, ...)</tt>.
 */
template <typename InputIterator1,
          typename InputIterator2,
          typename ValueType>
ValueType welford_inner_product(InputIterator1 first1,
                                InputIterator1 last1,
                                InputIterator2 first2,
                                ValueType      init)
{
    typename std::iterator_traits<InputIterator1>::value_type mean1;
    typename std::iterator_traits<InputIterator2>::value_type mean2;
    ValueType ncovar;
    const std::size_t N = welford_ncovariance(
            first1, last1, first2, mean1, mean2, ncovar);
    return init + (ncovar + N*(mean1*mean2));
}

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Algorithms for autoregressive parameter estimation and manipulation.
 *
 * @{
 */

#if _MSC_VER > 1400
# pragma float_control(push)
# pragma float_control(precise, on)
#endif

/**
 * Robustly compute negative one half the reflection coefficient assuming
 * \f$\vec{a}\f$ and \f$\vec{b}\f$ contain real-valued backward and forward
 * prediction error sequences, respectively.
 *
 * @param[in] a_first Beginning of the first input range \f$\vec{a}\f$.
 * @param[in] a_last  Exclusive end of first input range \f$\vec{a}\f$.
 * @param[in] b_first Beginning of the second input range \f$\vec{b}\f$.
 *
 * @return \f$\frac{\vec{a}\cdot\vec{b}}
 *                 {\vec{a}\cdot\vec{a} + \vec{b}\cdot\vec{b}}\f$.
 *
 * @see Wikipedia's article on <a href="">Kahan summation</a> for
 *      background on how the accumulation error is reduced in the result.
 */
template <typename ValueType,
          typename InputIterator1,
          typename InputIterator2>
ValueType
#if (AR_GCC_VERSION > 40305)
    __attribute__((__optimize__("no-associative-math")))
#endif
negative_half_reflection_coefficient(InputIterator1 a_first,
                                     InputIterator1 a_last,
                                     InputIterator2 b_first)
#if (AR_GCC_VERSION > 40305) ||  (_MSC_VER > 1400)
{
    ValueType ns = 0, nt, nc = 0, ny;  // Kahan numerator accumulation
    ValueType ds = 0, dt, dc = 0, dy;  // Kahan denominator accumulation

    while (a_first != a_last)
    {
        ValueType xa = *a_first++;     // Denominator: \vec{a}\cdot\vec{a}
        dy = (xa * xa) - dc;
        dt = ds + dy;
        dc = (dt - ds) - dy;
        ds = dt;

        ValueType xb = *b_first++;     // Denominator: \vec{b}\cdot\vec{b}
        dy = (xb * xb) - dc;
        dt = ds + dy;
        dc = (dt - ds) - dy;
        ds = dt;

        ny = (xa * xb) - nc;           // Numerator:   \vec{a}\cdot\vec{b}
        nt = ns + ny;
        nc = (nt - ns) - ny;
        ns = nt;
    }

    return (ns + nc) / (ds + dc);      // Correct final sums and form ratio
}
#else
#warning Using Non-Kahan version of ar::negative_half_reflection_coefficient.
{
    ValueType ns = 0;
    ValueType ds = 0;

    while (a_first != a_last)
    {
        ValueType xa  = *a_first++;
        ValueType xb  = *b_first++;
        ns           += xa * xb;
        ds           += xa * xa + xb * xb;
    }

    return ns / ds;
}
#endif

#if _MSC_VER > 1400
# pragma float_control(pop)
#endif

/**
 * Fit an autoregressive model to stationary time series data using %Burg's
 * method.  That is, find coefficients \f$a_i\f$ such that the sum of the
 * squared errors in the forward predictions \f$x_n = -a_1 x_{n-1} - \dots -
 * a_p x_{n-p}\f$ and backward predictions \f$x_n = -a_1 x_{n+1} - \dots - a_p
 * x_{n+p}\f$ are both minimized.  Either a single model of given order or a
 * hierarchy of models up to and including a maximum order may fit.
 *
 * The input data \f$\vec{x}\f$ are read from <tt>[data_first, data_last)</tt>
 * in a single pass.  The mean is computed, returned in \c mean, and \e
 * removed from further consideration whenever \c subtract_mean is true.
 * The estimated model parameters \f$a_i\f$ are output using \c params_first
 * with the behavior determined by the amount of data read, <tt>maxorder</tt>,
 * and the \c hierarchy flag:
 * <ul>
 *     <li>If \c hierarchy is \c false, only the \f$a_1, \dots,
 *         a_\text{maxorder}\f$ parameters for an AR(<tt>maxorder</tt>) process
 *         are output.</li>
 *     <li>If \c hierarchy is \c true, the <tt>maxorder*(maxorder+1)/2</tt>
 *         parameters \f$a_1, \dots, a_m\f$ for models AR(0), AR(1), AR(2),
 *         ..., AR(maxorder) are output.  Notice AR(0) has no parameters.
 *         </li>
 * </ul>
 * Note that the latter case is \e always computed; the \c hierarchy flag
 * merely controls what is output.  In both cases, the maximum order is limited
 * by the number of data samples provided and is output to \c maxorder.
 *
 * One mean squared discrepancy \f$\sigma^2_\epsilon\f$, also called the
 * innovation variance, and gain, defined as \f$\sigma^2_x /
 * \sigma^2_\epsilon\f$, are output for each model, including the trivial
 * zeroth order model when \c maxorder is zero or \c hierarchy is \c true,
 * using \c sigma2e_first and \c gain_first.  The autocorrelations for lags
 * <tt>[0,k]</tt> are output using \c autocor_first.  When \c hierarchy is \c
 * true, only lags <tt>[0,m]</tt> should be applied for some AR(<tt>m</tt>)
 * model.  Outputting the lag \c k autocorrelation is technically redundant as
 * it may be computed from \f$a_i\f$ and lags <tt>0, ..., k-1</tt>.
 * Autocovariances may be computed by multiplying the autocorrelations by the
 * gain times \f$\sigma^2_\epsilon\f$.
 *
 * The software aspects of the implementation differs from many other sources.
 * In particular,
 * <ul>
 *     <li>iterators are employed,</li>
 *     <li>the working precision is selectable using \c mean,</li>
 *     <li>the mean squared discrepancy calculation has been added,</li>
 *     <li>some loop index transformations have been performed,</li>
 *     <li>working storage may be passed into the method to reduce allocations
 *     across many invocations, and</li>
 *     <li>and all lower order models may be output during the recursion using
 *     \c hierarchy.</li>
 * </ul>
 * Gain and autocorrelation calculations have been added based on sections 5.2
 * and 5.3 of Broersen, P.  M.  T. Automatic autocorrelation and spectral
 * analysis.  Springer, 2006.  http://dx.doi.org/10.1007/1-84628-329-9.  The
 * classical algorithm, rather than the variant using denominator recursion due
 * to Andersen (http://dx.doi.org/10.1109/PROC.1978.11160), has been chosen as
 * the latter can be numerically unstable.
 *
 * @param[in]     data_first    Beginning of the input data range.
 * @param[in]     data_last     Exclusive end of the input data range.
 * @param[out]    mean          Mean of data.
 * @param[in,out] maxorder      On input, the maximum model order desired.
 *                              On output, the maximum model order computed.
 * @param[out]    params_first  Model parameters for a single model or
 *                              for an entire hierarchy of models.  At most
 *                              <tt>!hierarchy ? maxorder :
 *                              maxorder*(maxorder+1)/2</tt> values will be
 *                              output.
 * @param[out]    sigma2e_first The mean squared discrepancy for only
 *                              AR(<tt>maxorder</tt>) or for an entire
 *                              hierarchy.  Either one or at most
 *                              <tt>maxorder + 1</tt> values will be output.
 * @param[out]    gain_first    The model gain for only AR(<tt>maxorder</tt>)
 *                              or an entire hierarchy.  Either one or at most
 *                              <tt>maxorder + 1</tt> values will be output.
 * @param[out]    autocor_first Lag one through lag maxorder autocorrelations.
 *                              At most <tt>maxorder + 1</tt> values will be
 *                              output.
 * @param[in]     subtract_mean Should \c mean be subtracted from the data?
 * @param[in]     hierarchy     Should the entire hierarchy of estimated
 *                              models be output?
 * @param[in]     f             Working storage.  Reuse across invocations
 *                              may speed execution by avoiding allocations.
 * @param[in]     b             Working storage similar to \c f.
 * @param[in]     Ak            Working storage similar to \c f.
 * @param[in]     ac            Working storage similar to \c f.
 *
 * @returns the number data values processed within
 *          <tt>[data_first, data_last)</tt>.
 */
template <class InputIterator,
          class Value,
          class OutputIterator1,
          class OutputIterator2,
          class OutputIterator3,
          class OutputIterator4,
          class Vector>
std::size_t burg_method(InputIterator   data_first,
                        InputIterator   data_last,
                        Value&          mean,
                        std::size_t&    maxorder,
                        OutputIterator1 params_first,
                        OutputIterator2 sigma2e_first,
                        OutputIterator3 gain_first,
                        OutputIterator4 autocor_first,
                        const bool      subtract_mean,
                        const bool      hierarchy,
                        Vector&         f,
                        Vector&         b,
                        Vector&         Ak,
                        Vector&         ac)
{
    using std::bind2nd;
    using std::copy;
    using std::distance;
    using std::fill;
    using std::inner_product;
    using std::min;
    using std::minus;
    using std::size_t;

    // Initialize f from [data_first, data_last) and fix number of samples
    f.assign(data_first, data_last);
    const size_t N = f.size();

    // Stably compute the incoming data's mean and population variance
    mean = 0;
    Value sigma2e = 0;
    welford_variance_population(f.begin(), f.end(), mean, sigma2e);

    // When requested, subtract the just-computed mean from the data.
    // Adjust, if necessary, to make sigma2e the second moment.
    if (subtract_mean)
    {
        transform(f.begin(), f.end(), f.begin(), bind2nd(minus<Value>(), mean));
    }
    else
    {
        sigma2e += mean*mean;
    }

    // At most maxorder N-1 can be fit from N samples.  Beware N is unsigned.
    maxorder = (N == 0) ? 0 : min(static_cast<size_t>(maxorder), N-1);

    // Output sigma2e and gain for a zeroth order model, if requested.
    Value gain = 1;
    if (hierarchy || maxorder == 0)
    {
        *sigma2e_first++ = sigma2e;
        *gain_first++    = gain;
    }

    // Initialize and perform Burg recursion
    if (maxorder) b = f;  // Copy iff non-trivial work required
    Ak.assign(maxorder + 1, Value(0));
    Ak[0] = 1;
    ac.clear();
    ac.reserve(maxorder);
    for (size_t kp1 = 1; kp1 <= maxorder; ++kp1)
    {
        // Compute mu from f, b, and Dk and then update sigma2e and Ak using mu
        // Afterwards, Ak[1:kp1] contains AR(k) coefficients by the recurrence
        // By the recurrence, Ak[kp1] will also be the reflection coefficient
        Value mu = -2 * negative_half_reflection_coefficient<Value>(
                f.begin() + kp1, f.end(), b.begin());

        sigma2e *= (1 - mu*mu);
        for (size_t n = 0; n <= kp1/2; ++n)
        {
            Value t1 = Ak[n] + mu*Ak[kp1 - n];
            Value t2 = Ak[kp1 - n] + mu*Ak[n];
            Ak[n] = t1;
            Ak[kp1 - n] = t2;
        }

        // Update the gain per Broersen 2006 equation (5.25)
        gain *= 1 / (1 - Ak[kp1]*Ak[kp1]);

        // Compute and output the next autocorrelation coefficient
        // See Broersen 2006 equations (5.28) and (5.31) for details
        ac.push_back(-inner_product(ac.rbegin(), ac.rend(),
                                    Ak.begin() + 1, Ak[kp1]));

        // Output parameters and the input and output variances when requested
        if (hierarchy || kp1 == maxorder)
        {
            params_first = copy(Ak.begin() + 1, Ak.begin() + kp1 + 1,
                                params_first);
            *sigma2e_first++ = sigma2e;
            *gain_first++    = gain;
        }

        // Update f and b for the next iteration if another remains
        if (kp1 < maxorder)
        {
            for (size_t n = 0; n < N - kp1; ++n)
            {
                Value t1 = f[n + kp1] + mu*b[n];
                Value t2 = b[n] + mu*f[n + kp1];
                f[n + kp1] = t1;
                b[n] = t2;
            }
        }
    }

    // Output the lag [0,maxorder] autocorrelation coefficients in single pass
    *autocor_first++ = 1;
    copy(ac.begin(), ac.end(), autocor_first);

    // Return the number of values processed in [data_first, data_last)
    return N;
}

/** \copydoc burg_method(InputIterator,InputIterator,Value&,std::size_t&,OutputIterator1,OutputIterator2,OutputIterator3,OutputIterator4,const bool,const bool,Vector&,Vector&,Vector&,Vector&) */
template <class InputIterator,
          class Value,
          class OutputIterator1,
          class OutputIterator2,
          class OutputIterator3,
          class OutputIterator4>
std::size_t burg_method(InputIterator     data_first,
                        InputIterator     data_last,
                        Value&            mean,
                        std::size_t&      maxorder,
                        OutputIterator1   params_first,
                        OutputIterator2   sigma2e_first,
                        OutputIterator3   gain_first,
                        OutputIterator4   autocor_first,
                        const bool        subtract_mean = false,
                        const bool        hierarchy     = false)
{
    using std::vector;
    vector<Value> f, b, Ak, ac; // Working storage

    return burg_method(data_first, data_last, mean, maxorder,
                       params_first, sigma2e_first, gain_first,
                       autocor_first, subtract_mean, hierarchy,
                       f, b, Ak, ac);
}

// Type erasure for NoiseGenerator parameters within predictor.
// Either std::tr1::function or boost::function would better provide the
// desired capability but both add additional, undesired dependencies.
namespace
{

/** Abstract base class for NoiseGenerator-related type erasure. */
template <typename Value>
struct nullary
{
    virtual ~nullary() {}
    virtual Value operator()() = 0;
    virtual nullary* clone()   = 0;
};

/** A nullary function always returning zero. */
template<typename Value>
struct nullary_impl0 : public nullary<Value>
{
    Value operator()()
    {
        return 0;
    }
    nullary_impl0* clone()
    {
        return new nullary_impl0();
    }
};

/** A nullary function always invoking t(). */
template<typename Value, class T>
struct nullary_impl1 : public nullary<Value>
{
    nullary_impl1(T t) : t(t) {}
    Value operator()()
    {
        return t();
    }
    nullary_impl1* clone()
    {
        return new nullary_impl1(t);
    }
    T t;
};

}

/**
 * Simulate an autoregressive model process with an InputIterator interface.
 */
template <typename Value, typename Index = std::size_t>
class predictor
    : public std::iterator<std::input_iterator_tag, Value,
      std::ptrdiff_t, const Value*, const Value&>
{
private:
    typedef std::iterator<std::input_iterator_tag, Value,
            std::ptrdiff_t, const Value*, const Value&> base;

public:
    typedef typename base::difference_type   difference_type;
    typedef typename base::iterator_category iterator_category;
    typedef typename base::pointer           pointer;
    typedef typename base::reference         reference;
    typedef typename base::value_type        value_type;

    /** Singular instance marking prediction index \c n. */
    explicit predictor(Index n = 0) : n(n), d(), g(0), xn()
    {
#ifndef NDEBUG
        using std::numeric_limits;
        if (numeric_limits<Value>::has_quiet_NaN)
            xn = numeric_limits<Value>::quiet_NaN();
#endif
    }

    /**
     * Iterate on the process \f$x_n + a_1 x_{n - 1} + \dots + a_p x_{n - p} =
     * 0\f$.  Presumably \ref initial_conditions will be used to specify some
     * initial state as otherwise the process is identically zero.  The process
     * order \f$p\f$ is set by <tt>std::distance(params_first,
     * params_last)</tt>.
     *
     * @param params_first  Beginning of the process parameter range
     *                      starting with \f$a_1\f$.
     * @param params_last   End of the process parameter range.
     */
    template <class RandomAccessIterator>
    predictor(RandomAccessIterator params_first,
              RandomAccessIterator params_last)
        : n(0),
          d(2*std::distance(params_first, params_last), 0),
          g(new nullary_impl0<Value>()),
          xn((*g)())
    {
        // Finish preparing d = [ a_p, ..., a_1, 0, ..., 0 ]
        using std::vector;
        typename vector<Value>::size_type i = d.size() / 2;
        while (i --> 0) d[i] = *params_first++;

        // Now x_n = 0 because x_{n-p} = ... = x_{n-1} = 0 by construction.
    }

    /**
     * Iterate on the process \f$x_n + a_1 x_{n - 1} + \dots + a_p x_{n - p} =
     * \epsilon_n\f$ given zero initial conditions.  The process order \f$p\f$
     * is set by <tt>std::distance(params_first,params_last)</tt>.  The class
     * <tt>std::tr1::variate_generator</tt> may be helpful in constructing
     * normally distributed input.
     *
     * @param params_first  Beginning of the process parameter range
     *                      starting with \f$a_1\f$.
     * @param params_last   End of the process parameter range.
     * @param generator     A nullary callback for generating \f$\epsilon_n\f$.
     *                      For example, a random number generator distributed
     *                      like \f$N\left(0, \sigma^2_\epsilon\right)\f$.
     */
    template <class RandomAccessIterator,
              class NoiseGenerator>
    predictor(RandomAccessIterator params_first,
              RandomAccessIterator params_last,
              NoiseGenerator generator)
        : n(0),
          d(2*std::distance(params_first, params_last), 0),
          g(new nullary_impl1<Value,NoiseGenerator>(generator)),
          xn((*g)())
    {
        // Finish preparing d = [ a_p, ..., a_1, 0, ..., 0 ]
        using std::vector;
        typename vector<Value>::size_type i = d.size() / 2;
        while (i --> 0) d[i] = *params_first++;

        // Here x_0 = \epsilon_0 because x_{0-p} = ... = x_{0-1} = 0.
    }

    /** Copy constructor */
    predictor(const predictor& other)
        : n(other.n),
          d(other.d),
          g(other.g ? other.g->clone() : 0),
          xn(other.xn)
    {}

    /** Assignment operator */
    predictor& operator=(const predictor& other)
    {
        if (this != &other)
        {
            nullary<Value> *tmp = 0;
            try
            {
                tmp = other.g ? other.g->clone() : 0;
            }
            catch (...)
            {
                delete tmp;
                throw;
            }
            base::operator=(other);
            n = other.n;
            d = other.d;
            delete g;
            g = tmp;
            xn = other.xn;
        }
        return *this;
    }

    /** Destructor */
    ~predictor()
    {
        delete g;
    }

    /**
     * Specify process initial conditions \f$x_{n-1}, \dots, x_{n-p}\f$ where
     * \f$p\f$ is the process order fixed by the constructor.  The simulation
     * index \f$n\f$ is reset to zero and, optionally, \f$x_0\f$ is additively
     * adjusted by \c x0adjust.
     *
     * @param initial_first Beginning of the initial condition range
     *                      \f$x_{n-1}, \dots, x_{n-p}\f$
     *                      which must contain \f$p\f$ values.
     * @param x0adjust      An additive adjustment made to \f$\epsilon_0\f$.
     */
    template <class InputIterator>
    predictor& initial_conditions(InputIterator initial_first,
                                  const Value x0adjust = 0)
    {
        // Zero the simulation time.
        n = 0;

        // Set d = [ a_p, ..., a_1, x_{n-p}, ..., x_{n-1} ]
        using std::vector;
        typename vector<Value>::size_type i = d.size();
        typename vector<Value>::size_type p = i / 2;
        while (i --> p) d[i] = *initial_first++;

        // Make x_n := - a_p*x_{n-p} - ... - a_1*x_{n-1} + x_n + x0adjust.
        // By design, x_n was whatever it happened to be.
        using std::inner_product;
        xn += x0adjust;
        xn  = -inner_product(d.begin(), d.begin() + p, d.begin() + p, -xn);

        return *this;
    }

    // Concept: InputIterator

    /** Prefix increment. */
    predictor& operator++()
    {
        using std::distance;
        using std::inner_product;
        using std::vector;

        if (g)
        {
            typename vector<Value>::size_type p = d.size() / 2;
            if (p)
            {
                // Make x_n = - a_p*x_{n-p} - ... - a_1*x_{n-1} + \epsilon_n
                // by (conceptually) storing previously computed x_n into
                // circular buffer, updating ++n, and computing x_{n+1}.
                typename vector<Value>::iterator ab = d.begin();
                typename vector<Value>::iterator xb = ab + p;
                typename vector<Value>::iterator c  = xb + n % p;
                typename vector<Value>::iterator xe = d.end();
                *c++ =  xn;
                xn   =  inner_product(c,  xe, ab,                   -(*g)());
                xn   = -inner_product(xb,  c, ab + distance(c, xe),   xn   );
            }
            else
            {
                xn = (*g)();
            }
        }
        else
        {
#ifndef NDEBUG
            using std::numeric_limits;
            if (numeric_limits<Value>::has_quiet_NaN)
                xn = numeric_limits<Value>::quiet_NaN();
#endif
        }

        ++n;

        return *this;
    }

    /** Postfix increment. */
    predictor operator++(int) // Postfix increment
    {
        predictor t(*this);
        ++*this;
        return t;
    }

    /** Obtain the process prediction \f$x_n\f$. */
    reference operator* () const
    {
        return xn;
    }

    // Concept: EqualityComparable

    /** Check if two iterators represent the same simulation time. */
    bool operator== (const predictor& other) const
    {
        return n == other.n;
    }

    /** Check if two iterators represent different simulation times. */
    bool operator!= (const predictor& other) const
    {
        return !(*this == other);
    }

private:

    /** Running prediction index. */
    Index n;

    /**
     * State vector keeping \f$\tilde{a}\f$ in <tt>[d, d+p)</tt> and
     * \f$x_{n-p},\dots,x_{n-2},x_{n-1}\f$ in <tt>[d+p,d+p+p)/<tt> in circular
     * buffer fashion where <tt>p = d.size()/2</tt> is the autoregressive
     * process order.  The current location in circular buffer is maintained
     * using <tt>d + p + n % p</tt>.
     */
    typename std::vector<Value> d;

    /**
     * Noise generator used at every step.
     * An instance is singular whenever <tt>!g</tt>.
     */
    nullary<Value> *g;

    /**
     * Prediction at current index \c n.  Computed on advance to permit
     * repeated inexpensive dereferencing and <tt>*i++</tt> usage.
     */
    Value xn;
};

/**
 * Construct an iterator over the autocorrelation function \f$\rho_k\f$ given
 * process parameters and initial conditions.
 *
 * @param params_first  Beginning of range containing \f$a_1,\dots,a_p\f$.
 * @param params_last   Exclusive ending of the parameter range.
 * @param gain          The model gain \f$\sigma^2_x / \sigma^2_\epsilon\f$.
 * @param autocor_first Beginning of range containing \f$\rho_0,...\rho_p\f$.
 *
 * @return An InputIterator across the autocorrelation function starting with
 *         \f$\rho_0\f$.
 */
template <class RandomAccessIterator,
          class InputIterator,
          class Value>
predictor<typename std::iterator_traits<RandomAccessIterator>::value_type>
autocorrelation(RandomAccessIterator params_first,
                RandomAccessIterator params_last,
                Value                gain,
                InputIterator        autocor_first)
{
    predictor<typename std::iterator_traits<
            RandomAccessIterator
        >::value_type> p(params_first, params_last);
    p.initial_conditions(++autocor_first, 1 / gain);
    return p;
}

/**
 * Compute the decorrelation time for variance of the mean given
 * autocorrelation details.  That is, compute
 * \f{align}{
 *     T_0 &= 1 + 2 \sum_{i=1}^{N} \left(1 - \frac{i}{N}\right) \rho_i
 * \f}
 * following Trenberth, K. E. "Some effects of finite sample size and
 * persistence on meteorological statistics. Part I: Autocorrelations." Monthly
 * Weather Review 112 (1984).
 * http://dx.doi.org/10.1175/1520-0493(1984)112%3C2359:SEOFSS%3E2.0.CO;2
 *
 * Rather than \f$\rho\f$, \f$\left|\rho\right|\f$ may be used in the
 * definition of \f$T_0\f$ to better approximate the "decay of the correlation
 * envelope" according to section 17.1.5 of Hans von Storch and Francis W.
 * Zwiers.  Statistical analysis in climate research. Cambridge University
 * Press, March 2001. ISBN 978-0521012300.  Doing so is more robust for
 * oscillatory processes and always provides a larger, more conservative
 * estimate of \f$T_0\f$.
 *
 * @param N        Maximum lag used to compute the autocorrelation.
 * @param rho      A \ref predictor iterating over the \ref autocorrelation.
 * @param abs_rho  Use \f$\left|\rho\right|\f$ when calculating \f$T_0\f$?
 *
 * @return The decorrelation time \f$T_0\f$ assuming \f$\Delta{}t=1\f$.
 */
template <class Value>
Value decorrelation_time(const std::size_t N,
                         predictor<Value> rho,
                         const bool abs_rho = false)
{
    using std::abs;
    using std::size_t;

    Value T0 = *rho++;

    const Value twoinvN = Value(2) / N;
    if (abs_rho) {
        for (size_t i = 1; i <= N; ++i, ++rho)
            T0 += (2 - i*twoinvN) * abs(*rho);
    } else {
        for (size_t i = 1; i <= N; ++i, ++rho)
            T0 += (2 - i*twoinvN) * (*rho);
    }

    return T0;
}

/**
 * Compute the decorrelation time for a covariance given autocorrelation
 * details from two processes.  That is, compute
 * \f{align}{
 *     T_0 &= 1 + 2 \sum_{i=1}^{N} \left(1 - \frac{i}{N}\right)
 *                                 \rho_{1,i} \rho_{2,i}
 * \f}
 * following Trenberth, K. E. "Some effects of finite sample size and
 * persistence on meteorological statistics. Part I: Autocorrelations." Monthly
 * Weather Review 112 (1984).
 * http://dx.doi.org/10.1175/1520-0493(1984)112%3C2359:SEOFSS%3E2.0.CO;2
 *
 * Rather than \f$\rho\f$, \f$\left|\rho\right|\f$ may be used in the
 * definition of \f$T_0\f$ to better approximate the "decay of the correlation
 * envelope" according to section 17.1.5 of Hans von Storch and Francis W.
 * Zwiers.  Statistical analysis in climate research. Cambridge University
 * Press, March 2001. ISBN 978-0521012300.  Doing so is more robust for
 * oscillatory processes and always provides a larger, more conservative
 * estimate of \f$T_0\f$.
 *
 * @param N        Maximum lag used to compute the autocorrelation.
 * @param rho1     A \ref predictor iterating over the \ref autocorrelation
 *                 for the first process.
 * @param rho2     A \ref predictor iterating over the \ref autocorrelation
 *                 for the first process.
 * @param abs_rho  Use \f$\left|\rho\right|\f$ when calculating \f$T_0\f$?
 *
 * @return The decorrelation time \f$T_0\f$ assuming \f$\Delta{}t=1\f$.
 */
template <class Value>
Value decorrelation_time(const std::size_t N,
                         predictor<Value> rho1,
                         predictor<Value> rho2,
                         const bool abs_rho = false)
{
    using std::abs;
    using std::size_t;

    Value T0 = 1;
    ++rho1;
    ++rho2;

    const Value twoinvN = Value(2) / N;
    if (abs_rho) {
        for (size_t i = 1; i <= N; ++i, ++rho1, ++rho2)
            T0 += (2 - i*twoinvN) * abs((*rho1) * (*rho2));
    } else {
        for (size_t i = 1; i <= N; ++i, ++rho1, ++rho2)
            T0 += (2 - i*twoinvN) * ((*rho1) * (*rho2));
    }

    return T0;
}

/**
 * Solve a Toeplitz set of linear equations.  That is, find \f$s_{n+1}\f$
 * satisfying
 * \f[
 *      L_{n+1} s_{n+1} = d_{n+1}
 *      \mbox{ where }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of the
 * problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  A symmetric
 * Toeplitz solve can be performed by having \f$\vec{a}\f$ and \f$\vec{r}\f$
 * iterate over the same data.  The Hermitian case requires two buffers with
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.  The working precision
 * is fixed by the \c value_type of \c d_first.
 *
 * The algorithm is from Zohar, Shalhav. "The Solution of a Toeplitz Set of
 * Linear Equations." J. ACM 21 (April 1974): 272-276.
 * http://dx.doi.org/10.1145/321812.321822.  It has complexity like
 * <tt>O(2*(n+1)^2)</tt>.  Zohar improved upon earlier work from Page 1504 from
 * Trench, William F. "Weighting Coefficients for the Prediction of Stationary
 * Time Series from the Finite Past." SIAM Journal on Applied Mathematics 15
 * (November 1967): 1502-1510.  http://www.jstor.org/stable/2099503.See
 * Bunch, James R. "Stability of Methods for Solving Toeplitz Systems of
 * Equations." SIAM Journal on Scientific and Statistical Computing 6 (1985):
 * 349-364. http://dx.doi.org/10.1137/0906025 for a discussion of the
 * algorithm's stability characteristics.
 *
 * @param[in]  a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]  a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in]  r_first Beginning of the range containing \f$\vec{r}\f$.
 * @param[in]  d_first Beginning of the range containing \f$\vec{d}\f$
 *                     which should have <tt>n+1</tt> entries available.
 * @param[out] s_first Beginning of the output range to which
 *                     <tt>n+1</tt> entries will be written.
 */
template<class RandomAccessIterator,
         class InputIterator,
         class OutputIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        RandomAccessIterator r_first,
                        InputIterator        d_first,
                        OutputIterator       s_first)
{
    using std::copy;
    using std::distance;
    using std::inner_product;
    using std::invalid_argument;
    using std::iterator_traits;
    using std::reverse_iterator;
    using std::vector;

    // Tildes indicate transposes while hats indicate reversed vectors.

    // InputIterator::value_type determines the working precision
    typedef typename iterator_traits<InputIterator>::value_type value_type;
    typedef vector<value_type> vector_type;
    typedef typename vector_type::size_type size_type;

    // Determine problem size using [a_first,a_last) and ensure nontrivial
    typename iterator_traits<RandomAccessIterator>::difference_type dist
            = distance(a_first, a_last);
    if (dist < 1) throw invalid_argument("distance(a_first, a_last) < 1");
    const size_type n = static_cast<size_type>(dist);

    // Allocate working storage and set initial values for recursion:
    vector_type s;    s   .reserve(n+1); s   .push_back( *d_first);
    vector_type ehat; ehat.reserve(n  ); ehat.push_back(-a_first[0]);
    vector_type g;    g   .reserve(n  ); g   .push_back(-r_first[0]);
    value_type lambda  = 1 - a_first[0]*r_first[0];

    // Though recursive updates to s and g can be done in-place, updates to
    // ehat seemingly require one additional vector for storage:
    //
    // "This sequence of computations is economical of storage.  It is only
    // necessary to retain quantities computed at level m - 1 until the
    // computations at level m are complete." [Trench1967, page 1504]
    vector_type next_ehat; next_ehat.reserve(n);

    // Recursion for i = {1, 2, ..., n - 1}:
    for (size_type i = 1; i < n; ++i)
    {

        reverse_iterator<RandomAccessIterator> rhat_first(r_first + i);

        // \theta_i =  \delta_{i+1}  - \tilde{s}_i \hat{r}_i
        const value_type neg_theta = inner_product(
                s.begin(), s.end(), rhat_first, value_type(-(*++d_first)));

        // \eta_i   = -\rho_{-(i+1)} - \tilde{a}_i \hat{e}_i
        const value_type neg_eta   = inner_product(
                ehat.begin(), ehat.end(), a_first, value_type(a_first[i]));

        // \gamma_i = -\rho_{i+1}    - \tilde{g}_i \hat{r}_i
        const value_type neg_gamma = inner_product(
                g.begin(), g.end(), rhat_first, value_type(r_first[i]));

        /*
         * s_{i+1} = \bigl(\begin{smallmatrix}
         *              s_i + (\theta_i/\lambda_i) \hat{e}_i \\
         *              \theta_i/\lambda_i
         *          \end{smallmatrix}\bigr)
         *
         * \hat{e}_{i+1} = \bigl(\begin{smallmatrix}
         *                     \eta_i/\lambda_i \\
         *                     \hat{e}_i + (\eta_i/\lambda_i) g_i
         *                 \end{smallmatrix}\bigr)
         *
         * g_{i+1} = \bigl(\begin{smallmatrix}
         *               g_i + (\gamma_i/\lambda_i) \hat{e}_i \\
         *               \gamma_i/\lambda_i
         *           \end{smallmatrix}\bigr)
         */
        const value_type theta_by_lambda = -neg_theta/lambda;
        const value_type   eta_by_lambda = -neg_eta  /lambda;
        const value_type gamma_by_lambda = -neg_gamma/lambda;
        next_ehat.clear();
        next_ehat.push_back(eta_by_lambda);
        for (size_type j = 0; j < i; ++j)
        {
            s[j] += theta_by_lambda*ehat[j];
            next_ehat.push_back(ehat[j] + eta_by_lambda*g[j]);
            g[j] += gamma_by_lambda*ehat[j];
        }
        s.push_back(theta_by_lambda);
        g.push_back(gamma_by_lambda);
        ehat.swap(next_ehat);

        // \lambda_{i+1} = \lambda_i - \eta_i \gamma_i / \lambda_i
        lambda -= neg_eta*neg_gamma/lambda;
    }

    // Recursion for i = n differs slightly per Zohar's "Last Computed Values"
    // Computing g_n above was unnecessary but the incremental expense is small
    {
        reverse_iterator<RandomAccessIterator> rhat_first(r_first + n);

        // \theta_n =  \delta_{n+1}  - \tilde{s}_n \hat{r}_n
        const value_type neg_theta = inner_product(
                s.begin(), s.end(), rhat_first, value_type(-(*++d_first)));

        /*
         * s_{n+1} = \bigl(\begin{smallmatrix}
         *              s_n + (\theta_n/\lambda_n) \hat{e}_n \\
         *              \theta_n/\lambda_n
         *          \end{smallmatrix}\bigr)
         */
        const value_type theta_by_lambda = -neg_theta/lambda;
        for (size_type j = 0; j < n; ++j)
        {
            s[j] += theta_by_lambda*ehat[j];
        }
        s.push_back(theta_by_lambda);
    }

    // Output solution
    copy(s.begin(), s.end(), s_first);
}

/**
 * Solve a Toeplitz set of linear equations in-place.  That is, compute
 * \f[
 *      L_{n+1}^{-1} d_{n+1}
 *      \mbox{ for }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    r_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$, \f$\vec{r}\f$, and \f$\vec{d}\f$.  The dimension of the
 * problem is fixed by <tt>n = distance(a_first, a_last)</tt>.  A symmetric
 * Toeplitz solve can be performed by having \f$\vec{a}\f$ and \f$\vec{r}\f$
 * iterate over the same data.  The Hermitian case requires two buffers with
 * \f$vec{r}\f$ being the conjugate of \f$\vec{a}\f$.  The working precision
 * is fixed by the \c value_type of \c d_first.
 *
 * @param[in]     a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]     a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in]     r_first Beginning of the range containing \f$\vec{r}\f$.
 * @param[in,out] d_first Beginning of the range containing \f$\vec{d}\f$.
 *                        Also the beginning of the output range to which
 *                        <strong><tt>n+1</tt></strong> entries will be
 *                        written.
 */
template<class RandomAccessIterator,
         class ForwardIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        RandomAccessIterator r_first,
                        ForwardIterator      d_first)
{
    return zohar_linear_solve(a_first, a_last, r_first, d_first, d_first);
}


/**
 * Solve a real-valued, symmetric Toeplitz set of linear equations in-place.
 * That is, compute
 * \f[
 *      L_{n+1}^{-1} d_{n+1}
 *      \mbox{ for }
 *      L_{n+1} = \bigl(\begin{smallmatrix}
 *                    1   & \tilde{a}_n \\
 *                    a_n & L_n
 *                \end{smallmatrix}\bigr)
 * \f]
 * given \f$\vec{a}\f$ and \f$\vec{d}\f$.  The dimension of the problem is
 * fixed by <tt>n = distance(a_first, a_last)</tt>.  The working precision is
 * fixed by the \c value_type of \c d_first.
 *
 * @param[in]     a_first Beginning of the range containing \f$\vec{a}\f$.
 * @param[in]     a_last  End of the range containing \f$\vec{a}\f$.
 * @param[in,out] d_first Beginning of the range containing \f$\vec{d}\f$.
 *                        Also the beginning of the output range to which
 *                        <strong><tt>n+1</tt></strong> entries will be
 *                        written.
 */
template<class RandomAccessIterator,
         class ForwardIterator>
void zohar_linear_solve(RandomAccessIterator a_first,
                        RandomAccessIterator a_last,
                        ForwardIterator      d_first)
{
    return zohar_linear_solve(a_first, a_last, a_first, d_first);
}

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Method-specific estimation variance routines following Broersen.
 *
 * For details see either the FiniteSampleCriteria.tex write up or Broersen,
 * P. M. T. "Finite sample criteria for autoregressive order selection." IEEE
 * Transactions on Signal Processing 48 (December 2000): 3550-3558.
 * http://dx.doi.org/10.1109/78.887047.
 *
 * The selection criteria routines might be sped up for floating point
 * arguments given an appropriate digamma (psi) or Pochhammer symbol
 * implementation.  To do so with the GNU Scientific Library (GSL), e.g., try
 * @code
 *     #include <gsl/gsl_sf_psi.h>
 *     #define AR_DIGAMMA(x) gsl_sf_psi(x)
 *
 *     #include <gsl/gsl_sf_gamma.h>
 *     #define AR_POCHHAMMER(a,x) gsl_sf_poch(a,x)
 * @endcode
 * before including this header and link the GSL with your binary.
 * @{
 */

/** Denotes the sample mean was subtracted from a signal before estimation. */
struct mean_subtracted
{
    /** Computes the empirical variance estimate for order zero. */
    template <typename Result, typename Integer>
    static Result empirical_variance_zero(Integer N)
    {
        assert(N >= 1);

        Result den = N;
        return 1 / den;
    }
};

/** Denotes the sample mean was retained in a signal during estimation. */
struct mean_retained
{
    /** Computes the empirical variance estimate for order zero. */
    template <typename Result, typename Integer>
    static Result empirical_variance_zero(Integer)
    {
        return 0;
    }
};

namespace { // anonymous

// Oh the simultaneous love/hate because -Werror => -Werror=type-limits...

template <typename T, bool> struct is_nonnegative_helper;

template <typename T> struct is_nonnegative_helper<T, /*signed?*/ true>
{
    static bool check(const T& t) { return t >= 0; }
};

template <typename T> struct is_nonnegative_helper<T, /*signed?*/ false>
{
    static bool check(const T&) { return true; }
};

template <typename T> bool is_nonnegative(const T& t)
{
    using std::numeric_limits;
    return is_nonnegative_helper<T,numeric_limits<T>::is_signed>::check(t);
}

}


/**
 * A parent type for autoregressive process parameter estimation techniques.
 *
 * Each subclass should have an <tt>empirical_variance(N, i)</tt> method
 * following Broersen, P. M. and H. E. Wensink. "On Finite Sample Theory for
 * Autoregressive Model Order Selection." IEEE Transactions on Signal
 * Processing 41 (January 1993): 194+.
 * http://dx.doi.org/10.1109/TSP.1993.193138.
 */
struct estimation_method {};

/** Represents estimation by solving the Yule--Walker equations. */
template <class MeanHandling>
class YuleWalker : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(is_nonnegative(i));
        assert(static_cast<Integer1>(i) <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        Result num = N - i, den = N*(N + 2);
        return num / den;
    }
};

/** Represents estimation using %Burg's recursive method. */
template <class MeanHandling>
class Burg : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(is_nonnegative(i));
        assert(static_cast<Integer1>(i) <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        Result den = N + 1 - i;
        return 1 / den;
    }
};

/** Represents forward and backward prediction least squares minimization. */
template <class MeanHandling>
class LSFB : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(is_nonnegative(i));
        assert(static_cast<Integer1>(i) <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        // Factorizing expression will cause problems in unsigned arithmetic
        Result den = N + Result(3)/2 - Result(3)/2 * i;
        return 1 / den;
    }
};

/** Represents forward prediction least squares minimization. */
template <class MeanHandling>
class LSF : public estimation_method
{
public:

    /**
     * Approximates the empirical variance estimate.
     * @param N Number of observations.
     * @param i Variance order.
     */
    template <typename Result, typename Integer1, typename Integer2>
    static Result empirical_variance(Integer1 N, Integer2 i)
    {
        assert(N >= 1);
        assert(is_nonnegative(i));
        assert(static_cast<Integer1>(i) <= N);

        if (i == 0)
            return MeanHandling::template empirical_variance_zero<Result>(N);

        // Factorizing expression will cause problems in unsigned arithmetic
        Result den = N + 2 - 2*i;
        return 1 / den;
    }
};

/** An STL-ready binary_function for a given method's empirical variance. */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
struct empirical_variance_function
    : public std::binary_function<Result,Integer1,Integer2>
{
    Result operator() (Integer1 N, Integer2 i) const
    {
        return EstimationMethod::template empirical_variance<Result>(N, i);
    }
};

/**
 * An STL AdaptableGenerator for a given method's empirical variance.
 *
 * On each <tt>operator()</tt> invocation the method's empirical variance is
 * returned for the current model order.  The first invocation returns the
 * result for model order zero.
 */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
class empirical_variance_generator
{
private:

    Integer1 N;
    Integer2 i;

public:

    typedef Result result_type;

    empirical_variance_generator(Integer1 N) : N(N), i(0) {}

    Result operator() ()
    {
        return EstimationMethod::template empirical_variance<Result>(N, i++);
    }
};

/**
 * An immutable RandomAccessIterator over a method's empirical variance
 * sequence.
 *
 * Facilitates using algorithms like <tt>std::copy</tt>,
 * <tt>std::accumulate</tt>, and <tt>std::partial_sum</tt> when comparing a
 * hierarchy of models during model order selection.
 *
 * The (N+1)-length sequence of orders 0, 1, ..., N is iterated given sample
 * size N.  Default constructed instances represent past-end iterators.
 */
template <class EstimationMethod,
          typename Result,
          typename Integer1 = std::size_t,
          typename Integer2 = Integer1>
class empirical_variance_iterator
    : public std::iterator<std::random_access_iterator_tag, Result,
                           std::ptrdiff_t, const Result*, const Result&>
{
private:
    typedef std::iterator<std::random_access_iterator_tag, Result,
                          std::ptrdiff_t, const Result*, const Result&> base;

    Integer1 N;
    Integer2 i;

public:
    typedef typename base::difference_type   difference_type;
    typedef typename base::iterator_category iterator_category;
    typedef typename base::pointer           pointer;
    typedef typename base::reference         reference;
    typedef typename base::value_type        value_type;

    /** Construct a past-end iterator */
    empirical_variance_iterator() : N(0), i(1) {}  // Null

    /** Construct an iterator over sequence order 0, 1, ..., N (inclusive). */
    empirical_variance_iterator(Integer1 N) : N(N), i(0)
        { assert(N >= 1); }

    /** Construct an iterator over sequence order i, i+1, ..., N (inclusive). */
    empirical_variance_iterator(Integer1 N, Integer2 i) : N(N), i(i) {}

    // Forward traversal

    empirical_variance_iterator& operator++()
        { assert(N >= 1); ++i; return *this; }

    empirical_variance_iterator operator++(int)
        { assert(N >= 1); return empirical_variance_iterator(N, i++); }

    empirical_variance_iterator operator+(const difference_type& k)
        { assert(N >= 1); return empirical_variance_iterator(N, i + k); }

    empirical_variance_iterator& operator+=(const difference_type& k)
        { assert(N >= 1); i += k; return *this; }

    // Backward traversal

    empirical_variance_iterator& operator--()
        { assert(N >= 1); --i; return *this; }

    empirical_variance_iterator operator--(int)
        { assert(N >= 1); return empirical_variance_iterator(N, i--); }

    empirical_variance_iterator operator-(const difference_type& k)
        { assert(N >= 1); return empirical_variance_iterator(N, i - k); }

    empirical_variance_iterator& operator-=(const difference_type& k)
        { assert(N >= 1); i -= k; return *this; }

    // Distance support

    difference_type operator-(const empirical_variance_iterator& other) const
    {
        if (!this->N) {
            return 1 + other.N - other.i;
        } else if (!other.N) {
            return -static_cast<difference_type>(1 + this->N - this->i);
        } else {
            assert(this->N == other.N);
            return this->i - other.i;
        }
    }

    // EqualityComparable

    bool operator==(const empirical_variance_iterator& other) const
    {
        if (!this->N) {
            return other.i >= static_cast<Integer2>(other.N + 1);
        } else if (!other.N) {
            return this->i >= static_cast<Integer2>(this->N + 1);
        } else {
            return this->N == other.N && this->i == other.i;
        }
    }

    bool operator!=(const empirical_variance_iterator& other) const
        { return !(*this == other); }

    // LessThanComparable will trigger assertion on nonsense N cases

    bool operator<(const empirical_variance_iterator& other) const
        { return (*this - other) < 0; }

    bool operator<=(const empirical_variance_iterator& other) const
        { return (*this - other) <= 0; }

    bool operator>(const empirical_variance_iterator& other) const
        { return (*this - other) > 0; }

    bool operator>=(const empirical_variance_iterator& other) const
        { return (*this - other) >= 0; }

    // Dereference operations

    const value_type operator*() const
    {
        assert(is_nonnegative(i));
        assert(i <= static_cast<Integer2>(N));

        return EstimationMethod::template empirical_variance<Result>(N, i);
    }

    const value_type operator[](const difference_type &k) const
    {
        assert(is_nonnegative(i + k));
        assert(i + k <= static_cast<Integer2>(N));

        return EstimationMethod::template empirical_variance<Result>(N, i + k);
    }
};

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Criteria for autoregressive model order selection following Broersen.
 *
 * For details see either the FiniteSampleCriteria.tex write up or Broersen,
 * P. M. T. "Finite sample criteria for autoregressive order selection." IEEE
 * Transactions on Signal Processing 48 (December 2000): 3550-3558.
 * http://dx.doi.org/10.1109/78.887047.
 *
 * @{
 */

/**
 * A parent type for autoregressive model selection criterion.
 *
 * Each subclass should have an <tt>overfit_penalty(N, p)</tt> method
 * following Broersen, P. M. and H. E. Wensink. "On Finite Sample Theory for
 * Autoregressive Model Order Selection." IEEE Transactions on Signal
 * Processing 41 (January 1993): 194+.
 * http://dx.doi.org/10.1109/TSP.1993.193138.
 */
struct criterion
{
    /** Compute the underfit penalty given \f$\sigma^2_\epsilon\f$. */
    template <typename Result, typename Input>
    static Result underfit_penalty(Input sigma2e)
    {
        using std::log;
        return log(Result(sigma2e));
    }
};

/**
 * Evaluate a given \ref criterion for \c N samples and model order \c p.
 *
 * @param[in] sigma2e The residual \f$\sigma^2_\epsilon\f$
 * @param[in] N       Sample count used to compute \f$\sigma^2_\epsilon\f$
 * @param[in] p       The model order use to compute \f$sigma^2_\epsilon\f$
 *
 * @return the evaluated criterion.
 */
template <class    Criterion,
          typename Result,
          typename Integer1,
          typename Integer2>
Result evaluate(Result sigma2e, Integer1 N, Integer2 p)
{
    Result underfit = Criterion::template underfit_penalty<Result>(sigma2e);
    Result overfit  = Criterion::template overfit_penalty<Result>(N, p);
    return underfit + overfit;
}

/**
 * Represents the generalized information criterion (GIC).  The penalty factor
 * \f$\alpha\f$ is controlled by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <int AlphaNumerator = 3, int AlphaDenominator = 1>
struct GIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        return Result(AlphaNumerator) * p / (N * AlphaDenominator);
    }
};

/** Represents the Akaike information criterion (AIC). */
struct AIC : public GIC<2> {};

/** Represents the consistent criterion BIC. */
struct BIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::log;
        return log(Result(N)) * p / N;
    }
};

/** Represents the minimally consistent criterion (MCC). */
struct MCC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::log;
        return 2*log(log(Result(N))) * p / N;
    }
};

/**
 * Represents the asymptotically-corrected Akaike information criterion (AICC).
 */
struct AICC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        return 2 * Result(p) / (N - p - 1);
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to a particular
 * \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled by
 * <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class EstimationMethod,
          int AlphaNumerator = 3,
          int AlphaDenominator = 1 >
struct FIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        // This accumulate invocation is inefficient but always correct
        using std::accumulate;
        typedef empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
            > evi_type;
        Result sum = accumulate(evi_type(N), evi_type(N) + p, Result(0));

        return AlphaNumerator * sum / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * YuleWalker \ref estimation_method.  The penalty factor \f$\alpha\f$ is
 * controlled by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<YuleWalker<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = YuleWalker<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));

        Result num = (2*N - p - 1)*p;  // Avoids non-positive values
        Result den = 2*N*(N + 2);

        return AlphaNumerator * (v0 + num/den) / AlphaDenominator;
    }
};

// Specializations of the FIC for efficiency when digamma is available.
#ifdef AR_DIGAMMA

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * Burg \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<Burg<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = Burg<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        t -= AR_DIGAMMA(N + 1);
        t += AR_DIGAMMA(N + 1 - p);
        return AlphaNumerator * t / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * LSFB \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<LSFB<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = LSFB<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_DIGAMMA((3 + 2*N) / Result(3)    );
        Result b = AR_DIGAMMA((3 + 2*N) / Result(3) - p);
        return AlphaNumerator * (t - 2*(a - b)/3) / AlphaDenominator;
    }
};

/**
 * Represents the finite information criterion (FIC) as applied to the \ref
 * LSF \ref estimation_method.  The penalty factor \f$\alpha\f$ is controlled
 * by <tt>AlphaNumerator / AlphaDenominator</tt>.
 */
template <class MeanHandling,
          int AlphaNumerator,
          int AlphaDenominator>
struct FIC<LSF<MeanHandling>, AlphaNumerator, AlphaDenominator>
    : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result t  = LSF<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_DIGAMMA((2 + N) / Result(2)    );
        Result b = AR_DIGAMMA((2 + N) / Result(2) - p);
        return AlphaNumerator * (t - (a - b)/2) / AlphaDenominator;
    }
};

#endif /* AR_DIGAMMA */

/**
 * Represents the finite sample information criterion (FSIC) as applied to a
 * particular \ref estimation_method.
 */
template <class EstimationMethod>
struct FSIC : public criterion
{

private:

    /** A helper to compute \f$\frac{1+v}{1-v}\f$ for the method in use. */
    template <typename Result, typename Integer1, typename Integer2>
    class product_iterator
        : public empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
          >
    {
    private:
        typedef empirical_variance_iterator<
                EstimationMethod, Result, Integer1, Integer2
            > base;

    public:
        typedef typename base::difference_type   difference_type;
        typedef typename base::iterator_category iterator_category;
        typedef typename base::pointer           pointer;
        typedef typename base::reference         reference;
        typedef typename base::value_type        value_type;

        product_iterator(Integer1 N) : base(N) {}

        value_type operator*() const
        {
            const value_type v = this->base::operator*();
            return (1 + v) / (1 - v);
        }

        value_type operator[](const difference_type &k) const
        {
            const value_type v = this->base::operator[](k);
            return (1 + v) / (1 - v);
        }

    };

public:

    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        // This accumulate invocation is inefficient but always correct
        using std::multiplies;
        using std::accumulate;
        product_iterator<Result, Integer1, Integer2> first(N), last(N);
        last += p; // Avoids overloading requirements for (last + p)

        return accumulate(first, last, Result(1), multiplies<Result>()) - 1;
    }
};

// Specializations of the FIC for efficiency when Pochhammer is available.
#ifdef AR_POCHHAMMER

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * YuleWalker \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<YuleWalker<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = YuleWalker<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_POCHHAMMER(N*(N+3) - p, p);
        Result b = AR_POCHHAMMER(Result(1 + N) - N*N);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * Burg \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<Burg<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = Burg<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_POCHHAMMER(-Result(1) - N, p);
        Result b = AR_POCHHAMMER( Result(1) - N, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * LSFB \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<LSFB<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = LSFB<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_POCHHAMMER((-Result(2)-2*N)/3, p);
        Result b = AR_POCHHAMMER(( Result(2)-2*N)/3, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

/**
 * Represents the finite information criterion (FSIC) as applied to the \ref
 * LSF \ref estimation_method.
 */
template <class MeanHandling>
struct FSIC<LSF<MeanHandling> > : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        Result v0  = LSF<MeanHandling>
            ::template empirical_variance<Result>(N, Integer2(0));
        Result a = AR_POCHHAMMER((-Result(1)-N)/2, p);
        Result b = AR_POCHHAMMER(( Result(1)-N)/2, p);

        return (1 + v0) / (1 - v0) * (a / b) - 1;
    }
};

#endif /* AR_POCHHAMMER */


/**
 * Represents the combined information criterion (CIC) as applied to a
 * particular \ref estimation_method.
 */
template <class EstimationMethod>
struct CIC : public criterion
{
    /** Compute overfit penalty given \c N observations at model order \c p. */
    template <typename Result, typename Integer1, typename Integer2>
    static Result overfit_penalty(Integer1 N, Integer2 p)
    {
        using std::max;
        return max(
                FSIC<EstimationMethod>::template overfit_penalty<Result>(N, p),
                FIC <EstimationMethod>::template overfit_penalty<Result>(N, p)
            );
    }
};

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Algorithmic helpers for autoregressive model order selection.
 *
 * @{
 */

/**
 * Evaluate a \ref criterion on a hierarchy of models given
 * \f$\sigma^2_\epsilon\f$ for each model.  The index of the best model, i.e.
 * the one with minimum criterion value, is returned.
 *
 * @param[in]  N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 * @param[in]  ordfirst The model order corresponding to \c first.
 *                      When \f$sigma^2_\epsilon\f$ is produced entirely by
 *                      \ref burg_method with <tt>hierarchy == true</tt>,
 *                      this should be \c 0.
 * @param[in]  first    Beginning of the range holding \f$\sigma^2_\epsilon\f$
 * @param[in]  last     Exclusive end of input range.
 * @param[out] crit     Value assigned to each model by the criterion.
 *
 * @return The distance from \c first to the best model.
 *         An obscenely negative value is returned on error.
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    InputIterator,
          class    OutputIterator>
typename std::iterator_traits<InputIterator>::difference_type
evaluate_models(Integer1       N,
                Integer2       ordfirst,
                InputIterator  first,
                InputIterator  last,
                OutputIterator crit)
{
    using std::iterator_traits;
    using std::numeric_limits;

    typedef InputIterator iterator;
    typedef typename iterator_traits<iterator>::difference_type difference_type;
    typedef typename iterator_traits<iterator>::value_type      value_type;

    // Short circuit on trivial input
    if (first == last)
        return numeric_limits<difference_type>::min();

    // Handle first iteration without comparison as AICC blows up on N == 1
    value_type best_val = evaluate<Criterion>(*first++, N, ordfirst);
    difference_type best_pos = 0, dist = 0;

    // Scan through rest of candidates (up to order N-1) updating best as we go
    while (first != last && static_cast<Integer1>(++dist) < N)
    {
        value_type candidate = evaluate<Criterion>(*first++, N, dist+ordfirst);
        *crit++ = candidate;

        if (candidate < best_val) {
            best_val = candidate;
            best_pos = dist;
        }
    }

    return best_pos;
}

/**
 * Obtain the best model according to \ref criterion applied to
 * \f$\sigma^2_\epsilon\f$ given a hierarchy of candidates.
 *
 * On input, \c params, \c sigma2e, \c gain, and \c autocor should be <a
 * href="http://www.sgi.com/tech/stl/Sequence.html">Sequence</a>s which were
 * populated by \ref burg_method when \c hierarchy is \c true (or in some other
 * equivalent manner).  On output, these arguments will contain only values
 * relevant to the best model.
 *
 * @param[in]     N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 *                         For example, the return value of \ref burg_method.
 * @param[in]     minorder Constrain the best model to be at least this order.
 *                         Supplying zero specifies no constraint.
 * @param[in,out] params   Model parameters
 * @param[in,out] sigma2e  \f$\sigma^2_\epsilon\f$
 * @param[in,out] gain     Model gain
 * @param[in,out] autocor  Model autocorrelations
 * @param[out]    crit     Value assigned to each model by the criterion.
 *
 * @return The index of the best criterion value within \c crit.
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    Sequence1,
          class    Sequence2,
          class    Sequence3,
          class    Sequence4,
          class    OutputIterator>
typename Sequence1::difference_type
best_model(Integer1       N,
           Integer2       minorder,
           Sequence1&     params,
           Sequence2&     sigma2e,
           Sequence3&     gain,
           Sequence4&     autocor,
           OutputIterator crit)
{
    using std::advance;
    using std::copy_backward;
    using std::invalid_argument;
    using std::size_t;

    // Ensure all inputs have conformant sizes
    // Sanity checks written to provide order messages readable by end users
    AR_ENSURE_ARG(sigma2e.size() > 0);
    const size_t maxorder = sigma2e.size() - 1;
    AR_ENSURE_ARG(maxorder > 0);
    AR_ENSURE_ARG(is_nonnegative(minorder));
    AR_ENSURE_ARG(static_cast<size_t>(minorder) <= maxorder);
    AR_ENSURE_ARG(params .size() == (sigma2e.size()-1)*sigma2e.size()/2);
    AR_ENSURE_ARG(gain   .size() == sigma2e.size());
    AR_ENSURE_ARG(autocor.size() == sigma2e.size());

    // Find the best model index according to the given minorder and criterion
    typename Sequence1::difference_type best = minorder;
    {
        typename Sequence2::iterator iter = sigma2e.begin();
        advance(iter, minorder);
        best += evaluate_models<Criterion>(N, minorder, iter,
                                           sigma2e.end(), crit);
    }

    // Now trim away everything and leaving only the best model in Sequences...

    // ...first in params...
    {
        // Best parameters might overlap beginning of params so copy_backwards.
        // AR(0) is trivial and AR(1) starts at params.begin(), hence off by 1.
        typename Sequence1::iterator first  = params.begin();
        typename Sequence1::iterator last   = params.begin();
        typename Sequence1::iterator result = params.begin();
        advance(first,  ((best-1)*best)/2       );
        advance(last,   ((best-1)*best)/2 + best);
        advance(result,                     best);
        copy_backward(first, last, result);
        params.resize(best);
    }

    // ...next in sigma2e...
    {
        typename Sequence2::iterator cursor = sigma2e.begin();
        advance(cursor, best);
        *sigma2e.begin() = *cursor;
        sigma2e.resize(1);
    }

    // ...next in gain...
    {
        typename Sequence2::iterator cursor = gain.begin();
        advance(cursor, best);
        *gain.begin() = *cursor;
        gain.resize(1);
    }

    // ...and last in autocor...
    autocor.resize(best + 1);

    // ...but notice that resizing does not free capacity for std::vectors.

    return best;
}

namespace { // anonymous

// Used to discard output per http://stackoverflow.com/questions/335930/
struct null_output : std::iterator< std::output_iterator_tag, null_output >
{
    template <typename T> void operator=(const T&) {}

    null_output& operator++()
        { return *this; }

    null_output operator++(int)
        { null_output it(*this); ++*this; return it; }

    null_output& operator*()
        { return *this; }
};

}

/**
 * Find the index of the best model from a hierarchy of candidates according to
 * a \ref criterion given \f$\sigma^2_\epsilon\f$ for each model.
 *
 * @param[in]  N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 * @param[in]  ordfirst The model order corresponding to \c first.
 *                      When \f$sigma^2_\epsilon\f$ is produced entirely by
 *                      \ref burg_method, this should be \c 0u.
 * @param[in]  first    Beginning of the range holding \f$\sigma^2_\epsilon\f$
 * @param[in]  last     Exclusive end of input range.
 *
 * @return The distance from \c first to the best model.
 *
 * @see evaluate_models(Criterion,Integer1,Integer2,InputIterator,OutputIterator)
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    InputIterator>
typename std::iterator_traits<InputIterator>::difference_type
evaluate_models(Integer1      N,
                Integer2      ordfirst,
                InputIterator first,
                InputIterator last)
{
    return evaluate_models<Criterion>(N, ordfirst, first, last, null_output());
}

/**
 * Obtain the best model according to \ref criterion applied to
 * \f$\sigma^2_\epsilon\f$ given a hierarchy of candidates.
 *
 * On input, \c params, \c sigma2e, \c gain, and \c autocor should be <a
 * href="http://www.sgi.com/tech/stl/Sequence.html">Sequence</a>s which were
 * populated by \ref burg_method when \c hierarchy is \c true (or in some other
 * equivalent manner).  On output, these arguments will contain only values
 * relevant to the best model.
 *
 * @param[in]     N        Sample count used to compute \f$\sigma^2_\epsilon\f$.
 *                         For example, the return value of \ref burg_method.
 * @param[in]     minorder Constrain the best model to be at least this order.
 *                         Supplying zero specifies no constraint.
 * @param[in,out] params   Model parameters
 * @param[in,out] sigma2e  \f$\sigma^2_\epsilon\f$
 * @param[in,out] gain     Model gain
 * @param[in,out] autocor  Model autocorrelations
 *
 * @return The index of the best model within the inputs.
 *
 * @see best_model(Integer1,Integer2,Sequence1,Sequence2,Sequence3,Sequence4,OutputIterator)
 */
template <class    Criterion,
          typename Integer1,
          typename Integer2,
          class    Sequence1,
          class    Sequence2,
          class    Sequence3,
          class    Sequence4>
typename Sequence1::difference_type
best_model(Integer1       N,
           Integer2       minorder,
           Sequence1&     params,
           Sequence2&     sigma2e,
           Sequence3&     gain,
           Sequence4&     autocor)
{
    return best_model<Criterion>(N, minorder, params, sigma2e, gain,
                                 autocor, null_output());
}

/**
 * A template typedef and helper method returning a \ref best_model
 * implementation matching a model selection criterion provided at runtime.
 * Intended for use within interactive APIs, this method handles much of the
 * ugliness of resolving overloaded criterion logic and accounting for which
 * methods care as to whether or not the sample mean has been subtracted.
 * Class template parameters should be set to match the desired overload
 * returned by the \ref lookup method.
 *
 * Methods are looked up using their abbreviations as follows:
 * <dl>
 * <dt>AIC </dt><dd>Akaike information criterion</dd>
 * <dt>AICC</dt><dd>asymptotically-corrected Akaike information criterion</dd>
 * <dt>BIC </dt><dd>consistent criterion BIC</dd>
 * <dt>CIC </dt><dd>combined information criterion</dd>
 * <dt>FIC </dt><dd>finite information criterion</dd>
 * <dt>FSIC</dt><dd>finite sample information criterion</dd>
 * <dt>GIC </dt><dd>generalized information criterion</dd>
 * <dt>MCC </dt><dd>minimally consistent criterion</dd>
 * </dl>
 * Leading or trailing whitespace as well as capitalization differences are
 * ignored.
 *
 * @tparam EstimationMethod One of Burg, YuleWalker, LSFB, or LSF.
 * @tparam Integer1         Per the \c Integer1 parameter of \ref best_model
 * @tparam Integer2         Per the \c Integer2 parameter of \ref best_model
 * @tparam Sequence1        Per the \c Sequence1 parameter of \ref best_model
 * @tparam Sequence2        Per the \c Sequence2 parameter of \ref best_model
 * @tparam Sequence3        Per the \c Sequence3 parameter of \ref best_model
 * @tparam Sequence4        Per the \c Sequence4 parameter of \ref best_model
 */
template <template <class> class EstimationMethod,
          typename Integer1,
          typename Integer2,
          class    Sequence1,
          class    Sequence2 = Sequence1,
          class    Sequence3 = Sequence2,
          class    Sequence4 = Sequence3>
struct best_model_function
{
    /** The type returned by the \ref lookup function. */
    typedef typename Sequence1::difference_type (*type)(Integer1   N,
                                                        Integer2   minorder,
                                                        Sequence1& params,
                                                        Sequence2& sigma2e,
                                                        Sequence3& gain,
                                                        Sequence4& autocor);

    /**
     * Lookup a \ref best_model function pointer matching \c EstimationMethod,
     * the specified criterion c abbreviation, and whether or not the sample
     * mean has been subtracted from the data.
     *
     * @param abbreviation  Known abbreviations per \ref best_model_function.
     *                      If only whitespace, a reasonable default is used.
     * @param subtract_mean Has the sample mean been subtracted from the data?
     *
     * @tparam CharT     Permits use by any \c std::basic_string instantiation
     * @tparam Traits    Permits use by any \c std::basic_string instantiation
     * @tparam Allocator Permits use by any \c std::basic_string instantiation
     *
     * @return A function pointer which, when invoked, calls \ref best_model.
     *         When no known criterion matches \c abbreviation, \c NULL
     *         is returned.
     */
    template <class CharT, class Traits, class Allocator>
    static type lookup(std::basic_string<CharT,Traits,Allocator> abbreviation,
                       const bool subtract_mean);

    /**
     * Lookup a \ref best_model function pointer matching \c EstimationMethod,
     * the specified criterion c abbreviation, and whether or not the sample
     * mean has been subtracted from the data.
     *
     * @param abbreviation  Known abbreviations per \ref best_model_function.
     *                      If only whitespace, a reasonable default is used.
     * @param subtract_mean Has the sample mean been subtracted from the data?
     *
     * @return A function pointer which, when invoked, calls \ref best_model.
     *         When no known criterion matches \c abbreviation, \c NULL
     *         is returned.
     */
    static type lookup(const char *abbreviation, const bool subtract_mean)
    {
        std::string s(abbreviation);
        return lookup(s, subtract_mean);
    }
};

template <template <class> class EstimationMethod,
          typename Integer1,
          typename Integer2,
          class    Sequence1,
          class    Sequence2,
          class    Sequence3,
          class    Sequence4>
template <class CharT, class Traits, class Allocator>
typename best_model_function<
    EstimationMethod,Integer1,Integer2,Sequence1,Sequence2,Sequence3,Sequence4
>::type
best_model_function<
    EstimationMethod,Integer1,Integer2,Sequence1,Sequence2,Sequence3,Sequence4
>::lookup(std::basic_string<CharT,Traits,Allocator> abbrev,
          const bool subtract_mean)
{
    using std::basic_string;
    using std::toupper;

    // Canonicalize the abbreviation by making it uppercase and trimming it
    // For nothing but this reason the method accepts 'abbrev' by value
    typedef basic_string<CharT,Traits,Allocator> string_type;
    for (typename string_type::iterator p = abbrev.begin();
         abbrev.end() != p; ++p)
    {
        *p = toupper(*p);
    }
    abbrev.erase(0, abbrev.find_first_not_of(" \n\r\t"));
    abbrev.erase(1 + abbrev.find_last_not_of(" \n\r\t"));

    // Obtain best_model(...) per abbrev, subtract_mean, EstimationMethod
    // and assign to retval to provide unambiguous overload resolution.
    type retval;

    if      (abbrev.empty() || 0 == abbrev.compare("CIC" ))  // Default
    {
        if (subtract_mean)
            retval = best_model<CIC<EstimationMethod<mean_subtracted> > >;
        else
            retval = best_model<CIC<EstimationMethod<mean_retained  > > >;
    }
    else if (0 == abbrev.compare("AIC" ))
    {
        retval = best_model<AIC>;
    }
    else if (0 == abbrev.compare("AICC"))
    {
        retval = best_model<AICC>;
    }
    else if (0 == abbrev.compare("BIC" ))
    {
        retval = best_model<BIC>;
    }
    else if (0 == abbrev.compare("FIC" ))
    {
        if (subtract_mean)
            retval = best_model<FIC<EstimationMethod<mean_subtracted> > >;
        else
            retval = best_model<FIC<EstimationMethod<mean_retained  > > >;
    }
    else if (0 == abbrev.compare("FSIC"))
    {
        if (subtract_mean)
            retval = best_model<FSIC<EstimationMethod<mean_subtracted> > >;
        else
            retval = best_model<FSIC<EstimationMethod<mean_retained  > > >;
    }
    else if (0 == abbrev.compare("GIC" ))
    {
        retval = best_model<GIC<> >;
    }
    else if (0 == abbrev.compare("MCC" ))
    {
        retval = best_model<MCC>;
    } else
    {
        retval = NULL;
    }

    return retval;
}

/**
 * An adapter to add striding over another (usually random access) iterator.
 *
 * Modified from "C++ Cookbook" by Stephens, Diggins, Turkanis, and Cogswell.
 * Copyright 2006 O'Reilly Media, Inc., ISBN 0-596-00761-2.
 */
template<class Iterator>
class strided_adaptor
{
public:
    typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
    typedef typename std::iterator_traits<Iterator>::value_type        value_type;
    typedef typename std::iterator_traits<Iterator>::difference_type   difference_type;
    typedef typename std::iterator_traits<Iterator>::pointer           pointer;
    typedef typename std::iterator_traits<Iterator>::reference         reference;

    strided_adaptor() : i(NULL), step(0) {};

    strided_adaptor(const strided_adaptor& o) : i(o.i), step(o.step) {}

    strided_adaptor(Iterator x, difference_type n) : i(x), step(n) {}

    strided_adaptor& operator++()
    {
        using std::advance;
        advance(i, step);
        return *this;
    }

    strided_adaptor operator++(int)
    {
        strided_adaptor tmp(*this);
        using std::advance;
        advance(i, step);
        return tmp;
    }

    strided_adaptor& operator+=(difference_type x)
    {
        using std::advance;
        advance(i, x*step);
        return *this;
    }

    strided_adaptor& operator--( )
    {
        using std::advance;
        advance(i, -step);
        return *this;
    }

    strided_adaptor operator--(int)
    {
        strided_adaptor tmp(*this);
        using std::advance;
        advance(i, -step);
        return tmp;
    }

    strided_adaptor& operator-=(difference_type x)
    {
        using std::advance;
        advance(i, -x*step);
        return *this;
    }

    reference operator[](difference_type n)
    {
        return i[n * step];
    }

    reference operator*()
    {
        return *i;
    }

    bool operator==(const strided_adaptor& o) const
    {
        return i == o.i;
    }

    bool operator!=(const strided_adaptor& o) const
    {
        return i != o.i;
    }

    bool operator<(const strided_adaptor& o) const
    {
        return i < o.i;
    }

    difference_type operator-(const strided_adaptor& o) const
    {
        return (i - o.i) / step;
    }

    strided_adaptor operator+(difference_type x) const
    {
        strided_adaptor tmp(*this);
        tmp += x;
        return tmp;
    }

private:
    Iterator i;
    difference_type step;
};

/**
 * @}
 */

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

} // end namespace ar

#endif /* AR_HPP */
