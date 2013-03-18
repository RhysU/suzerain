//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_RUNNING_STATISTICS_HPP
#define SUZERAIN_RUNNING_STATISTICS_HPP

/** @file
 * Provides statistical accumulators for special cases of interest.
 * Though Boost Accumulators is generally preferred, some use cases can
 * benefit from custom logic.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace suzerain {

// TODO Assertions in running_statistics complain loudly when N == 0

/**
 * Accumulates running minimum, maximum, mean, and variance details from a data
 * stream.  Adapted from http://www.johndcook.com/standard_deviation.html.
 * Extended to track a fixed number of quantities with concurrently provided
 * samples.
 *
 * @tparam Real Floating point type used for input and accumulation
 * @tparam N    Number of distinct quantities to simultaneously track
 */
template <typename Real, std::size_t N>
class running_statistics
{
public:

    /** Default constructor. */
    running_statistics();

    /** Provide quantity samples in locations <tt>x[0], ..., x[N-1]</tt>. */
    void operator()(const Real x[N]);

    /** Obtain the running number of samples provided thus far. */
    inline std::size_t count() const;

    /** Obtain the minimum sample observed for quantity \c i. */
    inline Real min(std::size_t i) const;

    /** Obtain the maximum sample observed for quantity \c i. */
    inline Real max(std::size_t i) const;

    /** Obtain the running mean for quantity \c i. */
    inline Real avg(std::size_t i) const;

    /** Obtain the running variance for quantity \c i. */
    inline Real var(std::size_t i) const;

    /** Obtain the running standard deviation for quantity \c i. */
    inline Real std(std::size_t i) const;

    /** Reset the instance to its newly constructed state. */
    void clear();

private:

    Real oldM_[N], newM_[N], oldS_[N], newS_[N], min_[N], max_[N];

    std::size_t n_;
};

template <typename Real, std::size_t N>
running_statistics<Real,N>::running_statistics()
{
    clear();
}

template <typename Real, std::size_t N>
void running_statistics<Real,N>::operator()(const Real x[N])
{
    // Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232
    using std::min;
    using std::max;
    ++n_;

    if (n_ > 1) {  // Second and subsequent iteration

        // Process moments for current iteration
        const Real inv_n = static_cast<Real>(1) / n_;
        for (std::size_t i = 0; i < N; ++i)
            newM_[i] = oldM_[i] + (x[i] - oldM_[i]) * inv_n;
        for (std::size_t i = 0; i < N; ++i)
            newS_[i] = oldS_[i] + (x[i] - oldM_[i])*(x[i] - newM_[i]);

        // Prepare to process moments for next iteration
        for (std::size_t i = 0; i < N; ++i) oldM_[i] = newM_[i];
        for (std::size_t i = 0; i < N; ++i) oldS_[i] = newS_[i];

        // Process extrema for current iteration
        for (std::size_t i = 0; i < N; ++i) min_[i] = min(min_[i], x[i]);
        for (std::size_t i = 0; i < N; ++i) max_[i] = max(max_[i], x[i]);

    } else {       // Initial iteration

        // Initialize moments
        for (std::size_t i = 0; i < N; ++i) oldM_[i] = x[i];
        for (std::size_t i = 0; i < N; ++i) newM_[i] = x[i];
        for (std::size_t i = 0; i < N; ++i) oldS_[i] = 0;

        // Initialize extrema
        for (std::size_t i = 0; i < N; ++i) min_[i] = x[i];
        for (std::size_t i = 0; i < N; ++i) max_[i] = x[i];

    }
}

template <typename Real, std::size_t N>
std::size_t running_statistics<Real,N>::count() const
{
    return n_;
}

template <typename Real, std::size_t N>
Real running_statistics<Real,N>::min(std::size_t i) const
{
    assert(i < N);
    return min_[i];
}

template <typename Real, std::size_t N>
Real running_statistics<Real,N>::max(std::size_t i) const
{
    assert(i < N);
    return max_[i];
}

template <typename Real, std::size_t N>
Real running_statistics<Real,N>::avg(std::size_t i) const
{
    assert(i < N);
    using std::numeric_limits;
    return (n_ > 0) ? newM_[i] : numeric_limits<Real>::quiet_NaN();
}

template <typename Real, std::size_t N>
Real running_statistics<Real,N>::var(std::size_t i) const
{
    assert(i < N);
    using std::numeric_limits;
    return  (n_ >  1) ? newS_[i]/(n_ - 1)
          : (n_ == 1) ? 0
          :             numeric_limits<Real>::quiet_NaN();
}

template <typename Real, std::size_t N>
Real running_statistics<Real,N>::std(std::size_t i) const
{
    assert(i < N);
    using std::sqrt;
    return sqrt(var(i));
}

template <typename Real, std::size_t N>
void running_statistics<Real,N>::clear()
{
    n_ = 0;
    using std::numeric_limits;
    for (std::size_t i = 0; i < N; ++i)
        min_[i] = max_[i] = numeric_limits<Real>::quiet_NaN();
}

} // namespace suzerain

#endif // SUZERAIN_RUNNING_STATISTICS_HPP
