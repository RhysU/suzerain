//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
// TODO Permit serialization and reconstruction using a raw buffer
// TODO Provide MPI reduction operator using serialization capability
// TODO Defend against overflowing the counter, especially for MPI reduction

/**
 * Accumulates running minimum, maximum, mean, and variance details from a data
 * stream.  Adapted from http://www.johndcook.com/standard_deviation.html.
 * Extended to track a fixed number of distinct quantities with concurrently
 * provided samples and to permit merging statistics from multiple instances.
 * Storage overhead reduced relative to Cook's presentation.
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
    running_statistics& operator()(const Real x[N]);

    /** Incorporate running information from another instance. */
    running_statistics& operator()(const running_statistics& o);

    /** Obtain the running number of samples provided thus far. */
    std::size_t count() const
    { return n_; }

    /** Obtain the minimum sample observed for quantity \c i. */
    Real min(std::size_t i) const
    { assert(i < N); return min_[i]; }

    /** Obtain the maximum sample observed for quantity \c i. */
    Real max(std::size_t i) const
    { assert(i < N); return max_[i]; }

    /** Obtain the running mean for quantity \c i. */
    Real avg(std::size_t i) const
    { assert(i < N); return M_[i]; }

    /** Obtain the running sample variance for quantity \c i. */
    Real var(std::size_t i) const
    { assert(i < N); return n_ == 1 ? 0 : S_[i] / (n_ - 1); }

    /** Obtain the running sample standard deviation for quantity \c i. */
    Real std(std::size_t i) const
    { using std::sqrt; return sqrt(var(i)); }

    /** Reset the instance to its newly constructed state. */
    void clear();

    /** Number of quantities per <tt>operator()(const T*)</tt> invocation */
    static const std::size_t static_size = N;        // Compile time

    /** Number of quantities per <tt>operator()(const T*)</tt> invocation */
    std::size_t size() const { return static_size; } // Run time

private:

    Real M_[N], S_[N], min_[N], max_[N];

    std::size_t n_;

};

template <typename Real, std::size_t N>
running_statistics<Real,N>::running_statistics()
{
    clear();
}

template <typename Real, std::size_t N>
running_statistics<Real,N>& running_statistics<Real,N>::operator()(
    const Real x[N])
{
    // Algorithm from Knuth TAOCP vol 2, 3rd edition, page 232.
    // Knuth shows better behavior than Welford 1962 on test data.

    using std::max;
    using std::min;

    ++n_;

    if (n_ > 1) {  // Second and subsequent invocation

        const Real inv_n = Real(1) / n_;
        for (std::size_t i = 0; i < N; ++i) {
            const Real d  = x[i] - M_[i];
            M_[i]        += d * inv_n;
            S_[i]        += d * (x[i] - M_[i]);
            min_[i]       = min(min_[i], x[i]);
            max_[i]       = max(max_[i], x[i]);
        }

    } else {       // First invocation requires special treatment

        for (std::size_t i = 0; i < N; ++i) {
            M_[i]   = x[i];
            S_[i]   = 0;
            min_[i] = x[i];
            max_[i] = x[i];
        }
    }

    return *this;
}

template <typename Real, std::size_t N>
running_statistics<Real,N>& running_statistics<Real,N>::operator()(
    const running_statistics& o)
{
    using std::max;
    using std::min;

    const std::size_t total = n_ + o.n_;  // How many samples in combined data?
    if (o.n_ == 0) return *this;          // Other contains no data; run away
    if (n_   == 0) return *this = o;      // We contain no data; simply copy
    assert(total >= 1);                   // Voila, degeneracy sidestepped

    // Combine statistics from both instances into this
    for (std::size_t i = 0; i < N; ++i) {
        const Real dM = M_[i] - o.M_[i];  // Cancellation issues?
        M_[i]         = (n_ * avg(i) + o.n_ * o.avg(i)) / total;
        S_[i]         = (n_   == 1 ? 0 :   S_[i])
                      + (o.n_ == 1 ? 0 : o.S_[i])
                      + ((dM * dM) * (n_ * o.n_)) / total;
        min_[i]       = min(min_[i], o.min_[i]);
        max_[i]       = max(max_[i], o.max_[i]);
    }
    n_ = total;

    return *this;
}

template <typename Real, std::size_t N>
void running_statistics<Real,N>::clear()
{
    using std::numeric_limits;
    for (std::size_t i = 0; i < N; ++i)
        M_[i] = S_[i] = min_[i] = max_[i] = numeric_limits<Real>::quiet_NaN();
    n_ = 0;
}

} // namespace suzerain

#endif // SUZERAIN_RUNNING_STATISTICS_HPP
