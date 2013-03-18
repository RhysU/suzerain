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

#include <cmath>
#include <limits>

namespace suzerain {

/**
 * Accumulates running minimum, maximum, mean, and variance details from a data
 * stream.  Adapted from http://www.johndcook.com/standard_deviation.html.
 * Extended to track a fixed number of quantities sampled together.
 *
 * @tparam N       Number of distinct statistics to track
 * @tparam Real    Floating point type used for input and accumulation
 * @tparam Integer Integral type used for accumulation
 */
template <std::size_t N,
          typename Real = double>
class running_statistics
{
public:

    running_statistics()
    {
        reset();
    }

    Integer count() const
    {
        return n_;
    }

    Real minimum(std::size_t i) const
    {
        return min_[i];
    };

    Real maximum(std::size_t i) const
    {
        return max_[i];
    };

    Real mean(std::size_t i) const
    {
        using std::numeric_limits;
        return (n_ > 0) ? newM_[i] : numeric_limits<Real>::quiet_NaN();
    }

    Real variance() const
    {
        using std::numeric_limits;
        return  (n_ > 1) ? newS_/(n_ - 1)
              : (n_ = 1) ? 0
              :            numeric_limits<Real>::quiet_NaN();
    }

    Real stddev() const
    {
        using std::sqrt;
        return sqrt(variance());
    }

    void reset()
    {
        n_ = 0;
        using std::numeric_limits;
        for (std::size_t i = 0; i < N; ++i)
            min_[i] = max_[i] = numeric_limits<Real>::quiet_NaN();
    }

    void operator()(const Real* x)
    {
        using std::min;
        using std::max;
        using std::size_t;

        // See Knuth TAOCP vol 2, 3rd edition, page 232
        ++n_;

        if (n_ > 1) {  // Second and subsequent iteration

            // Process moments for current iteration
            const Real inv_n = static_cast<Real>(1) / n_;
            for (size_t i = 0; i < N; ++i)
                newM_[i] = oldM_[i] + (x[i] - oldM_[i]) * inv_n;
            for (size_t i = 0; i < N; ++i)
                newS_[i] = oldS_[i] + (x[i] - oldM_[i])*(x[i] - newM_[i]);

            // Prepare to process moments for next iteration
            for (size_t i = 0; i < N; ++i)
                oldM_[i] = newM_[i];
            for (size_t i = 0; i < N; ++i)
                oldS_[i] = newS_[i];

            // Process extrema for current iteration
            for (size_t i = 0; i < N; ++i)
                min_[i] = min(min_[i], x[i]);
            for (size_t i = 0; i < N; ++i)
                max_[i] = max(max_[i], x[i]);

        } else {       // Initial iteration

            // Initialize moments
            for (size_t i = 0; i < N; ++i)
                oldM_[i] = newM_[i] = x[i];
            for (size_t i = 0; i < N; ++i)
                oldS_[i] = 0;

            // Initialize extrama
            for (size_t i = 0; i < N; ++i)
                min_[i] = x[i];
            for (size_t i = 0; i < N; ++i)
                max_[i] = x[i];

        }
    }

private:

    Real oldM_[N], newM_[N], oldS_[N], newS_[N], min_[N], max_[N];

    std::size_t n_;
};

} // namespace suzerain

#endif // SUZERAIN_RUNNING_STATISTICS_HPP
