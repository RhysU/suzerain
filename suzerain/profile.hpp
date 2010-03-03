/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * profile.hpp: miscellaneous utilities for profiling parallel transposes
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_PROFILE_HPP
#define __SUZERAIN_PROFILE_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides miscellaneous utilities for profiling parallel transposes.
 */

namespace suzerain {

/** Provides miscellaneous parallel transpose profiling flags.  */
class ProfileDefinition : public ::suzerain::problem::IDefinition
{
public:

    /** Default constructor */
    ProfileDefinition();

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

    /**
     * Retrieve howmany scalar values make up a single state field.
     * Real-valued fields equate to <tt>howmany = 1</tt> while complex fields
     * require <tt>howmany = 2</tt>.  This multiplier should be used when
     * determining the amount of scalar state is transformed using backward(),
     * forward(), or roundtrip().
     *
     * @return How many scalar values comprise a single state field.
     */
    int howmany() const { return howmany_; };

    /**
     * Retrieve the number of state fields to transpose backward.
     * This typically means from wave space to physical space.
     *
     */
    int backward() const { return backward_; };

    /**
     * Retrieve the number of state fields to transpose forward.
     * This typically means from physical space to wave space.
     */
    int forward() const { return forward_; };

    /**
     * Retrieve the number of state fields to transpose both backward and
     * forward.
     */
    int roundtrip() const { return roundtrip_; };

    /**
     * Retrieve the number of repetitions to execute in the profiling test.
     */
    int nrep() const { return nrep_; };

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

    int howmany_;   /**< Stores howmany parameter value */
    int backward_;  /**< Stores backward parameter value */
    int forward_;   /**< Stores forward parameter value */
    int roundtrip_; /**< Stores roundtrip parameter value */
    int nrep_;      /**< Stores nrep parameter value */
};

} // namespace suzerain

#endif // __SUZERAIN_PROFILE_HPP
