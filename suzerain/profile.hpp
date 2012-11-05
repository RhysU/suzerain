//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// profile.hpp: miscellaneous utilities for profiling parallel transposes
// $Id$

#ifndef SUZERAIN_PROFILE_HPP
#define SUZERAIN_PROFILE_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides miscellaneous utilities for profiling parallel transposes.
 */

namespace suzerain {

/** Provides miscellaneous parallel transpose profiling flags.  */
class ProfileDefinition : public problem::definition_base
{
public:

    /** Default constructor */
    ProfileDefinition();

    /*! @copydoc problem::definition_base::options */
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

#endif // SUZERAIN_PROFILE_HPP
