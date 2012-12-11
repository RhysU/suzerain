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
// noise_specification.hpp: classes handling specifying additive noise
// $Id$

#ifndef SUZERAIN_NOISE_SPECIFICATION_HPP
#define SUZERAIN_NOISE_SPECIFICATION_HPP

#include <suzerain/common.hpp>

// TODO Distinguish between two- versus one-sided stretching in grid_specification

/** @file
 * Provides classes handling specifying additive noise
 */

namespace suzerain {

/**
 * Options for adding random noise to momentum fields.
 */
class noise_specification
{

public:

    /** Construct an instance with the given default values */
    explicit noise_specification(real_t fluct_percent = 0,
                                 unsigned long fluct_seed = 12345);

    /**
     * Maximum fluctuation magnitude to add as a percentage
     * of centerline streamwise momentum.
     */
    real_t percent;

    /**
     * Fraction of the X direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kxfrac_min;

    /**
     * Fraction of the X direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kxfrac_max;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kzfrac_min;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kzfrac_max;

    /** \ref rngstream generator seed (see L'Ecuyer et al. 2002) */
    unsigned long seed;

};

} // namespace suzerain

#endif // SUZERAIN_NOISE_SPECIFICATION_HPP
