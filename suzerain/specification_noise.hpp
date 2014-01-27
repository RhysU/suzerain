//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_NOISE_SPECIFICATION_HPP
#define SUZERAIN_NOISE_SPECIFICATION_HPP

/** @file
 * Provides classes handling specifying additive noise
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Options for adding random noise to momentum fields.
 */
class noise_specification
{

public:

    /** Construct an instance with the given default values */
    explicit noise_specification(
            real_t        fluct_percent = 0,
            unsigned long fluct_seed    = 12345,
            real_t        kxfrac_min    = 0,
            real_t        kxfrac_max    = 1,
            real_t        kzfrac_min    = 0,
            real_t        kzfrac_max    = 1);

    /**
     * Maximum fluctuation magnitude to add as a percentage
     * of centerline streamwise momentum.
     */
    real_t percent;

    /** \ref rngstream generator seed (see L'Ecuyer et al. 2002) */
    unsigned long seed;

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

};

} // namespace suzerain

#endif // SUZERAIN_NOISE_SPECIFICATION_HPP
