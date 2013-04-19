//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_NOISE_DEFINITION_HPP
#define SUZERAIN_SUPPORT_NOISE_DEFINITION_HPP

/** @file
 * Routines for adding noise/perturbations to state fields
 */

#include <suzerain/common.hpp>
#include <suzerain/noise_specification.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/**
 * Upgrades a \ref noise_specification with \ref definition_base behavior.
 * This permits using the instance with \ref program_options.
 */
class noise_definition
    : public noise_specification
    , public virtual definition_base
{
public:

    /** Construct an instance with the given default values */
    explicit noise_definition(
            real_t fluct_percent     = 0,
            unsigned long fluct_seed = 12345,
            real_t kxfrac_min        = 0,
            real_t kxfrac_max        = 1,
            real_t kzfrac_min        = 0,
            real_t kzfrac_max        = 1);

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_NOISE_DEFINITION_HPP
