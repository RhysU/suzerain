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

/** @file
 * @copydoc noise_specification.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/noise_specification.hpp>

namespace suzerain {

noise_specification::noise_specification(
        real_t        fluct_percent,
        unsigned long fluct_seed,
        real_t        kxfrac_min,
        real_t        kxfrac_max,
        real_t        kzfrac_min,
        real_t        kzfrac_max)
    : percent   (fluct_percent)
    , seed      (fluct_seed)
    , kxfrac_min(kxfrac_min)
    , kxfrac_max(kxfrac_max)
    , kzfrac_min(kzfrac_min)
    , kzfrac_max(kzfrac_max)
{
}

} // namespace suzerain
