//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc definition_base.hpp
 */

#include <suzerain/support/overridable.hpp>

namespace suzerain {

namespace support {

// For maybe_XXX_impl to indicates a \c real_t is a default value
static bool default_value_real_t(const real_t& v)
{ return (boost::math::isnan)(v); }

// For maybe_XXX_impl to indicates an \c int is a default value
static bool default_value_int(const int& v)
{ return v == 0; }

bool
maybe_override(const char*   name,
               const char*   description,
                     real_t& destination,
               const real_t& source,
               const bool    verbose)
{
    return internal::maybe_override_impl(
            name, description, destination, source,
            verbose, &default_value_real_t);
}

bool
maybe_override(const char* name,
               const char* description,
                     int&  destination,
               const int&  source,
               const bool  verbose)
{
    return internal::maybe_override_impl(
            name, description, destination, source,
            verbose, &default_value_int);
}

} // namespace support

} // namespace suzerain
