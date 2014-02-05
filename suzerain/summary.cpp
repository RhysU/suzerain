//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

/** @file
 * @copydoc summary.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/summary.hpp>

namespace suzerain {

summary::summary()
{
}

summary::summary(summary::storage_type::Index Ny)
    : storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
{
}

// An implementation consistency check.  Suspect the macros when it fails.
BOOST_STATIC_ASSERT(summary::nscalars::total == SUZERAIN_SUMMARY_COUNT);

// An implementation consistency check.  Again, the macros probably broke.
#define CHECK(data, name, description, for_each_offset) \
        BOOST_STATIC_ASSERT( for_each_offset == summary::offset::name );
SUZERAIN_SUMMARY_FOR_EACH(CHECK,)
#undef CHECK

#define DESCRIPTION(s, data, name_description_tuple) \
        BOOST_PP_TUPLE_ELEM(2, 1, name_description_tuple)
const char * summary::description[summary::nscalars::total] = {
SUZERAIN_SUMMARY_ENUM_TRANSFORM(DESCRIPTION,)
};
#undef DESCRIPTION


} // namespace suzerain
