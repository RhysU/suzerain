//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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
 * @copydoc samples.hpp
 */

#include <suzerain/samples.hpp>

namespace suzerain {

samples::samples()
    : t(std::numeric_limits<real_t>::quiet_NaN())
{
}

samples::samples(real_t)
    : t(t)
{
}

samples::samples(real_t t,
                 samples::storage_type::Index Ny)
    : t(t)
    , storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
{
}

} // namespace suzerain
