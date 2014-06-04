//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
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
 * @copydoc shared_esio_handle.hpp
 */

#include <suzerain/support/shared_esio_handle.hpp>

#include <esio/esio.h>
#include <esio/error.h>

namespace suzerain {

namespace support {

shared_esio_handle::shared_esio_handle(MPI_Comm comm)
    : shared_ptr<boost::remove_pointer<esio_handle>::type>(
                esio_handle_initialize(comm),
                esio_handle_finalize)
{
}

} // end namespace support

} // end namespace suzerain
