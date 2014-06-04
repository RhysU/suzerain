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

#ifndef SUZERAIN_SUPPORT_SHARED_ESIO_HANDLE_HPP
#define SUZERAIN_SUPPORT_SHARED_ESIO_HANDLE_HPP

/** @file
 * Provides \ref shared_esio_handle.
 */

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/support/esio_fwd.hpp>

namespace suzerain {

namespace support {

/**
 * Smart pointer semantics for managing an \ref esio_handle.  Provided because
 * ESIO's declaration structure makes creating \c shared_ptr instances and to
 * ease setting up automatic allocation/deallocation through the C API.
 */
class shared_esio_handle
    : public shared_ptr<boost::remove_pointer<esio_handle>::type>
{
public:

    /** Constructor producing an empty handle evaluating to false. */
    shared_esio_handle();

    /** Constructor producing a valid handle across \c comm. */
    explicit shared_esio_handle(MPI_Comm comm);

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_SHARED_ESIO_HANDLE_HPP
