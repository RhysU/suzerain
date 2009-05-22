/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * exceptions.hpp: Exceptions used within Suzerain
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_EXCEPTIONS_H
#define PECOS_SUZERAIN_EXCEPTIONS_H

#include <boost/exception.hpp>

namespace pecos
{

namespace suzerain
{

/** Reports arguments to functions that are outside the valid input range.
 *
 * \internal Intended to have the same semantics as \c std::domain_error
 * but with the benefit of boost::exception as a base class.
 */
class domain_error: public boost::exception { };

} // namespace suzerain

} // namespace pecos

#endif // PECOS_SUZERAIN_EXCEPTIONS_H
