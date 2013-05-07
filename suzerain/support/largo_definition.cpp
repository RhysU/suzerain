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
 * @copydoc largo_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/largo_definition.hpp>

#include <esio/error.h>

// FIXME: Include only needed headers
//        Review after basic functionality is implemented

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

largo_definition::largo_definition()
{
}


boost::program_options::options_description
largo_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("Largo parameters");


    auto_ptr<typed_value<string> > p;

    return retval;
}

void
largo_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing largo_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
}

void
largo_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading largo_definition parameters");

    largo_definition t;

    // All ranks load

}

} // namespace support

} // namespace suzerain
