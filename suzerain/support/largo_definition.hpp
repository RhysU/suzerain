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

#ifndef SUZERAIN_SUPPORT_LARGO_DEFINITION_HPP
#define SUZERAIN_SUPPORT_LARGO_DEFINITION_HPP

/** @file
 * Provides \ref largo_definition.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

class largo_formulation
    : boost::noncopyable
{
private:

    largo_formulation(const int v, const char *n, const char *d)
        : v(v), n(n), d(d)
    {
        // NOP
    }

    const int         v;  ///< A quickly comparable value
    const std::string n;  ///< A brief, human-readable name
    const std::string d;  ///< A relatively complete description

public:

    /**
     * The known slow growth formulation types.
     * @{
     */

    static const largo_formulation disable;
    static const largo_formulation temporal;
    static const largo_formulation spatial;

    /**@}*/

    /** Is a slow growth formulation in use? */
    bool enabled() const
    { return v != 0; }

    /** What is the name of the given formulation? */
    const std::string& name() const
    { return n; }

    /** What is the description of the given formulation? */
    const std::string& description() const
    { return d; }

    /** Is \c this the same formulation as \c that? */
    bool operator==(const largo_formulation& that) const
    { return this->v == that.v; }

    /** Is \c this a different formulation from \c that? */
    bool operator!=(const largo_formulation& that) const
    { return this->v != that.v; }

};

/**
 * Holds parameters defining largo boundary cases.
 */
class largo_definition
    : public virtual definition_base
    , public virtual loadable
    , public virtual savable
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    largo_definition();

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_LARGO_DEFINITION_HPP
