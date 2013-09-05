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

// #ifdef SUZERAIN_HAVE_LARGO
// #include <largo/largo.h>
// #endif

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

class largo_formulation
{
private:

    int         v;  ///< A quickly comparable value
    std::string n;  ///< A brief, human-readable name
    std::string d;  ///< A relatively complete description

    /** Maintains map from name to pointer-to-static instances. */
    static std::map<std::string,const largo_formulation*> by_name;

    /** Create a new (static) instance and register it in \c by_name. */
    largo_formulation(const int v, const char *n, const char *d);

public:

    /**
     * The known slow growth formulation types.
     * @{
     */
    static const largo_formulation disable;
    static const largo_formulation temporal;
    static const largo_formulation spatial;
    static const largo_formulation temporal_tensor_consistent;
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

    /**
     * Return an instance with the given formulation name.
     * Leading and trailing white spaces are ignored.
     * @throws std::invalid_argument on unknown name.
     */
    static const largo_formulation& lookup(const std::string& name);

    /** Return the set of known formulation names. */
    static std::set<std::string> names();

};

template< typename CharT, typename Traits >
std::basic_ostream<CharT,Traits>& operator<<(
        std::basic_ostream<CharT,Traits> &os,
        const largo_formulation &f)
{
    return os << f.name();
}

/**
 * Holds parameters defining largo boundary cases.
 */
class largo_definition
    : public virtual definition_base
    , public virtual loadable
    , public virtual overridable<largo_definition>
    , public virtual populatable<largo_definition>
    , public virtual savable
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    largo_definition();

    /** @copydoc populatable::populate */
    virtual void populate(
            const largo_definition& that,
            const bool verbose = false);

    /** @copydoc overridable::override */
    virtual void override(
            const largo_definition& that,
            const bool verbose = false);

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /** Which \ref largo_formulation is in use? */
    largo_formulation formulation;

    /** Growth rate of reference thickness \f$\Delta\f$ */
    real_t grdelta;

    /** Pointer to largo workspace */
    void * workspace;

    /** Get baseflow field and derivatives from coefficients*/
    void get_baseflow(
         const real_t    y,
         real_t *     base,
         real_t *   dybase,
         real_t *   dxbase);

    /** Get baseflow pressure and derivatives from coefficients*/
    void get_baseflow_pressure(
         const real_t    y,
         real_t &    Pbase,
         real_t &  dyPbase,
         real_t &  dxPbase);

private:
    /** Baseflow coefficients for the field*/
    MatrixXXr x;

    /** Baseflow coefficient base */
    std::string x_base;

    /** Baseflow coefficients for the field derivative*/
    MatrixXXr dx;

    /** Baseflow derivative coefficient base */
    std::string dx_base;

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_LARGO_DEFINITION_HPP
