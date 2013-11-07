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

#ifndef SUZERAIN_LARGO_FORMULATION_HPP
#define SUZERAIN_LARGO_FORMULATION_HPP

/** @file
 * Provides \ref largo_formulation.
 */

#include <suzerain/common.hpp>

namespace suzerain {

class largo_formulation
{

public:

    /**
     * The known slow growth formulation types.
     * @{
     */
    static const largo_formulation disable;
    static const largo_formulation temporal;
    static const largo_formulation spatial;
    static const largo_formulation temporal_tensor_consistent;
    static const largo_formulation spatiotemporal;
    static const largo_formulation temporal_consistent;
    static const largo_formulation spatiotemporal_consistent;
    /**@}*/

    /** Is a slow growth formulation in use? */
    bool enabled() const
    { return v != 0; }

    /** What is the name of the given formulation? */
    const std::string& name() const
    { return n; }

    /**
     * Is the given formulation strictly-temporal?
     * That is, is spatial homogenization totally absent?
     */
    bool is_strictly_temporal() const
    { return t; }

    /** What is the description of the given formulation? */
    const std::string& description() const
    { return d; }

    /**
     * Does the model expect conserved state growth rates?  That is, provide
     * \f$\rho\f$, \f$\rho u\f$, \f$\rho v\f$, \f$\rho w\f$, \f$\rho E\f$ and
     * \f$p\f$ within <code>grDA</code> and <code>grDArms</code> when calling
     * <code>largo_init</code>.
     */
    bool expects_conserved_growth_rates() const
    { return g == grspec_conserved; }

    /**
     * Does the model expect specific state growth rates?  That is, provide
     * \f$\rho\f$, \f$u\f$, \f$v\f$, \f$w\f$, \f$E\f$ and \f$p\f$ within
     * <code>grDA</code> and <code>grDArms</code> when calling
     * <code>largo_init</code>.
     */
    bool expects_specific_growth_rates() const
    { return g == grspec_specific; }

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

private:

    /**
     * Precisely what growth rate information does the model
     * require when calling <code>largo_init</code>?
     */
    enum grspec_type {

        /** This detail is nonsensical or unknown. */
        grspec_unknown = 0,

        /** The model expects conserved state growth rates. */
        grspec_conserved,

        /** The model expects specific state growth rates. */
        grspec_specific

    };

    int         v;  ///< A quickly-comparable, unique value
    std::string n;  ///< A brief, human-readable name
    bool        t;  ///< Is the formulation strictly-temporal?
    std::string d;  ///< A relatively complete description
    grspec_type g;  ///< How are growth rates specified?

    /** Maintains map from name to pointer-to-static instances. */
    static std::map<std::string,const largo_formulation*> by_name;

    /** Register names within \c by_name. */
    static void register_name(const std::string& name,
                              const largo_formulation* instance);

    /** Create a new (static) instance and register it. */
    largo_formulation(const int         v,
                      const char*       n,
                      const bool        t,
                      const char*       d,
                      const grspec_type g);

    /**
     * Create a new (static) instance and register it.
     *
     * Register \c misspellings in \c by_name as well to permit lookup
     * by any of them in addition to \c n.
     */
    largo_formulation(const int         v,
                      const char*       n,
                      const bool        t,
                      const char*       d,
                      const grspec_type g,
                      const std::vector<std::string>& misspellings);

};

template< typename CharT, typename Traits >
std::basic_ostream<CharT,Traits>& operator<<(
        std::basic_ostream<CharT,Traits> &os,
        const largo_formulation &f)
{
    return os << f.name();
}

} // namespace suzerain

#endif // SUZERAIN_LARGO_FORMULATION_HPP
