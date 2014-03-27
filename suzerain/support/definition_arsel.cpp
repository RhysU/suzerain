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
 * @copydoc definition_noise.hpp
 */

#include <suzerain/support/definition_arsel.hpp>

#include <boost/version.hpp>
#include <esio/esio.h>

#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

/** Helper used to parse size_t-based options */
static void parse_size_t(const std::string& s,
                         std::size_t* value,
                         const char* name)
{
    using std::floor;
#pragma warning(push,disable:2259)
    const real_t t = floor(exprparse<real_t>(s, name) + real_t(1)/2);
#pragma warning(pop)
    validation::ensure_nonnegative(t, name);
    *value = t;
}

/** Helper used to parse string-based options */
template<typename T>
static void parse_option(const std::string& s,
                         T* value,
                         void (*validator)(T, const char*),
                         const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

definition_arsel::definition_arsel()
    : specification_arsel()
{
}

// Strings used in options_description and populate/override/save/load.
static const char location      [] = "arsel";
static const char desc_location []
    = "Autoregressive autocorrelation analysis settings governing"
      " eff_N, eff_var, mu, mu_sigma, p, and T attributes on"
      " Reynolds averaged quantities. See 'Estimating Uncertainties"
      " in Statistics Computed from DNS' by Oliver et al. in PoF.";

static const char name_absrho   [] = "arsel_absrho";
static const char name_criterion[] = "arsel_criterion";
static const char name_maxorder [] = "arsel_maxorder";
static const char name_minorder [] = "arsel_minorder";
static const char name_wlenT0   [] = "arsel_wlenT0";

static const char * const attr_absrho    = name_absrho    + sizeof(location);
static const char * const attr_criterion = name_criterion + sizeof(location);
static const char * const attr_maxorder  = name_maxorder  + sizeof(location);
static const char * const attr_minorder  = name_minorder  + sizeof(location);
static const char * const attr_wlenT0    = name_wlenT0    + sizeof(location);

static const char desc_absrho   []
    = "Integrate absolute autocorrelation when determining T0";
static const char desc_criterion[]
    = "Employ the specified model selection criterion";
static const char desc_maxorder []
    = "Consider only models of at most order AR(p=MAX)";
static const char desc_minorder []
    = "Consider only models of at least order AR(p=MIN)";
static const char desc_wlenT0   []
    = "Integrate for T0 until WLEN times the input length";

boost::program_options::options_description
definition_arsel::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_positive;

    boost::program_options::options_description retval(desc_location);

    retval.add_options()
        (name_minorder, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("MIN")
#endif
         ->notifier(bind(&parse_size_t, _1, &minorder, name_minorder))
         ->default_value(lexical_cast<string>(minorder)),
         desc_minorder)
        (name_maxorder, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("MAX")
#endif
         ->notifier(bind(&parse_size_t, _1, &maxorder, name_maxorder))
         ->default_value(lexical_cast<string>(maxorder)),
         desc_maxorder)
        (name_criterion,
         value<std::string>()
         ->notifier(bind(&definition_arsel::criterion, this, _1))
         ->default_value(criterion()),
         desc_criterion)
        (name_absrho,
         bool_switch(&absrho)->default_value(absrho),
         desc_absrho)
        (name_wlenT0, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("WLEN")
#endif
         ->notifier(bind(&parse_option<real_t>, _1, &wlenT0,
                         &ensure_positive<real_t>, name_wlenT0))
         ->default_value(lexical_cast<string>(wlenT0)),
         desc_wlenT0)
        ;

    return retval;
}

void
definition_arsel::save(
        const esio_handle h) const
{
    DEBUG0("Storing definition_arsel parameters");

    // Only root writes the containing location
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0, desc_location);

    // Everyone writes the metadata
    const int absrho   = this->absrho;
    const int minorder = this->minorder;
    const int maxorder = this->maxorder;
    esio_attribute_write(h, location, attr_absrho,    &absrho);
    esio_attribute_write(h, location, attr_minorder,  &minorder);
    esio_attribute_write(h, location, attr_maxorder,  &maxorder);
    esio_attribute_write(h, location, attr_wlenT0,    &wlenT0);
    esio_string_set     (h, location, attr_criterion, criterion().c_str());

}

} // end namespace support

} // end namespace suzerain
