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
#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

static void parse_nonnegative(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_nonnegative(v, n);
    *t = v;
}

namespace support {

static void parse_formulation(const std::string& s, largo_formulation *t)
{
    *t = largo_formulation::lookup(s);
}

std::map<std::string,const largo_formulation*> largo_formulation::by_name;

largo_formulation::largo_formulation(
        const int v,
        const char *n,
        const char *d)
    : v(v), n(n), d(d)
{
    // Register for lookup of instances by whitespace trimmed name
    by_name[boost::algorithm::trim_copy(std::string(this->n))] = this;
}

const largo_formulation&
largo_formulation::lookup(const std::string& name)
{
    using namespace std;
    using namespace boost::algorithm;
    map<string,const largo_formulation*>::const_iterator i = by_name.find(trim_copy(name));
    if (i == by_name.end()) {
        ostringstream oss;
        oss << "Unknown largo_formulation '" << name << "'";
        throw invalid_argument(oss.str());
    } else {
        return *((*i).second);
    }
}

// This is ugly and wasteful.
std::set<std::string>
largo_formulation::names()
{
    using namespace std;
    set<string> retval;
    map<string,const largo_formulation*>::const_iterator i   = by_name.begin();
    map<string,const largo_formulation*>::const_iterator end = by_name.end();
    while (i != end) {
        retval.insert((*i).first);
        ++i;
    }
    return retval;
}

// BEGIN Add known formulations here
const largo_formulation largo_formulation::disable(
        0, "disable", "No slow growth formulation is in use");

const largo_formulation largo_formulation::temporal(
        1, "temporal", "Original temporal formulation by Topalian et al.");

const largo_formulation largo_formulation::spatial(
        2, "spatial", "Full spatial formulation by Topalian et al.");
// END Add known formulations here

largo_definition::largo_definition()
    : formulation(largo_formulation::disable)
    , grdelta    (std::numeric_limits<real_t>::quiet_NaN())
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_grdelta[]            = "grdelta";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_grdelta[]            = "Growth rate of reference thickness (Delta)";


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

    options_description retval("Largo-based slow growth parameters");

    // largo_formulation
    // Build help message so it contains list of known formulations
    std::set<string> names = largo_formulation::names();
    std::ostringstream largo_formulation_help;
    largo_formulation_help
        << "Which, if any, slow growth forcing should be added during time advance?"
        << " { ";
    copy(names.begin(), names.end(),
         std::ostream_iterator<string>(largo_formulation_help, " "));
    largo_formulation_help
        << "}";
    retval.add_options()
    ("largo_formulation",
     value<string>()
        ->default_value(largo_formulation::disable.name())
        ->notifier(bind(&parse_formulation, _1, &formulation)),
     largo_formulation_help.str().c_str())
    ;

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.
    auto_ptr<typed_value<string> > p;

    // grdelta
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &grdelta, "largo_grdelta"));
    if (!(boost::math::isnan)(grdelta)) {
        p->default_value(lexical_cast<string>(grdelta));
    }
    retval.add_options()("largo_grdelta", p.release(), desc_grdelta);

    return retval;
}

static const char location[] = "largo";

void
largo_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing largo_definition parameters");

    if (!formulation.enabled()) { return; }  // Shortcircuit on no formulation

    // Write out the "container" holding all other settings
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
            "Is a largo-based slow growth formulation in use?");

    // Write out the formulation name
    esio_string_set(h, location, "formulation", formulation.name().c_str());

    // scalars
    if (formulation == largo_formulation::disable) {
        // Nothing else to save
    } else if (formulation == largo_formulation::temporal) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::spatial) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else {
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

}

void
largo_definition::load(
        const esio_handle h,
        const bool verbose)
{
    // Overwrite instance members with a "disable" formulation with all NaNs
    *this = largo_definition();
    assert(formulation == largo_formulation::disable);

    // Only proceed if a largo definition is active in the restart
    int in_use = 0;
    esio_line_establish(h, 1, 0, 1); // All ranks load any data
    if (ESIO_NOTFOUND != esio_line_size(h, location, NULL)) {
        esio_line_read(h, location, &in_use, 0);
    }
    if (!in_use) return;

    DEBUG0("Loading largo_definition parameters");

    // Load formulation name and look it up in the static instance map
    {
        char *name = esio_string_get(h, location, "formulation");
        this->formulation = largo_formulation::lookup(name);
        free(name);
    }

    if (formulation == largo_formulation::disable) {
        // Nothing else to load
    } else if (formulation == largo_formulation::temporal) {
        esio_attribute_read(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::spatial) {
        esio_attribute_read(h, location, name_grdelta, &grdelta);
    } else {
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }
}

} // namespace support

} // namespace suzerain
