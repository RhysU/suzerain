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

#include <boost/assign/list_of.hpp>

#include <esio/error.h>
#include <esio/esio.h>

#include <gsl/gsl_poly.h>

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

void
largo_formulation::register_name(const std::string& name,
                                 const largo_formulation* instance)
{
    const std::string& trimmed = boost::algorithm::trim_copy(name);
    if (by_name.find(trimmed) != by_name.end()) {
        throw std::logic_error(std::string("Name collision on '")
            + trimmed + "' when registering with largo_formulation::by_name");
    }
    by_name[trimmed] = instance;
}

largo_formulation::largo_formulation(
        const int v,
        const char *n,
        const char *d)
    : v(v), n(n), d(d)
{
    register_name(this->n, this);
}

largo_formulation::largo_formulation(
        const int v,
        const char *n,
        const char *d,
        const std::vector<std::string>& misspellings)
    : v(v), n(n), d(d)
{
    register_name(this->n, this);

    for (std::size_t i = 0; i < misspellings.size(); ++i) {
        register_name(misspellings[i], this);
    }
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
        1, "bl_temporal", "Original temporal formulation by Topalian et al.",
        boost::assign::list_of("temporal")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::spatial(
        2, "bl_spatial", "Full spatial formulation by Topalian et al.",
        boost::assign::list_of("spatial")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::temporal_tensor_consistent(
        3, "bl_temporal_tensor-consistent",
           "Temporal tensor-consistent formulation by Topalian et al.",
        boost::assign::list_of("bl_temporal_tensor_consistent")
                              ("temporal_tensor_consistent")
                              ("temporal_tensor-consistent")
            .convert_to_container<std::vector<std::string> >());
// END Add known formulations here

largo_definition::largo_definition()
    : formulation(largo_formulation::disable)
    , grdelta    (std::numeric_limits<real_t>::quiet_NaN())
{
    (this->workspace) = NULL;
    (this->x).resize(0,0);
    (this->dx).resize(0,0);
}

// Strings used in options_description and populate/override/save/load.
static const char name_formulation[] = "formulation";
static const char name_grdelta[]            = "grdelta";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_formulation[] = "Name of the slow growth formulation";
static const char desc_grdelta[]     = "Growth rate of reference thickness (Delta)";


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

static const char location[]            = "largo";
static const char location_baseflow[]   = "largo_baseflow";
static const char location_baseflow_d[] = "largo_baseflow_derivative";

// For maybe_XXX_impl to indicates a \c largo_definition is a default value
static bool default_value_formulation(const largo_formulation& v)
{ return !v.enabled(); }

static bool maybe_populate(const char*               name,
                           const char*               description,
                                 largo_formulation&  destination,
                           const largo_formulation&  source,
                           const bool verbose)
{
    return internal::maybe_populate_impl(
            name, description, destination, source,
            verbose, &default_value_formulation);
}

void
largo_definition::populate(
        const largo_definition& that,
        const bool verbose)
{
    using support::maybe_populate;
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(formulation);
    CALL_MAYBE_POPULATE(grdelta);
#undef CALL_MAYBE_POPULATE
}

static bool maybe_override(const char*               name,
                           const char*               description,
                                 largo_formulation&  destination,
                           const largo_formulation&  source,
                           const bool verbose)
{
    return internal::maybe_override_impl(
            name, description, destination, source,
            verbose, &default_value_formulation);
}

void
largo_definition::override(
        const largo_definition& that,
        const bool verbose)
{
    using support::maybe_override;
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(formulation);
    CALL_MAYBE_OVERRIDE(grdelta);
#undef CALL_MAYBE_OVERRIDE
}

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
    } else if (formulation == largo_formulation::temporal_tensor_consistent) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else {
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

    // Write baseflow field coefficients (only master process)
    if (this->x.outerSize() > 0) {
        esio_plane_establish(h,
            this->x.outerSize(), 0, procid == 0 ?  this->x.outerSize() : 0,
            this->x.innerSize(), 0, procid == 0 ?  this->x.innerSize() : 0);
        esio_plane_write(h, location_baseflow, this->x.data(),
            this->x.outerStride(), this->x.innerStride(),
            "Baseflow field coefficients");

        // Write out the baseflow coefficient base
        esio_string_set(h, location_baseflow, "coefficient_base",
                        this->x_base.c_str());
    }

    // Write baseflow derivative coefficients (only master process)
    if (this->dx.outerSize() > 0) {
        esio_plane_establish(h,
            this->dx.outerSize(), 0, procid == 0 ?  this->dx.outerSize() : 0,
            this->dx.innerSize(), 0, procid == 0 ?  this->dx.innerSize() : 0);
        esio_plane_write(h, location_baseflow_d, this->dx.data(),
            this->dx.outerStride(), this->dx.innerStride(),
            "Baseflow derivative coefficients");

        // Write out the baseflow coefficient base
        esio_string_set(h, location_baseflow_d, "coefficient_base",
                        this->dx_base.c_str());
    }
}

void
largo_definition::load(
        const esio_handle h,
        const bool verbose)
{

    largo_definition t;
    assert(t.formulation == largo_formulation::disable);

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
        t.formulation = largo_formulation::lookup(name);
        free(name);
    }

    if (t.formulation == largo_formulation::disable) {
        // Nothing else to load
    } else if (t.formulation == largo_formulation::temporal) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::spatial) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::temporal_tensor_consistent) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else {
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

    this->populate(t, verbose);  // Prefer this to incoming

    // Load baseflow coefficients
    // Only proceed if largo_baseflow is present in the restart
    int neqns, ncoeffs;
    INFO0("Looking for baseflow coefficients");
    if (ESIO_SUCCESS == esio_plane_size(h, location_baseflow,
                                           &neqns, &ncoeffs)) {
        // Load baseflow coefficients
        x.resize(ncoeffs, neqns); //Using column-major storage

        esio_plane_establish(h,
            this->x.outerSize(), 0, this->x.outerSize(),
            this->x.innerSize(), 0, this->x.innerSize());
        esio_plane_read(h, location_baseflow, this->x.data(),
            this->x.outerStride(), this->x.innerStride());

        // Load baseflow coefficient base
//         char *name  = esio_string_get(h, location_baseflow_d, "coefficient_base");
//         this->dx_base = name;
//         free(name);

        WARN0("Assuming baseflow coefficients are for a polynomial base");
        this->x_base = "polynomial";
    }

    // Load baseflow derivative coefficients
    // Only proceed if largo_baseflow is present in the restart
    INFO0("Looking for baseflow derivative coefficients");
    if (ESIO_SUCCESS == esio_plane_size(h, location_baseflow_d,
                                           &neqns, &ncoeffs)) {
        // Load baseflow coefficients
        dx.resize(ncoeffs, neqns); //Using column-major storage

        esio_plane_establish(h,
            this->dx.outerSize(), 0, this->dx.outerSize(),
            this->dx.innerSize(), 0, this->dx.innerSize());
        esio_plane_read(h, location_baseflow_d, this->dx.data(),
            this->dx.outerStride(), this->dx.innerStride());

        // Load baseflow coefficient base
//         char *name  = esio_string_get(h, location_baseflow_d, "coefficient_base");
//         this->dx_base = name;
//         free(name);

        WARN0("Assuming baseflow derivative coefficients are for a polynomial base");
        this->dx_base = "polynomial";
    }
    return;
}

void
largo_definition::get_baseflow(
        const real_t      y,
        real_t *       base,
        real_t *     dybase,
        real_t *     dhbase)
{
   if (x.rows()) {
      real_t res[2];
      for (int j=0; j<x.cols()-1; j++){
          // Compute baseflow and y-derivative
          gsl_poly_eval_derivs(x.col(j).data(), x.rows(), y, &res[0], 2);
          base  [j] = res[0];
          dybase[j] = res[1];
      }
   } else {
       // Do nothing
       // Assume that baseflow arrays are initialized to zero
   }

   // Compute x-derivative of baseflow
   if (dx.rows()) {
      for (int j=0; j<dx.cols()-1; j++){
          // Compute x-derivative of baseflow
          gsl_poly_eval_derivs(dx.col(j).data(), dx.rows(), y, &dhbase[j], 1);
      }
   } else {
       // Do nothing
       // Assume that baseflow arrays are initialized to zero
   }
}


void
largo_definition::get_baseflow_pressure(
        const real_t      y,
        real_t &      Pbase,
        real_t &    dyPbase,
        real_t &    dhPbase)
{
   if (x.rows()) {
      real_t res[2];
      int j = x.cols()-1;
      // Compute baseflow and y-derivative
      gsl_poly_eval_derivs(x.col(j).data(), x.rows(), y, &res[0], 2);
      Pbase   = res[0];
      dyPbase = res[1];
   } else {
       // Do nothing
       // Assume that baseflow arrays are initialized to zero
   }

   // Compute x-derivative of baseflow
   if (dx.rows()) {
      int j = x.cols()-1;
      // Compute x-derivative of baseflow
      gsl_poly_eval_derivs(dx.col(j).data(), dx.rows(), y, &dhPbase, 1);
   } else {
       // Do nothing
       // Assume that baseflow arrays are initialized to zero
   }
}

} // namespace support

} // namespace suzerain
