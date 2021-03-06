//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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
 * @copydoc definition_largo.hpp
 */

#include <suzerain/support/definition_largo.hpp>

#include <esio/error.h>
#include <esio/esio.h>

#include <gsl/gsl_poly.h>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/largo_formulation.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain
{

static void parse_nonnegative(const std::string& s, real_t* t, const char* n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_nonnegative(v, n);
    *t = v;
}

namespace support
{

static void parse_formulation(const std::string& s, largo_formulation* t)
{
    *t = largo_formulation::lookup(s);
}

definition_largo::definition_largo()
    : specification_largo()
{
}

// Strings used in options_description and populate/override/save/load.
static const char name_formulation[]          = "formulation";
static const char name_grdelta[]              = "grdelta";
static const char name_ignore_fluctuations[]  = "ignore_fluctuations";
static const char name_ignore_gramp_mean[]    = "ignore_gramp_mean";
static const char name_ignore_gramp_rms[]     = "ignore_gramp_rms";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_formulation[]
    = "Name of the slow growth formulation";
static const char desc_grdelta[]
    = "Growth rate of reference thickness (Delta)";
static const char desc_ignore_fluctuations[]
    = "Should slow growth forcing ignore all fluctuations?"
    "  Meant for debugging and not preserved across restarts.";
static const char desc_ignore_gramp_mean[]
    = "Should slow growth forcing ignore all mean growth rates?"
    "  Meant for debugging and not preserved across restarts.";
static const char desc_ignore_gramp_rms[]
    = "Should slow growth forcing ignore all RMS growth rates?"
    "  Meant for debugging and not preserved across restarts.";

boost::program_options::options_description
definition_largo::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
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
    ("largo_ignore_fluctuations",
     bool_switch(&ignore_fluctuations)->default_value(ignore_fluctuations),
     desc_ignore_fluctuations)
    ("largo_ignore_gramp_mean",
     bool_switch(&ignore_gramp_mean)->default_value(ignore_gramp_mean),
     desc_ignore_gramp_mean)
    ("largo_ignore_gramp_rms",
     bool_switch(&ignore_gramp_rms)->default_value(ignore_gramp_rms),
     desc_ignore_gramp_rms)
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

// For maybe_XXX_impl to determine a \c largo_formulation is a defaulted
static bool default_value_formulation(const largo_formulation& v)
{
    return !v.enabled();
}

static bool maybe_populate(const char*              name,
                           const char*              description,
                           largo_formulation&       destination,
                           const largo_formulation& source,
                           const bool               verbose)
{
    return internal::maybe_populate_impl(
               name, description, destination, source,
               verbose, &default_value_formulation);
}

void
definition_largo::populate(
    const definition_largo& that,
    const bool verbose)
{
    using support::maybe_populate;
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(formulation);
    CALL_MAYBE_POPULATE(grdelta);
    CALL_MAYBE_POPULATE(ignore_fluctuations);
    CALL_MAYBE_POPULATE(ignore_gramp_mean);
    CALL_MAYBE_POPULATE(ignore_gramp_rms);
#undef CALL_MAYBE_POPULATE
    if (!this->baseflow) {
        this->baseflow = that.baseflow;
    }

    if ((this->gramp_mean.size() == 0)) {
        TRACE0("Populating growth rate amplitude mean");
        this->gramp_mean = that.gramp_mean;
    }
    if ((this->gramp_rms.size() == 0)) {
        TRACE0("Populating growth rate amplitude rms");
        this->gramp_rms = that.gramp_rms;
    }
}

static bool maybe_override(const char*              name,
                           const char*              description,
                           largo_formulation&       destination,
                           const largo_formulation& source,
                           const bool               verbose)
{
    return internal::maybe_override_impl(
               name, description, destination, source,
               verbose, &default_value_formulation);
}

void
definition_largo::override(
    const definition_largo& that,
    const bool verbose)
{
    using support::maybe_override;
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(formulation);
    CALL_MAYBE_OVERRIDE(grdelta);
    CALL_MAYBE_OVERRIDE(ignore_fluctuations);
    CALL_MAYBE_OVERRIDE(ignore_gramp_mean);
    CALL_MAYBE_OVERRIDE(ignore_gramp_rms);
#undef CALL_MAYBE_OVERRIDE
    if (that.baseflow) {
        this->baseflow = that.baseflow;
    }

    if ((that.gramp_mean.size() > 0)) {
        TRACE0("Overriding growth rate amplitude mean");
        this->gramp_mean = that.gramp_mean;
    }
    if ((that.gramp_rms.size() > 0)) {
        TRACE0("Overriding growth rate amplitude rms");
        this->gramp_rms = that.gramp_rms;
    }
}

static const char attr_base[]            = "coefficient_base";
static const char attr_formulation[]     = "formulation";
static const char location_gramp_mean[]  = "largo_gramp_mean";
static const char location_gramp_rms[]   = "largo_gramp_rms";
static const char location_baseflow_dx[] = "largo_baseflow_dx";
static const char location_baseflow_dy[] = "largo_baseflow_dy";
static const char location_baseflow[]    = "largo_baseflow";
static const char location[]             = "largo";

// std::strings permit easy comparison for definition_largo::load()
static const std::string type_polynomial("polynomial");
static const std::string type_uniform   ("uniform");
static const std::string type_map       ("map");

void
definition_largo::save(
    const esio_handle h) const
{
    DEBUG0("Storing definition_largo parameters");

    if (!formulation.enabled()) {
        return;    // Shortcircuit on no formulation
    }

    // Write out the "container" holding all other settings
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
                    "Is a largo-based slow growth formulation in use?");

    // Write out the formulation name
    esio_string_set(h, location, attr_formulation, formulation.name().c_str());

    // scalars
    if (formulation == largo_formulation::disable) {
        return; // Nothing else to save, STOP!
    } else if (formulation == largo_formulation::temporal) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::spatial) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::temporal_tensor_consistent) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::spatiotemporal) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::temporal_consistent) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else if (formulation == largo_formulation::spatiotemporal_consistent) {
        esio_attribute_write(h, location, name_grdelta, &grdelta);
    } else {
        FATAL0("Attempt to save unknown largo_formulation");
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

    // My most sincere apologies for using dynamic_cast to dispatch (#2503)
    // A simple OOPish design would push ESIO concerns into baseflow.hpp.
    // A proper OOPish design would introduce the visitor pattern.
    if (!baseflow) {

        // No baseflow to save

    } else if (dynamic_cast<baseflow_uniform*>(baseflow.get())) {

        baseflow_uniform* const b
            = dynamic_cast<baseflow_uniform*>(baseflow.get());

        // Write baseflow field values (only master process)
        const MatrixXXr& x  = b->x;
        esio_plane_establish(
                h,
                x.outerSize(), 0, procid == 0 ? x.outerSize() : 0,
                x.innerSize(), 0, procid == 0 ? x.innerSize() : 0);
        esio_plane_write(
                h, location_baseflow, x.data(),
                x.outerStride(), x.innerStride(),
                "Baseflow field values");
        esio_string_set(
                h, location_baseflow,
                attr_base, type_uniform.c_str());

    } else if (dynamic_cast<baseflow_polynomial*>(baseflow.get())) {

        baseflow_polynomial* const b
            = dynamic_cast<baseflow_polynomial*>(baseflow.get());

        // Write baseflow field coefficients (only master process)
        const MatrixXXr& x  = b->x;
        esio_plane_establish(
                h,
                x.outerSize(), 0, procid == 0 ? x.outerSize() : 0,
                x.innerSize(), 0, procid == 0 ? x.innerSize() : 0);
        esio_plane_write(
                h, location_baseflow, x.data(),
                x.outerStride(), x.innerStride(),
                "Baseflow field coefficients");
        esio_string_set(
                h, location_baseflow,
                attr_base, type_polynomial.c_str());

        // Write baseflow derivative coefficients (only master process)
        const MatrixXXr& dx = b->dx;
        esio_plane_establish(
                h,
                dx.outerSize(), 0, procid == 0 ?  dx.outerSize() : 0,
                dx.innerSize(), 0, procid == 0 ?  dx.innerSize() : 0);
        esio_plane_write(
                h, location_baseflow_dx, dx.data(),
                dx.outerStride(), dx.innerStride(),
                "Baseflow streamwise derivative coefficients");
        esio_string_set(
                h, location_baseflow_dx,
                attr_base, type_polynomial.c_str());

    } else if (dynamic_cast<baseflow_map*>(baseflow.get())) {

        baseflow_map* const b = dynamic_cast<baseflow_map*>(baseflow.get());

        // Write baseflow pointwise locations and state (only master process)
        MatrixXAr tmp(b->table.size(), /*y + state + p + u + v + w*/ 10);
        esio_plane_establish(
                h,
                tmp.outerSize(), 0, procid == 0 ? tmp.outerSize() : 0,
                tmp.innerSize(), 0, procid == 0 ? tmp.innerSize() : 0);

        // Pack pointwise data and save it to location_baseflow
        // computing primitive velocity on the fly.
        // Changes here require syncing the logic within load(...).
        baseflow_map::table_type::const_iterator it = b->table.begin();
        for (int i = 0; i < tmp.rows(); ++i) {
            const largo_state& base = it->second.base;
            tmp.row(i) << it->first,
                          base.e,
                          base.mx,
                          base.my,
                          base.mz,
                          base.rho,
                          base.p,
                          base.mx / base.rho,
                          base.my / base.rho,
                          base.mz / base.rho;
            ++it;
        }
        esio_plane_write(
                h, location_baseflow, tmp.data(),
                tmp.outerStride(), tmp.innerStride(),
                "Non-normative baseflow pointwise state:"
                " y, rho_E, rho_u, rho_v, rho_w, rho, p, u, v, w");
        esio_string_set(
                h, location_baseflow,
                attr_base, type_map.c_str());

        // Pack pointwise x derivatives and save to location_baseflow_dx
        // computing primitive velocity derivatives on the fly.
        // Changes here require syncing the logic within load(...).
        it = b->table.begin();
        for (int i = 0; i < tmp.rows(); ++i) {
            const largo_state& base   = it->second.base;
            const largo_state& dxbase = it->second.dxbase;
            tmp.row(i) << it->first,
                          dxbase.e,
                          dxbase.mx,
                          dxbase.my,
                          dxbase.mz,
                          dxbase.rho,
                          dxbase.p,
                          (dxbase.mx - dxbase.rho*base.mx/base.rho)/base.rho,
                          (dxbase.my - dxbase.rho*base.my/base.rho)/base.rho,
                          (dxbase.mz - dxbase.rho*base.mz/base.rho)/base.rho;
            ++it;
        }
        esio_plane_write(
                h, location_baseflow_dx, tmp.data(),
                tmp.outerStride(), tmp.innerStride(),
                "Non-normative baseflow pointwise x derivatives:"
                " y, rho_E, rho_u, rho_v, rho_w, rho, p, u, v, w");
        esio_string_set(
                h, location_baseflow_dx,
                attr_base, type_map.c_str());

        // Pack pointwise y derivatives and save to location_baseflow_dy
        // computing primitive velocity derivatives on the fly.
        // Changes here require syncing the logic within load(...).
        it = b->table.begin();
        for (int i = 0; i < tmp.rows(); ++i) {
            const largo_state& base   = it->second.base;
            const largo_state& dybase = it->second.dybase;
            tmp.row(i) << it->first,
                          dybase.e,
                          dybase.mx,
                          dybase.my,
                          dybase.mz,
                          dybase.rho,
                          dybase.p,
                          (dybase.mx - dybase.rho*base.mx/base.rho)/base.rho,
                          (dybase.my - dybase.rho*base.my/base.rho)/base.rho,
                          (dybase.mz - dybase.rho*base.mz/base.rho)/base.rho;
            ++it;
        }
        esio_plane_write(
                h, location_baseflow_dy, tmp.data(),
                tmp.outerStride(), tmp.innerStride(),
                "Non-normative baseflow pointwise y derivatives:"
                " y, rho_E, rho_u, rho_v, rho_w, rho, p, u, v, w");
        esio_string_set(
                h, location_baseflow_dy,
                attr_base, type_map.c_str());

    } else {

        FATAL0("Attempt to save unknown baseflow description");
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    }

    if (!this->gramp_mean.size() == 0) {
        TRACE0("Saving growth rate amplitude mean");

        // Write amplitude growth rate values for the mean (only master process)
        const std::vector<real_t>& gm  = this->gramp_mean;
        esio_line_establish(
                h,
                gm.size(), 0, procid == 0 ? gm.size() : 0);
        esio_line_write(
                h, location_gramp_mean, gm.data(), 0,
                "Amplitude growth rates for the mean values"
                " stored per the formulation employed");
    }

    if (!this->gramp_rms.size() == 0) {
        TRACE0("Saving growth rate amplitude rms");

        // Write amplitude growth rate values for the rms (only master process)
        const std::vector<real_t>& gr  = this->gramp_rms;
        esio_line_establish(
                h,
                gr.size(), 0, procid == 0 ? gr.size() : 0);
        esio_line_write(
                h, location_gramp_rms, gr.data(), 0,
                "Amplitude growth rates for the rms values"
                " stored per the formulation employed");
    }

    TRACE0("Runtime flag " << name_ignore_fluctuations << "is never saved");
    TRACE0("Runtime flag " << name_ignore_gramp_mean   << "is never saved");
    TRACE0("Runtime flag " << name_ignore_gramp_rms    << "is never saved");
}

void
definition_largo::load(
    const esio_handle h,
    const bool verbose)
{
    using std::free;

    definition_largo t;
    assert(t.formulation == largo_formulation::disable);

    // Only proceed if a largo definition is active in the restart
    int in_use = 0;
    esio_line_establish(h, 1, 0, 1); // All ranks load any data
    if (ESIO_NOTFOUND != esio_line_size(h, location, NULL)) {
        esio_line_read(h, location, &in_use, 0);
    }
    if (!in_use) {
        return;
    }

    DEBUG0("Loading definition_largo parameters");

    // Load formulation name and look it up in the static instance map
    {
        char* name = esio_string_get(h, location, attr_formulation);
        t.formulation = largo_formulation::lookup(name);
        free(name);
    }

    if (t.formulation == largo_formulation::disable) {
        // Nothing else to load, STOP!
        return this->populate(t, verbose);
    } else if (t.formulation == largo_formulation::temporal) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::spatial) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::temporal_tensor_consistent) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::spatiotemporal) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::temporal_consistent) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else if (t.formulation == largo_formulation::spatiotemporal_consistent) {
        esio_attribute_read(h, location, name_grdelta, &t.grdelta);
    } else {
        FATAL0("Attempt to load unknown largo_formulation");
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    }

    // Strategy for loading definition_largo baseflow information:
    // Load first.  Sort it out second.
    // Otherwise interspersing IO with logic is a gnarly mess.
    int neqns, ncoeffs;

    DEBUG0("Probing for baseflow state");
    shared_ptr<char> base_x;
    MatrixXXr x;
    if (ESIO_SUCCESS == esio_plane_size(h, location_baseflow,
                                        &neqns, &ncoeffs)) {
        base_x.reset(esio_string_get(h, location_baseflow, attr_base),
                     free);
        DEBUG0("Loading baseflow information of type: "
               << (base_x ? base_x.get() : "NULL"));
        x.resize(ncoeffs, neqns);
        esio_plane_establish(
                h,
                x.outerSize(), 0, x.outerSize(),
                x.innerSize(), 0, x.innerSize());
        esio_plane_read(
                h, location_baseflow, x.data(),
                x.outerStride(), x.innerStride());
    }

    DEBUG0("Probing for baseflow streamwise derivatives");
    shared_ptr<char> base_dx;
    MatrixXXr dx;
    if (ESIO_SUCCESS == esio_plane_size(h, location_baseflow_dx,
                                        &neqns, &ncoeffs)) {
        base_dx.reset(esio_string_get(h, location_baseflow_dx, attr_base),
                      free);
        DEBUG0("Loading baseflow streamwise derivative information of type: "
               << (base_dx ? base_dx.get() : "NULL"));
        dx.resize(ncoeffs, neqns);
        esio_plane_establish(
                h,
                dx.outerSize(), 0, dx.outerSize(),
                dx.innerSize(), 0, dx.innerSize());
        esio_plane_read(
                h, location_baseflow_dx, dx.data(),
                dx.outerStride(), dx.innerStride());
    } else if (ESIO_SUCCESS == esio_plane_size(  // Pre-r42574 legacy handling
                h, "largo_baseflow_derivative", &neqns, &ncoeffs)) {
        WARN0("Detected legacy layout for '/largo_baseflow_derivative'");
        WARN0("Update code and scripts to use '/'" << location_baseflow_dx);
        base_dx.reset(esio_string_get(h, "largo_baseflow_derivative",
                      attr_base), free);
        DEBUG0("Loading baseflow streamwise derivative information of type: "
               << (base_dx ? base_dx.get() : "NULL"));
        dx.resize(ncoeffs, neqns);
        esio_plane_establish(
                h,
                dx.outerSize(), 0, dx.outerSize(),
                dx.innerSize(), 0, dx.innerSize());
        esio_plane_read(
                h, "largo_baseflow_derivative", dx.data(),
                dx.outerStride(), dx.innerStride());
    }

    DEBUG0("Probing for baseflow wall-normal derivatives");
    shared_ptr<char> base_dy;
    MatrixXXr dy;
    if (ESIO_SUCCESS == esio_plane_size(h, location_baseflow_dy,
                                        &neqns, &ncoeffs)) {
        base_dy.reset(esio_string_get(h, location_baseflow_dy, attr_base),
                      free);
        DEBUG0("Loading baseflow wall-normal derivative information of type: "
               << (base_dy ? base_dy.get() : "NULL"));
        dy.resize(ncoeffs, neqns);
        esio_plane_establish(
                h,
                dy.outerSize(), 0, dy.outerSize(),
                dy.innerSize(), 0, dy.innerSize());
        esio_plane_read(
                h, location_baseflow_dy, dy.data(),
                dy.outerStride(), dy.innerStride());
    }

    using boost::starts_with;
    if (!base_x && !base_dx && !base_dy) {

        // No baseflow to load

    } else if (    (base_x  && type_map == base_x .get())
                || (base_dx && type_map == base_dx.get())
                || (base_dy && type_map == base_dy.get())) {

        DEBUG0("Preparing table-driven baseflow description");
        SUZERAIN_ENSURE(x.rows() == dx.rows());
        SUZERAIN_ENSURE(x.rows() == dy.rows());
        SUZERAIN_ENSURE(x.cols() == dx.cols());
        SUZERAIN_ENSURE(x.cols() == dy.cols());

        // See save(...) method for the origin of the magic numbers below
        shared_ptr<baseflow_map> const p = make_shared<baseflow_map>();
        for (int i = 0; i < x.rows(); ++i) {
            baseflow_map::row& row = p->table[x(i, 0)];
            row.  base.e   =  x(i, 1);
            row.  base.mx  =  x(i, 2);
            row.  base.my  =  x(i, 3);
            row.  base.mz  =  x(i, 4);
            row.  base.rho =  x(i, 5);
            row.  base.p   =  x(i, 6);
            row.dxbase.e   = dx(i, 1);
            row.dxbase.mx  = dx(i, 2);
            row.dxbase.my  = dx(i, 3);
            row.dxbase.mz  = dx(i, 4);
            row.dxbase.rho = dx(i, 5);
            row.dxbase.p   = dx(i, 6);
            row.dybase.e   = dy(i, 1);
            row.dybase.mx  = dy(i, 2);
            row.dybase.my  = dy(i, 3);
            row.dybase.mz  = dy(i, 4);
            row.dybase.rho = dy(i, 5);
            row.dybase.p   = dy(i, 6);
        }
        t.baseflow = p;

    } else if (    (base_x  && starts_with(base_x .get(), type_polynomial))
                && (base_dx && starts_with(base_dx.get(), type_polynomial)) ){

        // Fuzziness added to the matching to wallpaper over Python/ESIO hiccup
        if (type_polynomial != base_x.get()) {
            WARN0("Polynomial baseflow from starts_with(base_x, 'polynomial')");
        }
        if (type_polynomial != base_dx.get()) {
            WARN0("Polynomial baseflow from starts_with(base_x, 'polynomial')");
        }

        DEBUG0("Preparing polynomial-based baseflow description");
        shared_ptr<baseflow_polynomial> p = make_shared<baseflow_polynomial>();
        p->x       = x;
        p->dx      = dx;
        t.baseflow = p;

    } else if (base_x && type_uniform == base_x.get()) {

        DEBUG0("Preparing uniform baseflow description");
        shared_ptr<baseflow_uniform> p = make_shared<baseflow_uniform>();
        p->x       = x;
        t.baseflow = p;

    } else {

        FATAL0("Attempt to load unknown baseflow description");
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    }

    TRACE0("Probing for amplitude growth rates for the mean");
    if (ESIO_SUCCESS == esio_line_size(h, location_gramp_mean,
                                        &neqns)) {
        TRACE0("Loading growth rate amplitude coefficients for the mean");
        t.gramp_mean.resize(neqns);
        esio_line_establish(
                h, t.gramp_mean.size(), 0, t.gramp_mean.size());
        esio_line_read(
                h, location_gramp_mean, t.gramp_mean.data(), 0);
    }

    TRACE0("Probing for amplitude growth rates for the RMS");
    if (ESIO_SUCCESS == esio_line_size(h, location_gramp_rms,
                                        &neqns)) {
        TRACE0("Loading growth rate amplitude coefficients for RMS");
        t.gramp_rms.resize(neqns);
        esio_line_establish(
                h, t.gramp_rms.size(), 0, t.gramp_rms.size());
        esio_line_read(
                h, location_gramp_rms, t.gramp_rms.data(), 0);
    }

    TRACE0("Runtime flag " << name_ignore_fluctuations << "is never loaded");
    TRACE0("Runtime flag " << name_ignore_gramp_mean   << "is never loaded");
    TRACE0("Runtime flag " << name_ignore_gramp_rms    << "is never loaded");

    return this->populate(t, verbose);
}

} // namespace support

} // namespace suzerain
