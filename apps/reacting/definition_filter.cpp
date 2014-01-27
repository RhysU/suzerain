//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc definition_filter.hpp
 */

#include "definition_filter.hpp"

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem defintion for filter
 * source term, e.g., strength coefficient.
 */

namespace suzerain {

static void parse_positive(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_positive(v, n);
    *t = v;
}

static void parse_nonnegative(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_nonnegative(v, n);
    *t = v;
}

namespace reacting {

filter_definition::filter_definition()
    : filter_phi(std::numeric_limits<real_t>::quiet_NaN())
{
}

filter_definition::filter_definition(
        const real_t filter_phi)
    : filter_phi(filter_phi)
{
}

// Strings used in options_description and populate/override/save/load.
static const char name_filter_phi[] = "filter_phi";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_filter_phi[] = "Filter source strength";

boost::program_options::options_description
filter_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("filter source parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // filter_phi
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &filter_phi, name_filter_phi));
    if (!(boost::math::isnan)(filter_phi)) {
        p->default_value(lexical_cast<string>(filter_phi));
    }
    retval.add_options()(name_filter_phi, p.release(), desc_filter_phi);

    return retval;
}

void
filter_definition::populate(
        const filter_definition& that,
        const bool verbose)
{
    using support::maybe_populate;
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(filter_phi);
#undef CALL_MAYBE_POPULATE
}

void
filter_definition::override(
        const filter_definition& that,
        const bool verbose)
{
    using support::maybe_override;
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(filter_phi);
#undef CALL_MAYBE_OVERRIDE
}

void
filter_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing filter_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_filter_phi, &this->filter_phi, 0, desc_filter_phi);
}

void
filter_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading filter_definition parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    filter_definition t;
    esio_line_read(h, name_filter_phi, &t.filter_phi, 0);
    this->populate(t, verbose);  // Prefer this to incoming
}

// FIXME Default method
static const suzerain_filterop_method method
    = SUZERAIN_FILTEROP_COOKCABOT2005;

// FIXME Use default parameters for the method
static const double *r_method_params
    = NULL;
//static const double r_method_params[1] = {6.55026621150074e-01};

// FIXME Use default parameters for the method
static const complex_double *z_method_params
    = NULL;

// FIXME Default boundary treatment on first edge
static const suzerain_filterop_boundary_treatment b_first
    = SUZERAIN_FILTEROP_BOUNDARY_NOFILTER;

// FIXME Default boundary treatment on second edge
static const suzerain_filterop_boundary_treatment b_last
    = SUZERAIN_FILTEROP_BOUNDARY_NOFILTER;

int
filter_definition::prepare_real(
        const int n) const
{
    r.reset(suzerain_filterop_alloc(n, method, r_method_params,
                                    b_first, b_last),
            suzerain_filterop_free);
    return suzerain_filterop_factorize(r.get());
}

int
filter_definition::prepare_complex(
        const int n) const
{
    z.reset(suzerain_filteropz_alloc(n, method, z_method_params,
                                     b_first, b_last),
            suzerain_filteropz_free);
    return suzerain_filteropz_factorize(z.get());
}

void
filter_definition::reset()
{
    r.reset();
    z.reset();
}

void filter_definition::save_filteropz(const esio_handle h,
                                       const int n,
                                       const char *location)
{
    // Only save filter op if hasn't been prepared yet
    if (SUZERAIN_UNLIKELY(!r)) {
        DEBUG0("Storing filter operators");

        // Create
        r.reset(suzerain_filterop_alloc(n, method, r_method_params,
                                        b_first, b_last),
                suzerain_filterop_free);

        // Write

        // Only root writes data
        int procid;
        esio_handle_comm_rank(h, &procid);

        char name[8] = {};
        char comment[127] = {};

        // A transpose
        snprintf(name, sizeof(name), "AT");
        snprintf(comment, sizeof(comment),
                 "Filter operator trans(A(i,j)) = AT[j,ku+i-j] for"
                 " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)");
        const int lda = r->ldat; //r->kuat + 1 + r->klat;// cop->ku(k) + 1 + cop->kl(k);
        //const int lda = r->kuat + 1 + r->klat;// cop->ku(k) + 1 + cop->kl(k);
        esio_plane_establish(h,
                             r->n, 0, (procid == 0 ? r->n : 0),
                             lda,  0, (procid == 0 ? lda  : 0));
        //esio_plane_write(h, name, r->A_T+r->klat, 0, 0, comment);
        esio_plane_write(h, name, r->A_T, 0, 0, comment);
        esio_attribute_write(h, name, "kl", r->klat);
        esio_attribute_write(h, name, "ku", r->kuat);
        esio_attribute_write(h, name, "m",  r->n);
        esio_attribute_write(h, name, "n",  r->n);

        // B transpose
        snprintf(name, sizeof(name), "BT");
        snprintf(comment, sizeof(comment),
                 "Filter operator trans(B(i,j)) = BT[j,ku+i-j] for"
                 " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)");
        const int ldb = r->kubt + 1 + r->klbt;// cop->ku(k) + 1 + cop->kl(k);
        esio_plane_establish(h,
                             r->n, 0, (procid == 0 ? r->n : 0),
                             ldb,  0, (procid == 0 ? ldb  : 0));
        esio_plane_write(h, name, r->B_T, 0, 0, comment);
        esio_attribute_write(h, name, "kl", r->klbt);
        esio_attribute_write(h, name, "ku", r->kubt);
        esio_attribute_write(h, name, "m",  r->n);
        esio_attribute_write(h, name, "n",  r->n);

        // Destroy
        r.reset();

    }

    // Otherwise, don't

}

} // namespace reacting

} // namespace suzerain
