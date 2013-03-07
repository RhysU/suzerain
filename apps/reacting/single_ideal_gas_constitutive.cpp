//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc single_ideal_gas_constitutive.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "single_ideal_gas_constitutive.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file 
 * Provides constitutive models for a single species ideal gas
 * with constant specific heats, constant Prandtl number, and power
 * law viscosity.
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

single_ideal_gas_constitutive::single_ideal_gas_constitutive()
    : Cp    (std::numeric_limits<real_t>::quiet_NaN())
    , Cv    (std::numeric_limits<real_t>::quiet_NaN())
    , Pr    (std::numeric_limits<real_t>::quiet_NaN())
    , T0    (std::numeric_limits<real_t>::quiet_NaN())
    , mu0   (std::numeric_limits<real_t>::quiet_NaN())
    , beta  (std::numeric_limits<real_t>::quiet_NaN())
    , alpha (std::numeric_limits<real_t>::quiet_NaN())
{
}

single_ideal_gas_constitutive::single_ideal_gas_constitutive(
        const real_t Cp,
	const real_t Cv,
	const real_t Pr,
	const real_t T0,
	const real_t mu0,
	const real_t beta,
	const real_t alpha)
    : Cp    (Cp)
    , Cv    (Cv)
    , Pr    (Pr)
    , T0    (T0)
    , mu0   (mu0)
    , beta  (beta)
    , alpha (alpha)
{
}

single_ideal_gas_constitutive::~single_ideal_gas_constitutive()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_Cp[]    = "Cp";
static const char name_Cv[]    = "Cv";
static const char name_Pr[]    = "Pr";
static const char name_T0[]    = "T0";
static const char name_mu0[]   = "mu0";
static const char name_beta[]  = "beta";
static const char name_alpha[] = "alpha";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_Cp[]    = "Specific heat at constant pressure";
static const char desc_Cv[]    = "Specific heat at constant volume";
static const char desc_Pr[]    = "Prandtl number";
static const char desc_T0[]    = "Reference temperature for viscosity power law";
static const char desc_mu0[]   = "Reference viscosity for viscosity power law";
static const char desc_beta[]  = "Exponent for viscosity power law";
static const char desc_alpha[] = "Ratio of bulk to dynamic viscosity";

boost::program_options::options_description
single_ideal_gas_constitutive::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("Single species ideal gas constitutive law parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // Cp
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Cp, name_Cp));
    if (!(boost::math::isnan)(Cp)) {
        p->default_value(lexical_cast<string>(Cp));
    }
    retval.add_options()(name_Cp, p.release(), desc_Cp);

    // Cv
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Cv, name_Cv));
    if (!(boost::math::isnan)(Cv)) {
        p->default_value(lexical_cast<string>(Cv));
    }
    retval.add_options()(name_Cv, p.release(), desc_Cv);

    // Pr
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Pr, name_Pr));
    if (!(boost::math::isnan)(Pr)) {
        p->default_value(lexical_cast<string>(Pr));
    }
    retval.add_options()(name_Pr, p.release(), desc_Pr);

    // T0
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &T0, name_T0));
    if (!(boost::math::isnan)(T0)) {
        p->default_value(lexical_cast<string>(T0));
    }
    retval.add_options()(name_T0, p.release(), desc_T0);

    // mu0
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &mu0, name_mu0));
    if (!(boost::math::isnan)(mu0)) {
        p->default_value(lexical_cast<string>(mu0));
    }
    retval.add_options()(name_mu0, p.release(), desc_mu0);

    // beta
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &beta, name_beta));
    if (!(boost::math::isnan)(beta)) {
        p->default_value(lexical_cast<string>(beta));
    }
    retval.add_options()(name_beta, p.release(), desc_beta);

    // alpha
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &alpha, name_alpha));
    if (!(boost::math::isnan)(alpha)) {
        p->default_value(lexical_cast<string>(alpha));
    }
    retval.add_options()(name_alpha, p.release(), desc_alpha);


    return retval;
}

void
single_ideal_gas_constitutive::populate(
        const single_ideal_gas_constitutive& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(Cp);
    CALL_MAYBE_POPULATE(Cv);
    CALL_MAYBE_POPULATE(Pr);
    CALL_MAYBE_POPULATE(T0);
    CALL_MAYBE_POPULATE(mu0);
    CALL_MAYBE_POPULATE(beta);
    CALL_MAYBE_POPULATE(alpha);
#undef CALL_MAYBE_POPULATE
}

void
single_ideal_gas_constitutive::override(
        const single_ideal_gas_constitutive& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(Cp);
    CALL_MAYBE_OVERRIDE(Cv);
    CALL_MAYBE_OVERRIDE(Pr);
    CALL_MAYBE_OVERRIDE(T0);
    CALL_MAYBE_OVERRIDE(mu0);
    CALL_MAYBE_OVERRIDE(beta);
    CALL_MAYBE_OVERRIDE(alpha);
#undef CALL_MAYBE_OVERRIDE
}

void
single_ideal_gas_constitutive::save(
        const esio_handle h) const
{
    DEBUG0("Storing single_ideal_gas_constitutive parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_Cp,    &this->Cp,    0, desc_Cp   );
    esio_line_write(h, name_Cv,    &this->Cv,    0, desc_Cv   );
    esio_line_write(h, name_Pr,    &this->Pr,    0, desc_Pr   );
    esio_line_write(h, name_T0,    &this->T0,    0, desc_T0   );
    esio_line_write(h, name_mu0,   &this->mu0,   0, desc_mu0  );
    esio_line_write(h, name_beta,  &this->beta,  0, desc_beta );
    esio_line_write(h, name_alpha, &this->alpha, 0, desc_alpha);
}

void
single_ideal_gas_constitutive::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading single_ideal_gas_constitutive parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    single_ideal_gas_constitutive t;
    esio_line_read(h, name_Cp,    &t.Cp,    0);
    esio_line_read(h, name_Cv,    &t.Cv,    0);
    esio_line_read(h, name_Pr,    &t.Pr,    0);
    esio_line_read(h, name_T0,    &t.T0,    0);
    esio_line_read(h, name_mu0,   &t.mu0,   0);
    esio_line_read(h, name_beta,  &t.beta,  0);
    esio_line_read(h, name_alpha, &t.alpha, 0);
    this->populate(t, verbose);  // Prefer this to incoming
}

// Evaluate: takes state and gives back everything we need from
// Cantera, including temp, pres, transport props, enthalpies, and
// reaction rates
void 
single_ideal_gas_constitutive::evaluate (const real_t  e,
					 const real_t* m,
					 const real_t  rho,
					 const real_t* species,
					 const real_t* cs,
					 real_t& T,
					 real_t& p,
					 real_t* Ds,
					 real_t& mu,
					 real_t& kap,
					 real_t* hs,
					 real_t* om) const
{
    // NOTE: Really only have e, m, rho input and T, p, mu, kap
    // output.  Everything else is used for the multispecies case.

    // convenience (some of which could be moved out to save a couple flops)
    const real_t irho = 1.0/rho;
    const real_t gam = Cp/Cv;
    const real_t gmi = gam-1.0;
    const real_t Rgas = Cp-Cv;

    p = gmi * (e - 0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]));

    T = irho * p / Rgas;

    mu = mu0 * pow( T/T0, beta );

    kap = mu * Cp / Pr;

    Ds = NULL;
    hs = NULL;
    om = NULL;
}

} // namespace reacting

} // namespace suzerain
