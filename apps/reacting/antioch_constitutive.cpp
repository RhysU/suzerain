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
 * @copydoc antioch_constitutive.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#ifdef HAVE_ANTIOCH

#include "antioch_constitutive.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Class wrapping libantioch functionality for chemically reacting
 * flow
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

antioch_constitutive::antioch_constitutive()
    : species_names(0)
    , chem_input_file()
    , Le   (std::numeric_limits<real_t>::quiet_NaN())
    , alpha(std::numeric_limits<real_t>::quiet_NaN())
{
}

antioch_constitutive::antioch_constitutive(
    const std::vector<std::string>& species_names,
    const std::string& chem_input_file,
    const real_t Le,
    const real_t alpha)
    : species_names  (species_names)
    , chem_input_file(chem_input_file)
    , Le    (Le)
    , alpha (alpha)
{
}

antioch_constitutive::~antioch_constitutive()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_species_names[]   = "species";
static const char name_chem_input_file[] = "chemfile";
static const char name_Le[]              = "Le";
static const char name_alpha[]           = "alpha";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_species_names[]   = "List of species in mixture";
static const char desc_chem_input_file[] = "Chemistry input file";
static const char desc_Le[]              = "Lewis number";
static const char desc_alpha[]           = "Ratio of bulk to dynamic viscosity";

boost::program_options::options_description
antioch_constitutive::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;
    using std::vector;

    options_description retval("libantioch constitutive law parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // species
    retval.add_options()(name_species_names, 
                         value<vector<string> >(), 
                         desc_species_names);

    // chemistry input file
    retval.add_options()(name_chem_input_file,
                         value<string>(),
                         desc_chem_input_file);
    // Le
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Le, name_Le));
    if (!(boost::math::isnan)(Le)) {
        p->default_value(lexical_cast<string>(Le));
    }
    retval.add_options()(name_Le, p.release(), desc_Le);

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
antioch_constitutive::populate(
        const antioch_constitutive& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(Le);
    CALL_MAYBE_POPULATE(alpha);
#undef CALL_MAYBE_POPULATE
}

void
antioch_constitutive::override(
        const antioch_constitutive& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(Le);
    CALL_MAYBE_OVERRIDE(alpha);
#undef CALL_MAYBE_OVERRIDE
}

void
antioch_constitutive::save(
        const esio_handle h) const
{
    DEBUG0("Storing antioch_constitutive parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    // NOTE: approach similar to that taken for "metadata_generated"
    // in driver_base::save_metadata
    static const char acd[] = "antioch_constitutive_data";

    // FIXME: Figure out how to set this to the antioch revision number
    real_t antioch_rev = 0.0;
    esio_line_write(h, acd, &antioch_rev, 0, 
                    "Place to stash all antioch_constitutive related data.");

    // number of species
    int Ns = this->species_names.size();
    esio_attribute_write(h, acd, "Ns", Ns);

    // species names
    for (size_t i=0; i<species_names.size(); ++i)
    {
        std::string speci("Species_"+i);
        esio_string_set(h, acd, speci.c_str(), this->species_names[i].c_str());
    }

    // chemistry input file
    esio_string_set(h, acd, "chem_input_file", this->chem_input_file.c_str());    

    // Lewis number
    esio_attribute_write(h, acd, name_Le   , &this->Le);

    // alpha
    esio_attribute_write(h, acd, name_alpha, &this->alpha);
}

void
antioch_constitutive::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading antioch_constitutive parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    antioch_constitutive t;

    static const char acd[] = "antioch_constitutive_data";

    // number of species
    int Ns;
    esio_attribute_read(h, acd, "Ns", &Ns);

    t.species_names.reserve(Ns);

    // species names
    for (int i=0; i<Ns; ++i)
    {
        std::string speci("Species_"+i);
        char* tmp = esio_string_get(h, acd, speci.c_str());
        t.species_names.push_back(tmp);
        free(tmp);
    }

    // chemistry input file
    char* tmp = esio_string_get(h, acd, "chem_input_file");
    t.chem_input_file = tmp;
    free(tmp);

    // Lewis number
    esio_attribute_read(h, acd, name_Le   , &t.Le);

    // alpha
    esio_attribute_read(h, acd, name_alpha, &t.alpha);

    // Prefer this to incoming (handle explictly for string data)
    if (this->species_names.size() > 0)
        INFO0("Clutching onto " << name_species_names << " over incoming");
    else
        this->species_names = t.species_names;

    if (this->chem_input_file.size() > 0)
        INFO0("Clutching onto " << name_chem_input_file << " over incoming");
    else
        this->chem_input_file = t.chem_input_file;


    // Prefer this to incoming
    this->populate(t, verbose);  
}

// Evaluate: takes state and gives back everything we need from
// Cantera, including temp, pres, transport props, enthalpies, and
// reaction rates
void 
antioch_constitutive::evaluate (const real_t  e,
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
    FATAL0("antioch_constitutive::evaluate is not implemented yet!");
}

} // namespace reacting

} // namespace suzerain

#endif // HAVE_ANTIOCH
