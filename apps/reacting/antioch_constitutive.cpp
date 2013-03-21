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

#include <antioch/antioch_version.h>
#include <antioch/read_reaction_set_data_xml.h>

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
    : antioch_ver(Antioch::get_antioch_version())
    , species_names(0)
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
    : antioch_ver(Antioch::get_antioch_version())
    , species_names  (species_names)
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
                         value(&species_names), 
                         desc_species_names);

    // chemistry input file
    retval.add_options()(name_chem_input_file,
                         value(&chem_input_file),
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

    using boost::lexical_cast;

    DEBUG0("Storing antioch_constitutive parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    // NOTE: approach similar to that taken for "metadata_generated"
    // in driver_base::save_metadata
    static const char acd[] = "antioch_constitutive_data";

    // FIXME: version is written by save but not read by load
    esio_line_write(h, acd, &antioch_ver, 0, 
                    "Antioch version number and antioch_contitutive data.");

    // number of species
    int Ns = this->species_names.size();
    esio_attribute_write(h, acd, "Ns", Ns);

    // species names
    for (size_t i=0; i<species_names.size(); ++i)
    {
        std::string speci("Species_" + lexical_cast<std::string>(i));
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
    using boost::lexical_cast;

    DEBUG0("Loading antioch_constitutive parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    antioch_constitutive t;

    static const char acd[] = "antioch_constitutive_data";

    // antioch version
    esio_line_read(h, acd, &t.antioch_ver, 0);

    // Will never overwrite antioch version info
    if (t.antioch_ver != this->antioch_ver) {
        //... but warn if it doesn't match
        WARN0("Antioch version has changed.  Was " << t.antioch_ver 
              << " when written.  Is now " << this->antioch_ver);
    }

    // number of species
    int Ns;
    esio_attribute_read(h, acd, "Ns", &Ns);

    t.species_names.reserve(Ns);

    // species names
    for (int i=0; i<Ns; ++i)
    {
        std::string speci("Species_" + lexical_cast<std::string>(i));
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


void antioch_constitutive::init_antioch()
{
    WARN0("antioch_constitutive is not fully functional yet!");

    mixture    = make_shared<Antioch::ChemicalMixture<real_t> >(species_names);
    reactions  = make_shared<Antioch::ReactionSet<real_t> >(*mixture);
    cea_thermo = make_shared<Antioch::CEAThermodynamics<real_t> >(*mixture);
    sm_thermo  = make_shared<Antioch::StatMechThermodynamics<real_t> >(*mixture);


    if (this->Ns() > 1) {
        // FIXME: This will do for now, but need to think about what
        // happens in parallel to avoid all ranks trying to read this
        // file.
        Antioch::read_reaction_set_data_xml<real_t>(chem_input_file, 
                                                    false /* verbose */, 
                                                    *reactions);
        kinetics = make_shared<Antioch::KineticsEvaluator<real_t> >(*reactions);
    }

    // TODO: Add consistency asserts... everybody has same number of
    // species, number of reactions is sane, valid curve fits,
    // anything else?
}

// Evaluate everything
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
    //WARN0("antioch_constitutive::evaluate is not fully functional yet!");

    const real_t irho = 1.0/rho;

    // FIXME: Incoming vectors are currently raw pointers (result of
    // .data calls on Eigen objects), but Antioch wants std::vector.
    // I'll kludge it here to get the thing compiling.  Soon, roystgnr
    // should have Antioch taking Eigen::Array objects (see issue
    // #2799), which we can pass directly from
    // apply_navier_stokes_spatial_operator.  If not, I will refactor
    // apply_navier_stokes_spatial_operator to use Eigen::Map to wrap
    // std::vector, allowing easy passing to Antioch.
    //
    // Either way, this whole function will be refactored, but the
    // basic steps are as follows:

    const size_t Ns = this->Ns();

    // Mass fractions
    std::vector<real_t> Y(cs, cs+Ns);

    // Mixture gas constant
    const real_t R_mix = this->mixture->R(Y);

    std::vector<real_t> molar_densities(Ns,0.0);
    this->mixture->molar_densities(rho,Y,molar_densities);

    // Compute temperature from internal energy (assuming thermal equilibrium)
    const real_t re_internal = e - 0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
    T = this->sm_thermo->T_from_e_tot(irho*re_internal, Y);

    // Compute pressure: ideal gas law with mixture gas constant
    p = rho*R_mix*T;

    // Use CEA thermo to compute h_RT_minus_s_R for reaction
    // calculations 
    //
    // NOTE: This is how FIN-S does it, so we follow for
    // complete consistency.  But, it might make more sense to use
    // stat mech based thermo to get this info... I'm not sure.
    std::vector<real_t> h_RT_minus_s_R(Ns);
    typedef typename Antioch::CEAThermodynamics<real_t>::Cache<real_t> Cache;
    Cache cea_cache(T);
    this->cea_thermo->h_RT_minus_s_R(cea_cache,h_RT_minus_s_R);


    // protect against calling kinetics when there are no kinetics
    // TODO: Set up antioch to avoid this if (i.e., make call to kinetics ok)
    if (Ns>1) {
        // Species eqn source terms
        std::vector<real_t> omega_dot(Ns);
        this->kinetics->compute_mass_sources(T, rho, R_mix, Y, molar_densities, 
                                             h_RT_minus_s_R, omega_dot);
        for (size_t i=0; i<Ns; ++i) om[i] = omega_dot[i];
    } else {
        om[0] = 0.0;
    }
    
    // Species enthalpies (assuming thermal equilibrium)
    for (unsigned int i=0; i<Ns; ++i)
        hs[i] = this->sm_thermo->h_tot(i, T);

    // Transport
    // FIXME: Add antioch support

    // For now, hardcode some numbers for air
    mu = 1.716e-5 * pow( T/273.0, 0.6667 );

    real_t Cp = this->sm_thermo->cp(T, T, Y);

    kap = mu * Cp / 0.7;

    // Is this right?  Copied from FIN-S but looks like inverse if Le
    // to me.
    real_t D0 = this->Le*kap*irho/Cp;

    for (unsigned int i=0; i<Ns; ++i)
        Ds[i] = D0;

    // TODO: assert that om sums to zero
    // TODO: assert that T and p are positive
    // TODO: assert transport props are positive
}


real_t 
antioch_constitutive::e_from_T (const real_t  T,
                                const std::vector<real_t> mass_fractions) const
{
    return this->sm_thermo->e_tot(T, mass_fractions);
}


} // namespace reacting

} // namespace suzerain

#endif // HAVE_ANTIOCH
