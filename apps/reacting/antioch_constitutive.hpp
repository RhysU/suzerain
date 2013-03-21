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

#ifndef SUZERAIN_ANTIOCH_CONSTITUTIVE_HPP
#define SUZERAIN_ANTIOCH_CONSTITUTIVE_HPP

/** @file
 * Class wrapping libantioch functionality for chemically reacting
 * flow
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#ifdef HAVE_ANTIOCH // only makes sense when antioch is available

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

#include "antioch/vector_utils.h" // should this be necessary?
#include <antioch/chemical_mixture.h>
#include <antioch/reaction_set.h>
#include <antioch/cea_thermo.h>
#include <antioch/stat_mech_thermo.h>
#include <antioch/kinetics_evaluator.h>


namespace suzerain {

namespace reacting {

/**
 * Provides constitutive models for a mixture of ideal gasses
 * leveraging libantioch
 */
class antioch_constitutive : public support::definition_base,
                             public boost::noncopyable
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    antioch_constitutive();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param Ns   Number of species
     * @param bulk_rho_u Bulk streamwise momentum target.
     */
    antioch_constitutive(const std::vector<std::string>& species_names,
                         const std::string& chem_input_file,
                         const real_t Le,
                         const real_t alpha);


    /** Virtual destructor to permit use as a base class */
    virtual ~antioch_constitutive();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const antioch_constitutive& that,
            const bool verbose = false);

    /**
     * Override members in \c this with non-NaN values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when an override occurs?
     */
    virtual void override(
            const antioch_constitutive& that,
            const bool verbose = false);

    /**
     * Save scenario into an ESIO-based file.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h Open, writable handle in which details will be saved.
     */
    virtual void save(
            const esio_handle h) const;

    /**
     * Populate scenario from an ESIO-based file.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param h       Open, readable handle from which details will be loaded.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc support::definition_base::options_description() */
    boost::program_options::options_description options_description();


    /**
     * Initialize antioch objects that persist over whole program
     * life.  Must be called after \ref species_names and \ref
     * chem_input_file have been properly initialized.
     */
    void init_antioch();

    /**
     * Given conserved state, compute required thermodynamic and
     * transport quantities for evaluating Navier-Stokes operator.
     *
     * @param[in] e       Total energy per unit volume.
     * @param[in] m       Pointer to momentum components, m[0] = x-momentum, etc.
     * @param[in] rho     Mixture density.
     * @param[in] species Species densities.
     * @param[in] cs      Species mass fractions. 
     * @param[out] T      Temperature.
     * @param[out] p      Pressure.
     * @param[out] Ds     Mass diffusivities. 
     * @param[out] mu     Dynamic viscosity.
     * @param[out] kap    Thermal conductivity.
     * @param[out] hs     Species enthalpies.
     * @param[out] om     Reaction source terms.
     */
    void evaluate (const real_t  e,
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
                   real_t* om) const;

    /**
     * Given temperature and mass fractions, compute internal energy
     *
     * @param[in] T              Temperature
     * @param[in] mass_fractions Mass fractions
     * @param[out] e             Internal energy
     */
    real_t e_from_T (const real_t  T,
                     const std::vector<real_t> mass_fractions) const;

    /**
     * Report the number of species in the mixture.
     */
    const std::size_t Ns() const { return species_names.size(); }

    /**
     * Version of libantioch reported by Antioch::get_antioch_version
     */
    int antioch_ver;

    /**
     * Vector of names of species in mixture
     */
    std::vector<std::string> species_names;

    /**
     * Kinetics input filename
     */
    std::string chem_input_file;

    /**
     * The Lewis number used to compute diffusivities
     */
    real_t Le;

    /**
     * The ratio of bulk viscosity to dynamic viscosity.
     */
    real_t alpha;


    /**
     * Antioch::ChemicalMixture object, required for anything with
     * antioch
     */
    shared_ptr<Antioch::ChemicalMixture<real_t> > mixture;

    /**
     * Antioch::ReactionSet object, used by Kinetics in evaluating
     * reaction srcs
     */
    shared_ptr<Antioch::ReactionSet<real_t> > reactions;

    /**
     * Antioch::CEAThermodynamics object, used to compute inputs for
     * reaction src evaluation
     */
    shared_ptr<Antioch::CEAThermodynamics<real_t> > cea_thermo;

    /**
     * Antioch::StatMechThermodynamics object, used to compute
     * temperature from internal energy
     */
    shared_ptr<Antioch::StatMechThermodynamics<real_t> > sm_thermo;

    /**
     * Antioch::KineticsEvaluator object, used to actually compute
     * mass sources due to chemical reactions.
     *
     * FIXME: This object used in this way probably makes this class
     * thread unsafe.  See comments in antioch/kinetics_evaluator.h
     * and issue #2798.
     */
    shared_ptr<Antioch::KineticsEvaluator<real_t> > kinetics;

    
    // TODO: Add antioch support for transport properties

};

} // namespace reacting

} // namespace suzerain

#endif // HAVE_ANTIOCH

#endif // SUZERAIN_ANTIOCH_CONSTITUTIVE_HPP
