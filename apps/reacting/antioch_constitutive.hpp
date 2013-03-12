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


namespace suzerain {

namespace reacting {

/**
 * Provides constitutive models for a mixture of ideal gasses
 * leveraging libantioch
 */
class antioch_constitutive : public support::definition_base
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
     * Report the number of species in the mixture.
     */
    const std::size_t Ns() const { return species_names.size(); }

protected:

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
};

} // namespace reacting

} // namespace suzerain

#endif // HAVE_ANTIOCH

#endif // SUZERAIN_ANTIOCH_CONSTITUTIVE_HPP
