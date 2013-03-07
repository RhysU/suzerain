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

#ifndef SUZERAIN_SINGLE_IDEAL_GAS_CONSTITUTIVE_HPP
#define SUZERAIN_SINGLE_IDEAL_GAS_CONSTITUTIVE_HPP

/** @file
 * Classes handling reacting channel flow problem scenario parameters.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

/** @file 
 * Provides constitutive models for a single species ideal gas
 * with constant specific heats, constant Prandtl number, and power
 * law viscosity.
 */

namespace suzerain {

namespace reacting {

/**
 * Holds parameters defining channel flow case.
 */
class single_ideal_gas_constitutive : public support::definition_base
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    single_ideal_gas_constitutive();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param bulk_rho   Bulk density target.
     * @param bulk_rho_u Bulk streamwise momentum target.
     */
    single_ideal_gas_constitutive(const real_t Cp,
				 const real_t Cv,
				 const real_t Pr,
				 const real_t T0,
				 const real_t mu0,
				 const real_t beta,
				 const real_t alpha);


    /** Virtual destructor to permit use as a base class */
    virtual ~single_ideal_gas_constitutive();

    /**
     * Populate any NaN members in \c this with values from \c that.
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const single_ideal_gas_constitutive& that,
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
            const single_ideal_gas_constitutive& that,
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
     * @param[in] species Species densities (should not contain anything here).
     * @param[in] cs      Species mass fractions (should not contain anything here).
     * @param[out] T      Temperature.
     * @param[out] p      Pressure.
     * @param[out] Ds     Mass diffusivities (won't have anything here).
     * @param[out] mu     Dynamic viscosity.
     * @param[out] kap    Thermal conductivity.
     * @param[out] hs     Species enthalpies (won't have anything here).
     * @param[out] om     Reaction source terms (won't have anything here).
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
                   real_t* om);

    /**
     * Report the number of species (just one here).
     */
    const std::size_t Ns() const { return 1; }

    /**
     * The specific heat at constant pressure.
     */
    real_t Cp;

    /**
     * The specific heat at constant volume.
     */
    real_t Cv;

    /**
     * The Prandtl number.
     */
    real_t Pr;

    /**
     * The reference temperature for the viscosity power law.
     */
    real_t T0;

    /**
     * The reference viscosity for the viscosity power law.
     */
    real_t mu0;

    /**
     * The power in the viscosity power law.
     */
    real_t beta;

    /**
     * The ratio of bulk viscosity to dynamic viscosity.
     */
    real_t alpha;

};

} // namespace reacting

} // namespace suzerain

#endif // SUZERAIN_SINGLE_IDEAL_GAS_CONSTITUTIVE_HPP
