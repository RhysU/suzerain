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

#ifndef SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP
#define SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP

/** @file
 * Classes handling perfect gas problem scenario parameters.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

// FIXME Break bulk_rho and bulk_rho_u into channel_definition class

/** @file
 * Provides classes handling problem scenario parameters which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain {

namespace perfect {

/**
 * Holds nondimensional parameters like the Reynolds and Prandtl numbers as
 * well as nondimensional problem geometry.  See the Suzerain model document's
 * nondimensionalization section for more information.
 */
class scenario_definition : public support::definition_base
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    scenario_definition();

    /**
     * Construct an instance with the given parameter values.
     *
     * @param Re         Reynolds number.
     * @param Ma         Mach number.
     * @param Pr         Prandtl number.
     * @param bulk_rho   Bulk density target.
     * @param bulk_rho_u Bulk streamwise momentum target.
     * @param alpha      Ratio of bulk to dynamic viscosity.
     * @param beta       Temperature power law exponent.
     * @param gamma      Ratio of specific heats.
     */
    scenario_definition(const real_t Re,
                        const real_t Ma,
                        const real_t Pr,
                        const real_t bulk_rho,
                        const real_t bulk_rho_u,
                        const real_t alpha,
                        const real_t beta,
                        const real_t gamma);

    /** @copydoc support::definition_base::options_description() */
    boost::program_options::options_description options_description();

    /**
     * The Reynolds number \f$\mbox{Re}=\frac{\rho_{0} u_{0}
     * l_{0}}{\mu_{0}}\f$.
     */
    real_t Re;

    /**
     * The Mach number \f$\mbox{Ma}=\frac{u_{0}}{a_{0}}\f$.
     */
    real_t Ma;

    /**
     * The Prandtl number \f$\mbox{Pr}=\frac{\mu_{0}
     * C_{p}}{\kappa_{0}}\f$.
     */
    real_t Pr;

    /**
     * The bulk density used as a target for integral constraints.
     */
    real_t bulk_rho;

    /**
     * The bulk streamwise momentum used as a target for integral constraints.
     */
    real_t bulk_rho_u;

    /**
     * The ratio of bulk viscosity to dynamic viscosity according to \f$
     * \mu_{B} = \alpha \mu \f$ or equivalently \f$ \lambda = \left( \alpha -
     * \frac{2}{3}\mu \right)\f$.
     */
    real_t alpha;

    /**
     * The temperature power law exponent \f$\beta\f$ where
     * \f$\frac{\mu}{\mu_0} = \left(\frac{T}{T_0}\right)^{\beta}\f$.
     */
    real_t beta;

    /**
     * The ratio of specific heats \f$\gamma=C_p/C_v\f$.
     */
    real_t gamma;

};

/**
 * Save a scenario_definition in an ESIO-based file.
 *
 * @param h        Open, writable handle in which details will be saved.
 * @param scenario Scenario to be saved.
 */
void save(const esio_handle h,
          const scenario_definition& scenario);

/**
 * Load a scenario_definition from an ESIO-based file.
 *
 * @param h        Open, readable handle from which details will be loaded.
 * @param scenario Scenario to be saved.
 */
void load(const esio_handle h,
          scenario_definition& scenario);

} // namespace perfect

} // namespace suzerain

#endif // SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP
