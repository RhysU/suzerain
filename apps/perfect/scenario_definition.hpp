//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// scenario_definition.hpp: classes handling problem scenario parameters
// $Id$

#ifndef SCENARIO_DEFINITION_HPP
#define SCENARIO_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling problem scenario parameters which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain { namespace perfect {

/**
 * Holds nondimensional parameters like the Reynolds and Prandtl numbers as
 * well as nondimensional problem geometry.  See the Suzerain model document's
 * nondimensionalization section for more information.
 */
class ScenarioDefinition : public suzerain::problem::IDefinition
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    ScenarioDefinition();

    /**
     * Construct an instance with the given parameter values.
     * Parameter values are evaluated via suzerain::exprparse().
     *
     * @param Re        Reynolds number.
     * @param Ma        Mach number.
     * @param Pr        Prandtl number.
     * @param bulk_rho  Bulk density target.
     * @param bulk_rhou Bulk streamwise momentum target.
     * @param alpha     Ratio of bulk to dynamic viscosity.
     * @param beta      Temperature power law exponent.
     * @param gamma     Ratio of specific heats.
     */
    ScenarioDefinition(const char * Re,
                       const char * Ma,
                       const char * Pr,
                       const char * bulk_rho,
                       const char * bulk_rhou,
                       const char * alpha,
                       const char * beta,
                       const char * gamma);

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
    real_t bulk_rhou;

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

private:
    /** Options initialization common to all constructors */
    void initialize_options(const char * default_Re,
                            const char * default_Ma,
                            const char * default_Pr,
                            const char * default_bulk_rho,
                            const char * default_bulk_rhou,
                            const char * default_alpha,
                            const char * default_beta,
                            const char * default_gamma);
};

} /* namespace perfect */ } /* namespace suzerain */

#endif // SCENARIO_DEFINITION_HPP
