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

#ifndef SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP
#define SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

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
class scenario_definition
    : public virtual support::definition_base
    , public virtual support::loadable
    , public virtual support::overridable<scenario_definition>
    , public virtual support::populatable<scenario_definition>
    , public virtual support::savable
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
     * @param bulk_rho_E Bulk streamwise total energy target.
     * @param alpha      Ratio of bulk to dynamic viscosity.
     * @param beta       Temperature power law exponent.
     * @param gamma      Ratio of specific heats.
     */
    scenario_definition(const real_t Re,
                        const real_t Ma,
                        const real_t Pr,
                        const real_t bulk_rho,
                        const real_t bulk_rho_u,
                        const real_t bulk_rho_E,
                        const real_t alpha,
                        const real_t beta,
                        const real_t gamma);

    /** Virtual destructor to permit use as a base class */
    virtual ~scenario_definition();

    /** @copydoc support::populatable::populate */
    virtual void populate(
            const scenario_definition& that,
            const bool verbose = false);

    /** @copydoc support::overridable::override */
    virtual void override(
            const scenario_definition& that,
            const bool verbose = false);

    /** @copydoc support::savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc support::loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

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
     * The value \c NaN implies the constraint should not be enforced.
     */
    real_t bulk_rho;

    /**
     * The bulk streamwise momentum used as a target for integral constraints.
     * The value \c NaN implies the constraint should not be enforced.
     */
    real_t bulk_rho_u;

    /**
     * The bulk total energy used as a target for integral constraints.
     * The value \c NaN implies the constraint should not be enforced.
     */
    real_t bulk_rho_E;

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

} // namespace perfect

} // namespace suzerain

#endif // SUZERAIN_PERFECT_SCENARIO_DEFINITION_HPP
