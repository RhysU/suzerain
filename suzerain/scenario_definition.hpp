/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * scenario_definition.hpp: classes handling problem scenario parameters
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_SCENARIO_DEFINITION_HPP
#define __SUZERAIN_SCENARIO_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem scenario parameters which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain {

namespace problem {

/**
 * Holds nondimensional parameters like the Reynolds and Prandtl numbers as
 * well as nondimensional problem geometry.  See the Suzerain model document's
 * nondimensionalization section for more information.
 */
template< typename FPT = double >
class ScenarioDefinition : public IDefinition
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    ScenarioDefinition()
        : IDefinition("Nondimensional scenario parameters"),
          Re(std::numeric_limits<FPT>::quiet_NaN()),
          Ma(std::numeric_limits<FPT>::quiet_NaN()),
          Pr(std::numeric_limits<FPT>::quiet_NaN()),
          bulk_rho(std::numeric_limits<FPT>::quiet_NaN()),
          bulk_rhou(std::numeric_limits<FPT>::quiet_NaN()),
          alpha(std::numeric_limits<FPT>::quiet_NaN()),
          beta(std::numeric_limits<FPT>::quiet_NaN()),
          gamma(std::numeric_limits<FPT>::quiet_NaN()),
          Lx(std::numeric_limits<FPT>::quiet_NaN()),
          Ly(std::numeric_limits<FPT>::quiet_NaN()),
          Lz(std::numeric_limits<FPT>::quiet_NaN())
    {
        initialize_options(NULL, NULL, NULL,
                           NULL, NULL,
                           NULL, NULL, NULL,
                           NULL, NULL, NULL);
    }

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
     * @param Lx        Domain length in the X direction.
     * @param Ly        Domain length in the Y direction.
     * @param Lz        Domain length in the Z direction.
     */
    ScenarioDefinition(const char * Re,
                       const char * Ma,
                       const char * Pr,
                       const char * bulk_rho,
                       const char * bulk_rhou,
                       const char * alpha,
                       const char * beta,
                       const char * gamma,
                       const char * Lx,
                       const char * Ly,
                       const char * Lz)
        : IDefinition("Nondimensional scenario parameters"),
          Re(std::numeric_limits<FPT>::quiet_NaN()),
          Ma(std::numeric_limits<FPT>::quiet_NaN()),
          Pr(std::numeric_limits<FPT>::quiet_NaN()),
          bulk_rho(std::numeric_limits<FPT>::quiet_NaN()),
          bulk_rhou(std::numeric_limits<FPT>::quiet_NaN()),
          alpha(std::numeric_limits<FPT>::quiet_NaN()),
          beta(std::numeric_limits<FPT>::quiet_NaN()),
          gamma(std::numeric_limits<FPT>::quiet_NaN()),
          Lx(std::numeric_limits<FPT>::quiet_NaN()),
          Ly(std::numeric_limits<FPT>::quiet_NaN()),
          Lz(std::numeric_limits<FPT>::quiet_NaN())
    {
        initialize_options(Re, Ma, Pr,
                           bulk_rho, bulk_rhou,
                           alpha, beta, gamma,
                           Lx, Ly, Lz);
    }

    /**
     * The Reynolds number \f$\mbox{Re}=\frac{\rho_{0} u_{0}
     * l_{0}}{\mu_{0}}\f$.
     */
    FPT Re;

    /**
     * The Mach number \f$\mbox{Ma}=\frac{u_{0}}{a_{0}}\f$.
     */
    FPT Ma;

    /**
     * The Prandtl number \f$\mbox{Pr}=\frac{\mu_{0}
     * C_{p}}{\kappa_{0}}\f$.
     */
    FPT Pr;

    /**
     * The bulk density used as a target for integral constraints.
     */
    FPT bulk_rho;

    /**
     * The bulk streamwise momentum used as a target for integral constraints.
     */
    FPT bulk_rhou;

    /**
     * The ratio of bulk viscosity to dynamic viscosity according to \f$
     * \mu_{B} = \alpha \mu \f$ or equivalently \f$ \lambda = \left( \alpha -
     * \frac{2}{3}\mu \right)\f$.
     */
    FPT alpha;

    /**
     * The temperature power law exponent \f$\beta\f$ where
     * \f$\frac{\mu}{\mu_0} = \left(\frac{T}{T_0}\right)^{\beta}\f$.
     */
    FPT beta;

    /**
     * The ratio of specific heats \f$\gamma=C_p/C_v\f$.
     */
    FPT gamma;

    /**
     * The domain length in the X direction.
     */
    FPT Lx;

    /**
     * The domain length in the Y direction.
     */
    FPT Ly;

    /**
     * The domain length in the Z direction.
     */
    FPT Lz;

private:
    /** Options initialization common to all constructors */
    void initialize_options(const char * default_Re,
                            const char * default_Ma,
                            const char * default_Pr,
                            const char * default_bulk_rho,
                            const char * default_bulk_rhou,
                            const char * default_alpha,
                            const char * default_beta,
                            const char * default_gamma,
                            const char * default_Lx,
                            const char * default_Ly,
                            const char * default_Lz);

    static void parse_positive(const std::string& s, FPT *t, const char *n) {
        const FPT v = suzerain::exprparse<FPT>(s, n);
        suzerain::validation::ensure_positive(v, n);
        *t = v;
    }

    static void parse_nonnegative(const std::string& s, FPT *t, const char *n) {
        const FPT v = suzerain::exprparse<FPT>(s, n);
        suzerain::validation::ensure_nonnegative(v, n);
        *t = v;
    }
};

template< typename FPT >
void ScenarioDefinition<FPT>::initialize_options(
        const char * default_Re,
        const char * default_Ma,
        const char * default_Pr,
        const char * default_bulk_rho,
        const char * default_bulk_rhou,
        const char * default_alpha,
        const char * default_beta,
        const char * default_gamma,
        const char * default_Lx,
        const char * default_Ly,
        const char * default_Lz)
{
    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;

    // Re
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Re, "Re"));
    if (default_Re) p->default_value(default_Re);
    this->add_options()("Re", p.release(), "Reynolds number");

    // Ma
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Ma, "Ma"));
    if (default_Ma) p->default_value(default_Ma);
    this->add_options()("Ma", p.release(), "Mach number");

    // Pr
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Pr, "Pr"));
    if (default_Pr) p->default_value(default_Pr);
    this->add_options()("Pr", p.release(), "Prandtl number");

    // bulk_rho
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &bulk_rho, "bulk_rho"));
    if (default_bulk_rho) p->default_value(default_bulk_rho);
    this->add_options()("bulk_rho", p.release(), "Bulk density target");

    // bulk_rhou
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &bulk_rhou, "bulk_rhou"));
    if (default_bulk_rhou) p->default_value(default_bulk_rhou);
    this->add_options()("bulk_rhou", p.release(), "Bulk momentum target");

    // alpha
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &alpha, "alpha"));
    if (default_alpha) p->default_value(default_alpha);
    this->add_options()("alpha", p.release(),
                        "Ratio of bulk to dynamic viscosity");

    // beta
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &beta, "beta"));
    if (default_beta) p->default_value(default_beta);
    this->add_options()("beta", p.release(), "Temperature power law exponent");

    // gamma
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &gamma, "gamma"));
    if (default_gamma) p->default_value(default_gamma);
    this->add_options()("gamma", p.release(), "Ratio of specific heats");

    // Lx
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Lx, "Lx"));
    if (default_Lx) p->default_value(default_Lx);
    this->add_options()("Lx", p.release(),
            "Nondimensional grid length in streamwise X direction");

    // Ly
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Ly, "Ly"));
    if (default_Ly) p->default_value(default_Ly);
    this->add_options()("Ly", p.release(),
            "Nondimensional grid length in wall normal Y direction");

    // Lz
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Lz, "Lz"));
    if (default_Lz) p->default_value(default_Lz);
    this->add_options()("Lz", p.release(),
            "Nondimensional grid length in spanwise Z direction");
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_SCENARIO_DEFINITION_HPP
