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
#include <suzerain/problem.hpp>
#include <suzerain/types.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem scenario parameters, which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain {

namespace problem {

/**
 * Holds nondimensional parameters like the Reynolds and Prandtl numbers.  See
 * the Suzerain model document's nondimensionalization section for more
 * information.
 *
 */
template< typename FPT = double >
class ScenarioDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default parameters.
     *
     * @param default_Re Default Reynolds number.
     * @param default_Pr Default Prandtl number.
     * @param default_gamma Default ratio of specific heats.
     * @param default_beta Default temperature power law exponent.
     */
    ScenarioDefinition(FPT default_Re,
                       FPT default_Pr    = 0.7,
                       FPT default_gamma = 1.4,
                       FPT default_beta  = 0.666 );

    /**
     * Retrieve the Reynolds number $\mbox{Re}=\frac{\rho_{0} u_{0}
     * l_{0}}{\mu_{0}}\f$.
     *
     * @return the Reynolds number.
     */
    FPT Re() const { return Re_; }

    /**
     * Retrieve the Prandtl number \f$\mbox{Pr}=\frac{\mu_{0}
     * C_{p}}{\kappa_{0}}\f$.
     *
     * @return the Prandtl number.
     */
    FPT Pr() const { return Pr_; }

    /**
     * Retrieve the ratio of specific heats \f$\gamma=C_p/C_v\f$.
     *
     * @return the ratio of specific heats.
     */
    FPT gamma() const { return gamma_; }

    /**
     * Retrieve the temperature power law exponent \f$\beta\f$ where
     * \f$\frac{\mu}{\mu_0} = \left(\frac{T}{T_0}\right)^{\beta}\f$.
     *
     * @return the temperature power law exponent.
     */
    FPT beta() const { return beta_; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

    FPT Re_;     /**< Stores the Reynolds number */
    FPT Pr_;     /**< Stores the Prandtl number */
    FPT gamma_;  /**< Stores the ratio of specific heats */
    FPT beta_;   /**< Stores the temperature power law exponent */
};

template< typename FPT >
ScenarioDefinition<FPT>::ScenarioDefinition(
        FPT default_Re,
        FPT default_Pr,
        FPT default_gamma,
        FPT default_beta)
    : options_("Nondimensional scenario parameters"),
      Re_(default_Re),
      Pr_(default_Pr),
      gamma_(default_gamma),
      beta_(default_beta)
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_positive;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_positive_FPT(ensure_positive<FPT>);

    options_.add_options()
        ("Re", po::value<FPT>(&Re_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"Re"))
            ->default_value(default_length),
        "Reynolds number")
        ("Pr", po::value<FPT>(&Pr_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"Pr"))
            ->default_value(default_length),
        "Prandtl number")
        ("gamma", po::value<FPT>(&gamma_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"gamma"))
            ->default_value(default_length),
        "Ratio of specific heats")
        ("beta", po::value<FPT>(&beta_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"beta"))
            ->default_value(default_length),
        "Temperature power law exponent")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_SCENARIO_DEFINITION_HPP
