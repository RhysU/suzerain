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
     * Construct an instance with the given default parameters.
     *
     * @param default_Re Default Reynolds number.
     * @param default_Pr Default Prandtl number.
     * @param default_gamma Default ratio of specific heats.
     * @param default_beta Default temperature power law exponent.
     * @param default_Lx Default domain length in the X direction.
     * @param default_Ly Default domain length in the Y direction.
     * @param default_Lz Default domain length in the Z direction.
     */
    ScenarioDefinition(FPT default_Re,
                       FPT default_Pr,
                       FPT default_gamma,
                       FPT default_beta,
                       FPT default_Lx,
                       FPT default_Ly,
                       FPT default_Lz);

    /**
     * Retrieve the Reynolds number \f$\mbox{Re}=\frac{\rho_{0} u_{0}
     * l_{0}}{\mu_{0}}\f$.
     *
     * @return the Reynolds number.
     */
    FPT Re() const { return Re_; }

    /** \copydoc Re() const */
    FPT& Re() { return Re_; }

    /**
     * Retrieve the Prandtl number \f$\mbox{Pr}=\frac{\mu_{0}
     * C_{p}}{\kappa_{0}}\f$.
     *
     * @return the Prandtl number.
     */
    FPT Pr() const { return Pr_; }

    /** \copydoc Pr() const */
    FPT& Pr() { return Pr_; }

    /**
     * Retrieve the ratio of specific heats \f$\gamma=C_p/C_v\f$.
     *
     * @return the ratio of specific heats.
     */
    FPT gamma() const { return gamma_; }

    /** \copydoc gamma() const */
    FPT& gamma() { return gamma_; }

    /**
     * Retrieve the temperature power law exponent \f$\beta\f$ where
     * \f$\frac{\mu}{\mu_0} = \left(\frac{T}{T_0}\right)^{\beta}\f$.
     *
     * @return the temperature power law exponent.
     */
    FPT beta() const { return beta_; }

    /** \copydoc beta() const */
    FPT& beta() { return beta_; }

    /**
     * Retrieve the domain length in the X direction.
     *
     * @return the domain's X length.
     */
    FPT Lx() const { return Lx_; }

    /** \copydoc Lx() const */
    FPT& Lx() { return Lx_; }

    /**
     * Retrieve the domain length in the Y direction.
     *
     * @return the domain's Y length.
     */
    FPT Ly() const { return Ly_; }

    /** \copydoc Ly() const */
    FPT& Ly() { return Ly_; }

    /**
     * Retrieve the domain length in the Z direction.
     *
     * @return the domain's Z length.
     */
    FPT Lz() const { return Lz_; }

    /** \copydoc Lz() const */
    FPT& Lz() { return Lz_; }

private:

    FPT Re_;     /**< Stores the Reynolds number */
    FPT Pr_;     /**< Stores the Prandtl number */
    FPT gamma_;  /**< Stores the ratio of specific heats */
    FPT beta_;   /**< Stores the temperature power law exponent */
    FPT Lx_;     /**< Stores the X direction length */
    FPT Ly_;     /**< Stores the Y direction length */
    FPT Lz_;     /**< Stores the Z direction length */

};

template< typename FPT >
ScenarioDefinition<FPT>::ScenarioDefinition(
        FPT default_Re,
        FPT default_Pr,
        FPT default_gamma,
        FPT default_beta,
        FPT default_Lx,
        FPT default_Ly,
        FPT default_Lz)
    : IDefinition("Nondimensional scenario parameters")
{
    using ::std::auto_ptr;
    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_positive_FPT(ensure_positive<FPT>);
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_nonnegative_FPT(ensure_nonnegative<FPT>);

    // Complicated add_options() calls done to allow changing the validation
    // routine in use when the default provided value is zero.  Zero is
    // generally used a NOP value by some client code.

    { // Re
        auto_ptr<typed_value<FPT> > v(value(&Re_));
        if (default_Re) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Re"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Re"));
        }
        v->default_value(default_Re);
        this->add_options()("Re", v.release(), "Reynolds number");
    }

    { // Pr
        auto_ptr<typed_value<FPT> > v(value(&Pr_));
        if (default_Pr) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Pr"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Pr"));
        }
        v->default_value(default_Pr);
        this->add_options()("Pr", v.release(), "Prandtl number");
    }

    { // gamma
        auto_ptr<typed_value<FPT> > v(value(&gamma_));
        if (default_gamma) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "gamma"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "gamma"));
        }
        v->default_value(default_gamma);
        this->add_options()("gamma", v.release(), "Ratio of specific heats");
    }

    { // beta
        auto_ptr<typed_value<FPT> > v(value(&beta_));
        if (default_beta) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "beta"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "beta"));
        }
        v->default_value(default_beta);
        this->add_options()("beta", v.release(),
                "Temperature power law exponent");
    }

    { // Lx
        auto_ptr<typed_value<FPT> > v(value(&Lx_));
        if (default_Lx) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Lx"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Lx"));
        }
        v->default_value(default_Lx);
        this->add_options()("Lx", v.release(),
                "Nondimensional grid length in streamwise X direction");
    }

    { // Ly
        auto_ptr<typed_value<FPT> > v(value(&Ly_));
        if (default_Ly) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Ly"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Ly"));
        }
        v->default_value(default_Ly);
        this->add_options()("Ly", v.release(),
                "Nondimensional grid length in wall normal Y direction");
    }

    { // Lz
        auto_ptr<typed_value<FPT> > v(value(&Lz_));
        if (default_Lz) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Lz"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Lz"));
        }
        v->default_value(default_Lz);
        this->add_options()("Lz", v.release(),
                "Nondimensional grid length in spanwise Z direction");
    }
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_SCENARIO_DEFINITION_HPP
