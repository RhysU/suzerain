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
     * Construct an instance with the given parameters.
     *
     * @param default_Re        Default Reynolds number.
     * @param default_Pr        Default Prandtl number.
     * @param default_Ma        Default Mach number.
     * @param default_bulk_rho  Default bulk density target.
     * @param default_bulk_rhou Default bulk streamwise momentum target.
     * @param default_gamma     Default ratio of specific heats.
     * @param default_beta      Default temperature power law exponent.
     * @param default_Lx        Default domain length in the X direction.
     * @param default_Ly        Default domain length in the Y direction.
     * @param default_Lz        Default domain length in the Z direction.
     */
    explicit ScenarioDefinition(FPT default_Re        = 0,
                                FPT default_Pr        = 0,
                                FPT default_Ma        = 0,
                                FPT default_bulk_rho  = 0,
                                FPT default_bulk_rhou = 0,
                                FPT default_gamma     = 0,
                                FPT default_beta      = 0,
                                FPT default_Lx        = 0,
                                FPT default_Ly        = 0,
                                FPT default_Lz        = 0);

    /**
     * The Reynolds number \f$\mbox{Re}=\frac{\rho_{0} u_{0}
     * l_{0}}{\mu_{0}}\f$.
     */
    FPT Re;

    /**
     * The Prandtl number \f$\mbox{Pr}=\frac{\mu_{0}
     * C_{p}}{\kappa_{0}}\f$.
     */
    FPT Pr;

    /**
     * The Mach number \f$\mbox{Ma}=\frac{u_{0}}{a_{0}}\f$.
     */
    FPT Ma;

    /**
     * The bulk density used as a target for integral constraints.
     */
    FPT bulk_rho;

    /**
     * The bulk streamwise momentum used as a target for integral constraints.
     */
    FPT bulk_rhou;

    /**
     * The ratio of specific heats \f$\gamma=C_p/C_v\f$.
     */
    FPT gamma;

    /**
     * The temperature power law exponent \f$\beta\f$ where
     * \f$\frac{\mu}{\mu_0} = \left(\frac{T}{T_0}\right)^{\beta}\f$.
     */
    FPT beta;

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
};

template< typename FPT >
ScenarioDefinition<FPT>::ScenarioDefinition(
        FPT default_Re,
        FPT default_Pr,
        FPT default_Ma,
        FPT default_bulk_rho,
        FPT default_bulk_rhou,
        FPT default_gamma,
        FPT default_beta,
        FPT default_Lx,
        FPT default_Ly,
        FPT default_Lz)
    : IDefinition("Nondimensional scenario parameters"),
      Re(default_Re),
      Pr(default_Pr),
      Ma(default_Ma),
      bulk_rho(default_bulk_rho),
      bulk_rhou(default_bulk_rhou),
      gamma(default_gamma),
      beta(default_beta),
      Lx(default_Lx),
      Ly(default_Ly),
      Lz(default_Lz)
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
        auto_ptr<typed_value<FPT> > v(value(&this->Re));
        if (default_Re) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Re"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Re"));
        }
        v->default_value(default_Re);
        this->add_options()("Re", v.release(), "Reynolds number");
    }

    { // Pr
        auto_ptr<typed_value<FPT> > v(value(&this->Pr));
        if (default_Pr) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Pr"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Pr"));
        }
        v->default_value(default_Pr);
        this->add_options()("Pr", v.release(), "Prandtl number");
    }

    { // Ma
        auto_ptr<typed_value<FPT> > v(value(&this->Ma));
        if (default_Ma) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "Ma"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "Ma"));
        }
        v->default_value(default_Ma);
        this->add_options()("Ma", v.release(), "Mach number");
    }

    { // bulk_rho
        auto_ptr<typed_value<FPT> > v(value(&this->bulk_rho));
        if (default_bulk_rho) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "bulk_rho"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "bulk_rho"));
        }
        v->default_value(default_bulk_rho);
        this->add_options()("bulk_rho", v.release(),
                "bulk density target");
    }

    { // bulk_rhou
        auto_ptr<typed_value<FPT> > v(value(&this->bulk_rhou));
        if (default_bulk_rhou) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "bulk_rhou"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "bulk_rhou"));
        }
        v->default_value(default_bulk_rhou);
        this->add_options()("bulk_rhou", v.release(),
                "bulk streamwise momentum target");
    }

    { // gamma
        auto_ptr<typed_value<FPT> > v(value(&this->gamma));
        if (default_gamma) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,    "gamma"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "gamma"));
        }
        v->default_value(default_gamma);
        this->add_options()("gamma", v.release(), "Ratio of specific heats");
    }

    { // beta
        auto_ptr<typed_value<FPT> > v(value(&this->beta));
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
        auto_ptr<typed_value<FPT> > v(value(&this->Lx));
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
        auto_ptr<typed_value<FPT> > v(value(&this->Ly));
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
        auto_ptr<typed_value<FPT> > v(value(&this->Lz));
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
