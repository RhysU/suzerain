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
        initialize_options();
    }

    /**
     * Construct an instance with the given parameter values.
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
    ScenarioDefinition(FPT Re,
                       FPT Ma,
                       FPT Pr,
                       FPT bulk_rho,
                       FPT bulk_rhou,
                       FPT alpha,
                       FPT beta,
                       FPT gamma,
                       FPT Lx,
                       FPT Ly,
                       FPT Lz)
        : IDefinition("Nondimensional scenario parameters"),
          Re(Re),
          Ma(Ma),
          Pr(Pr),
          bulk_rho(bulk_rho),
          bulk_rhou(bulk_rhou),
          alpha(alpha),
          beta(beta),
          gamma(gamma),
          Lx(Lx),
          Ly(Ly),
          Lz(Lz)
    {
        initialize_options();
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
    void initialize_options();
};

template< typename FPT >
void ScenarioDefinition<FPT>::initialize_options()
{

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ensure_positive(::suzerain::validation::ensure_positive<FPT>);
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ensure_nonnegative(::suzerain::validation::ensure_nonnegative<FPT>);

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    using ::boost::math::isnan;
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;
    using ::std::auto_ptr;
    using ::std::bind2nd;

    { // Re
        auto_ptr<typed_value<FPT> > v(value(&this->Re));
        v->notifier(bind2nd(ensure_positive, "Re"));
        if (!(isnan)(Re)) v->default_value(Re);
        this->add_options()("Re", v.release(), "Reynolds number");
    }

    { // Ma
        auto_ptr<typed_value<FPT> > v(value(&this->Ma));
        v->notifier(bind2nd(ensure_positive, "Ma"));
        if (!(isnan)(Ma)) v->default_value(Ma);
        this->add_options()("Ma", v.release(), "Mach number");
    }

    { // Pr
        auto_ptr<typed_value<FPT> > v(value(&this->Pr));
        v->notifier(bind2nd(ensure_positive, "Pr"));
        if (!(isnan)(Pr)) v->default_value(Pr);
        this->add_options()("Pr", v.release(), "Prandtl number");
    }

    { // bulk_rho
        auto_ptr<typed_value<FPT> > v(value(&this->bulk_rho));
        v->notifier(bind2nd(ensure_nonnegative, "bulk_rho"));
        if (!(isnan)(bulk_rho)) v->default_value(bulk_rho);
        this->add_options()("bulk_rho", v.release(), "bulk density target");
    }

    { // bulk_rhou
        auto_ptr<typed_value<FPT> > v(value(&this->bulk_rhou));
        v->notifier(bind2nd(ensure_nonnegative, "bulk_rhou"));
        if (!(isnan)(bulk_rhou)) v->default_value(bulk_rhou);
        this->add_options()("bulk_rhou", v.release(), "bulk momentum target");
    }

    { // alpha
        auto_ptr<typed_value<FPT> > v(value(&this->alpha));
        v->notifier(bind2nd(ensure_nonnegative, "alpha"));
        if (!(isnan)(alpha)) v->default_value(alpha);
        this->add_options()("alpha", v.release(),
                "Ratio of bulk to dynamic viscosity");
    }

    { // beta
        auto_ptr<typed_value<FPT> > v(value(&this->beta));
        v->notifier(bind2nd(ensure_nonnegative, "beta"));
        if (!(isnan)(beta)) v->default_value(beta);
        this->add_options()("beta", v.release(),
                "Temperature power law exponent");
    }

    { // gamma
        auto_ptr<typed_value<FPT> > v(value(&this->gamma));
        v->notifier(bind2nd(ensure_positive, "gamma"));
        if (!(isnan)(gamma)) v->default_value(gamma);
        this->add_options()("gamma", v.release(), "Ratio of specific heats");
    }

    { // Lx
        auto_ptr<typed_value<FPT> > v(value(&this->Lx));
        v->notifier(bind2nd(ensure_positive, "Lx"));
        if (!(isnan)(Lx)) v->default_value(Lx);
        this->add_options()("Lx", v.release(),
                "Nondimensional grid length in streamwise X direction");
    }

    { // Ly
        auto_ptr<typed_value<FPT> > v(value(&this->Ly));
        v->notifier(bind2nd(ensure_positive, "Ly"));
        if (!(isnan)(Ly)) v->default_value(Ly);
        this->add_options()("Ly", v.release(),
                "Nondimensional grid length in wall normal Y direction");
    }

    { // Lz
        auto_ptr<typed_value<FPT> > v(value(&this->Lz));
        v->notifier(bind2nd(ensure_positive, "Lz"));
        if (!(isnan)(Lz)) v->default_value(Lz);
        this->add_options()("Lz", v.release(),
                "Nondimensional grid length in spanwise Z direction");
    }
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_SCENARIO_DEFINITION_HPP
