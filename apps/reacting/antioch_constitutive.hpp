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

#ifndef SUZERAIN_REACTING_ANTIOCH_CONSTITUTIVE_HPP
#define SUZERAIN_REACTING_ANTIOCH_CONSTITUTIVE_HPP

/** @file
 * Class wrapping libantioch functionality for chemically reacting
 * flow
 */

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

#ifdef SUZERAIN_HAVE_ANTIOCH

#include <antioch/antioch_version.h>

// As of Antioch 0.0.7-ish the compiler-agnostic installation #defines details
// from its test environment.  If the compiler used to build this application
// does not enable C++11 support but the Antioch installation did, we can run
// afoul of decltype not being supported.  What follows is a grotesque
// workaround to permit a C++11 Antioch but a non-C++11 application build...
#if __cplusplus <= 199711L
# warning "Undef-ing ANTIOCH_HAVE_CXX11 as compiler does not claim C++11 conformance!"
# include <antioch_config.h>
# ifdef ANTIOCH_HAVE_CXX11
#  undef ANTIOCH_HAVE_CXX11
# endif
#endif

SUZERAIN_GCC_DIAG_OFF(sign-compare);

#include <antioch/vector_utils_decl.h>
#include <antioch/eigen_utils_decl.h>
#include <antioch/vector_utils.h>
#include <antioch/eigen_utils.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/reaction_set.h>
#include <antioch/cea_thermo.h>
#include <antioch/stat_mech_thermo.h>
#include <antioch/kinetics_evaluator.h>
#include <antioch/mixture_viscosity.h>
#include <antioch/blottner_viscosity.h>
#include <antioch/eucken_thermal_conductivity.h>
#include <antioch/wilke_mixture.h>
#include <antioch/wilke_evaluator.h>

SUZERAIN_GCC_DIAG_ON(sign-compare);

namespace suzerain {

namespace reacting {

/**
 * Provides constitutive models for a mixture of ideal gasses
 * leveraging libantioch
 */
class antioch_constitutive
    : public virtual support::definition_base
    , public virtual support::loadable
    , public virtual support::overridable<antioch_constitutive>
    , public virtual support::populatable<antioch_constitutive>
    , public virtual support::savable
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

    /** @copydoc support::populatable::populate */
    virtual void populate(
            const antioch_constitutive& that,
            const bool verbose = false);

    /** @copydoc support::override::override */
    virtual void override(
            const antioch_constitutive& that,
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
     * Initialize antioch objects that persist over whole program
     * life.  Must be called after \ref species_names and \ref
     * chem_input_file have been properly initialized.
     */
    void init_antioch();

    /**
     * Given conserved state, compute required thermodynamic and
     * transport quantities for evaluating Navier-Stokes operator.
     *
     * @param[in] e       Total energy per unit volume.
     * @param[in] m       Pointer to momentum components, m[0] = x-momentum, etc.
     * @param[in] rho     Mixture density.
     * @param[in] species Species densities.
     * @param[in] cs      Species mass fractions.
     * @param[in] Tinit  Initial guess for temperature.
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
                   const real_t  Tinit,
                   real_t& T,
                   real_t& p,
                   real_t* Ds,
                   real_t& mu,
                   real_t& kap,
                   real_t* hs,
                   real_t* om,
                   real_t& a,
                   real_t& Cv,
                   real_t& Cp) const;

    /**
     * Given conserved state, compute required thermodynamic and
     * transport quantities for evaluating Navier-Stokes operator.
     *
     * @param[in] e       Total energy per unit volume.
     * @param[in] m       Pointer to momentum components, m[0] = x-momentum, etc.
     * @param[in] rho     Mixture density.
     * @param[in] species Species densities.
     * @param[in] cs      Species mass fractions.
     * @param[in] Tinit  Initial guess for temperature.
     * @param[out] T      Temperature.
     * @param[out] p      Pressure.
     * @param[out] Ds     Mass diffusivities.
     * @param[out] mu     Dynamic viscosity.
     * @param[out] kap    Thermal conductivity.
     * @param[out] hs     Species enthalpies.
     * @param[out] om     Reaction source terms.
     */
    void evaluate (const real_t    e,
                   const Vector3r& m,
                   const real_t    rho,
                   const VectorXr& species,
                   const VectorXr& cs,
                   const real_t    Tinit,
                   real_t&   T,
                   real_t&   p,
                   VectorXr& Ds,
                   real_t&   mu,
                   real_t&   kap,
                   VectorXr& hs,
                   VectorXr& om,
                   real_t&   a,
                   real_t&   Cv,
                   real_t&   Cp) const;

    /**
     * Given conserved state, compute required thermodynamic and
     * transport quantities for evaluating Navier-Stokes operator.
     *
     * @param[in] e       Total energy per unit volume.
     * @param[in] m       Pointer to momentum components, m[0] = x-momentum, etc.
     * @param[in] rho     Mixture density.
     * @param[in] species Species densities.
     * @param[in] cs      Species mass fractions.
     * @param[in] Tinit  Initial guess for temperature.
     * @param[out] T      Temperature.
     * @param[out] p      Pressure.
     */
    void evaluate_pressure (const real_t    e,
                            const Vector3r& m,
                            const real_t    rho,
                            const VectorXr& species,
                            const VectorXr& cs,
                            const real_t    Tinit,
                            real_t&   T,
                            real_t&   p) const;

    /**
     * Given conserved state, compute derivatives of pressure wrt
     * state required for linear operator reference quantities.
     *
     * @param[in] e       Total energy per unit volume.
     * @param[in] m       Pointer to momentum components, m[0] = x-momentum, etc.
     * @param[in] rho     Mixture density.
     * @param[in] species Species densities.
     * @param[in] cs      Species mass fractions.
     * @param[in] Tinit  Initial guess for temperature.
     * @param[out] T      Temperature.
     * @param[out] p      Pressure.
     * @param[out] p_rho  Derivative of pressure wrt mixture density
     * @param[out] p_rsum Sum of mass frac(s) * dp/d(rho(s)) for s=2:Ns
     * @param[out] p_m    Derivatives of p wrt m
     * @param[out] p_e    Derivative of p wrt e (= rho*total energy)
     * @param[out] mu     Mixture viscosity
     * @param[out] kap    Mixture thermal conductivity
     * @param[out] Ds     Diffusivities
     * @param[out] gamma  Mixture heat capacity ratio
     * @param[out] a      Mixture speed of sound
     */
    void evaluate_pressure_derivs_and_trans (const real_t    e,
                                             const Vector3r& m,
                                             const real_t    rho,
                                             const VectorXr& species,
                                             const VectorXr& cs,
                                             const real_t    Tinit,
                                             real_t&   T,
                                             real_t&   p,
                                             real_t&   p_rho,
                                             real_t&   p_rsum,
                                             Vector3r& p_m,
                                             real_t&   p_e,
                                             real_t&   mu,
                                             real_t&   kaporCv,
                                             VectorXr& Ds,
                                             real_t&   gamma,
                                             real_t&   a) const;

    /**
     * Given temperature and mass fractions, compute internal energy
     *
     * @param[in] T              Temperature
     * @param[in] mass_fractions Mass fractions
     * @param[out] e             Internal energy
     */
    real_t e_from_T (const real_t  T,
                     const std::vector<real_t> mass_fractions) const;

    /**
     * Given temperature compute specific total energy for each
     * species.
     *
     * @param[in] T       Temparature.
     * @param[out] etots  Specific total energies
     */
    void etots_from_T (const real_t    T,
                       VectorXr&       etots) const;

    /**
     * Given temperature compute specific total enthalpy for each
     * species.
     *
     * @param[in] T       Temparature.
     * @param[out] etots  Specific total enthalpies
     */
    void htots_from_T (const real_t    T,
                       VectorXr&       htots) const;

    /**
     * Given temperature and species mass fractions compute
     * required thermodynamic quantities for the nonreflecting
     * boundary
     *
     * @param[in] T       Temperature
     * @param[in] cs      Species mass fractions
     * @param[out] a      Mixture speed of sound
     * @param[out] gamma  Mixture heat capacity ratio
     * @param[out] R_mix  Mixture gas constant
     * @param[out] etots  Specific total energies
     */
    void evaluate_for_nonreflecting (const real_t      T,
                                     VectorXr&        cs,
                                     real_t&           a,
                                     real_t&       gamma,
                                     real_t&       R_mix,
                                     VectorXr&     etots) const;

    /**
     * Report the number of species in the mixture.
     */
    std::size_t Ns() const { return species_names.size(); }

    /**
     * Version of libantioch reported by Antioch::get_antioch_version
     */
    int antioch_ver;

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


    /**
     * Antioch::ChemicalMixture object, required for anything with
     * antioch
     */
    shared_ptr<Antioch::ChemicalMixture<real_t> > mixture;

    /**
     * Antioch::ReactionSet object, used by Kinetics in evaluating
     * reaction srcs
     */
    shared_ptr<Antioch::ReactionSet<real_t> > reactions;

    /**
     * Antioch::CEAThermodynamics object, used to compute inputs for
     * reaction src evaluation
     */
    shared_ptr<Antioch::CEAThermodynamics<real_t> > cea_thermo;

    /**
     * Antioch::StatMechThermodynamics object, used to compute
     * temperature from internal energy
     */
    shared_ptr<Antioch::StatMechThermodynamics<real_t> > sm_thermo;

    /**
     * Antioch::KineticsEvaluator object, used to actually compute
     * mass sources due to chemical reactions.
     *
     * FIXME: This object used in this way probably makes this class
     * thread unsafe.  See comments in antioch/kinetics_evaluator.h
     * and issue #2798.
     */
    shared_ptr<Antioch::KineticsEvaluator<real_t> > kinetics;

    /**
     * Antioch::MixtureViscosity object, used to compute species
     * viscosities using Blottner curve fits.
     */
    shared_ptr<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<real_t>,
                                real_t> > mixture_mu;

    /**
     * Antioch::EuckenThermalConductivity object, used to compute
     * species thermal conductivity according to Eucken model.
     */
    shared_ptr<Antioch::EuckenThermalConductivity<
                   Antioch::StatMechThermodynamics<real_t> > > mixture_kappa;

    /**
     * Antioch::WilkeMixture object, used by Antioch::WilkeEvaluator
     */
    shared_ptr<Antioch::WilkeMixture<real_t> > wilke_mixture;

    /**
     * Antioch::WilkeMixture object, computes mixture transport
     * properties
     */
    shared_ptr<Antioch::WilkeEvaluator<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<real_t>, real_t>,
                                       Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<real_t> >,
                                       real_t> > wilke_evaluator;

private:

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    antioch_constitutive(const antioch_constitutive&);
    antioch_constitutive& operator=(const antioch_constitutive&);
};

} // namespace reacting

} // namespace suzerain

#endif // HAVE_ANTIOCH

#endif // SUZERAIN_REACTING_ANTIOCH_CONSTITUTIVE_HPP
