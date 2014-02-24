//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

/** @file
 * @copydoc manufactured_solution.hpp
 */

#include "manufactured_solution.hpp"

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/blas_et_al.hpp>
#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

namespace suzerain {

namespace reacting {

const std::string manufactured_solution::default_caption(
        "Manufactured solution parameters");

manufactured_solution::manufactured_solution()
    : caption(default_caption)
{
}

manufactured_solution::manufactured_solution(
        const std::string& caption)
    : caption(caption)
{
}

/** Helper for manufactured_solution::finish_construction(). */
static void option_adder(
        boost::program_options::options_description_easy_init easy_init,
        const char *description,
        const std::string &name,
        real_t &value)
{
    easy_init(name.c_str(),
              boost::program_options::value(&value)->default_value(value),
              description);
}

// Strings used in options_description and save/load.
static const char name_T_r[]  = "T_r";
static const char name_mu_r[] = "mu_r";
static const char name_beta[] = "beta";

// Descriptions used in options_description and save/load.
static const char desc_T_r[]  = "Reference temperature for viscosity power law (for MMS)";
static const char desc_mu_r[] = "Reference viscosity for viscosity power law (for MMS)";
static const char desc_beta[] = "Exponent for viscosity power law (for MMS)";

boost::program_options::options_description
manufactured_solution::options_description()
{
    using boost::bind;
    using boost::program_options::value;

    boost::program_options::options_description retval(caption);

    rho.foreach_parameter(bind(option_adder, retval.add_options(),
                               "Affects density field",     _1, _2));
    u.foreach_parameter  (bind(option_adder, retval.add_options(),
                               "Affects streamwise velocity field",  _1, _2));
    v.foreach_parameter  (bind(option_adder, retval.add_options(),
                               "Affects wall-normal velocity field",  _1, _2));
    w.foreach_parameter  (bind(option_adder, retval.add_options(),
                               "Affects spanwise velocity field",  _1, _2));
    T.foreach_parameter  (bind(option_adder, retval.add_options(),
                               "Affects temperature field", _1, _2));

    // TODO: Check validity of incoming values
    retval.add_options()(name_T_r,
                         value(&this->T_r)->default_value(273.0),
                         desc_T_r);

    retval.add_options()(name_mu_r,
                         value(&this->mu_r)->default_value(1852./100000000.),
                         desc_mu_r);

    retval.add_options()(name_beta,
                         value(&this->beta)->default_value(2./3.),
                         desc_beta);

    return retval;
}

manufactured_solution&
manufactured_solution::match(antioch_constitutive& cmods)
{
    // Check incoming cmods for consistency with available
    // manufactured solution
    if (cmods.Ns()!=1)
        throw std::invalid_argument("MMS only available for single species!");

    if (cmods.species_names[0]!="CPAir" && cmods.species_names[0]!="CPN2")
        throw std::invalid_argument(
            "MMS only available for calorically perfect gas!");

    // Tell user we're going to mess with their incoming cmods
    INFO0("Manipulating constitutive model parameters to "
          "match input manufactured solution.");

    // Since only allow single species, mass fractions are trivial
    std::vector<real_t> trivial_mass_fractions(1, 1.0);

    // Thermodynamics
    // Note: T, Tv, and mass_fractions inputs
    // Temperatures are irrelevant (calorically perfect) but mass
    // fraction input is necessary
    real_t Cv = cmods.sm_thermo->cv(1.0,1.0,trivial_mass_fractions);
    real_t Cp = cmods.sm_thermo->cp(1.0,1.0,trivial_mass_fractions);

    real_t Rmix = cmods.mixture->R(trivial_mass_fractions);

    this->gamma = Cp/Cv;
    this->R     = Rmix;

    // Transport
    // Note: Parameters of transport models are stored as part of MMS
    // class.  Here we update the constitutive laws class to be
    // consistent with what is stored by MMS.

    // Viscosity: Make Blottner match power law
    const real_t mu0 = this->mu_r;
    const real_t beta = this->beta;
    const real_t Tref = this->T_r;

    std::vector<real_t> blottner_coeffs(3);
    // Blottner says: mu = 0.1*exp( a*log(T)^2 + b*log(T) + c).  Thus,

    // for a = 0 reduces blottner to power law
    blottner_coeffs[0] = 0.0;

    // with b = power law exponent
    blottner_coeffs[1] = beta;

    // and 0.1*exp(c) = mu0/(Tref^beta)
    blottner_coeffs[2] = std::log(mu0/(0.1*std::pow(Tref, beta)));

    // Set blottner parameters
    cmods.mixture_mu->reset_coeffs(0, blottner_coeffs);

    // Thermal conductivity: Weird factor of 19/4 makes
    // man. soln. match Eucken for CPAir and CPN2
    this->kappa_r  = (real_t(19)/real_t(4))*this->mu_r*this->R;

    // Bulk viscosity
    this->lambda_r = (cmods.alpha-(real_t(2)/real_t(3)))*this->mu_r;

    return *this;
}

manufactured_solution&
manufactured_solution::match(const specification_grid& grid)
{
    this->Lx = grid.L.x();
    this->Ly = grid.L.y();
    this->Lz = grid.L.z();
    return *this;
}

manufactured_solution&
manufactured_solution::isothermal_channel()
{
    //nsctpl_rholut::isothermal_channel(*this);
    nsctpl::isothermal_channel(*this);

    // Rescale time frequencies and density coefficients
    const real_t rho_scale  = real_t(355)/real_t(10000000);
    const real_t itime_scale = 200;

    // density coefficients
    this->rho.a_0  = rho_scale * real_t(1);
    this->rho.a_xy = rho_scale * real_t(1) / real_t(11);
    this->rho.a_y  = rho_scale * real_t(1) / real_t(7);
    this->rho.a_yz = rho_scale * real_t(1) / real_t(31);

    // time frequency coefficients
    this->rho.f_xy = itime_scale * real_t(3);
    this->rho.f_y  = itime_scale * real_t(1);
    this->u.f_xy   = itime_scale * real_t(3);
    this->u.f_y    = itime_scale * real_t(1);
    this->u.f_yz   = itime_scale * real_t(2);
    this->v.f_y    = itime_scale * real_t(1);
    this->v.f_yz   = itime_scale * real_t(2);
    this->w.f_xy   = itime_scale * real_t(3);
    this->w.f_y    = itime_scale * real_t(1);
    this->w.f_yz   = itime_scale * real_t(2);
    this->T.f_xy   = itime_scale * real_t(3);
    this->T.f_y    = itime_scale * real_t(1);
    this->T.f_yz   = itime_scale * real_t(2);

    return *this;
}

manufactured_solution&
manufactured_solution::isothermal_flat_plate()
{
    //nsctpl_rholut::isothermal_flat_plate(*this);
    nsctpl::isothermal_flat_plate(*this);

    // Rescale time frequencies and density coefficients
    const real_t rho_scale  = real_t(355)/real_t(10000000);
    const real_t itime_scale = 200;

    // density coefficients
    this->rho.a_0  = rho_scale * real_t(1);
    this->rho.a_xy = rho_scale * real_t(1) / real_t(11);
    this->rho.a_y  = rho_scale * real_t(1) / real_t(7);
    this->rho.a_yz = rho_scale * real_t(1) / real_t(31);

    // time frequency coefficients
    this->rho.f_xy = itime_scale * real_t(3);
    this->rho.f_y  = itime_scale * real_t(1);
    this->u.f_xy   = itime_scale * real_t(3);
    this->u.f_y    = itime_scale * real_t(1);
    this->u.f_yz   = itime_scale * real_t(2);
    this->v.f_y    = itime_scale * real_t(1);
    this->v.f_yz   = itime_scale * real_t(2);
    this->w.f_xy   = itime_scale * real_t(3);
    this->w.f_y    = itime_scale * real_t(1);
    this->w.f_yz   = itime_scale * real_t(2);
    this->T.f_xy   = itime_scale * real_t(3);
    this->T.f_y    = itime_scale * real_t(1);
    this->T.f_yz   = itime_scale * real_t(2);

    return *this;
}

/** Helper for save(...,shared_ptr<manufactured_solution>&,...). */
static void
attribute_storer(const esio_handle &h,
                 const char *location,
                 const std::string &name,
                 const real_t &value)
{
    esio_attribute_write(h, location, name.c_str(), &value);
}

void save(const esio_handle h,
          const shared_ptr<manufactured_solution>& msoln,
          const antioch_constitutive& cmods,
          const specification_grid& grid,
          const char *location)
{
    SUZERAIN_UNUSED(cmods);  // TODO?

    // Only proceed if a manufactured solution is being provided
    if (!msoln) return;

    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
            "Is a suzerain::reacting::manufactured_solution in use?");

    DEBUG0("Storing manufactured_solution parameters");

    // Non-scenario solution parameters are stored as attributes under location
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->u  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->v  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->w  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->T  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));

    // Scenario params
    esio_line_write(h, name_T_r,  &msoln->T_r,  0, desc_T_r );
    esio_line_write(h, name_mu_r, &msoln->mu_r, 0, desc_mu_r);
    esio_line_write(h, name_beta, &msoln->beta, 0, desc_beta);

#pragma warning(push,disable:1572)
    // Check parameters stored with the scenario not the manufactured solution
    // because scenario parameters should be loaded from definition_scenario

    // FIXME: Add this check back!!!

    // Check parameters stored with the grid not the manufactured solution
    // because grid parameters should be loaded from specification_grid
    if (msoln->Lx    != grid.L.x())
        WARN0("Manufactured solution Lx mismatches with grid");
    if (msoln->Ly    != grid.L.y())
        WARN0("Manufactured solution Ly mismatches with grid");
    if (msoln->Lz    != grid.L.z())
        WARN0("Manufactured solution Lz mismatches with grid");
#pragma warning(pop)
}

/** Helper for load(...,shared_ptr<manufactured_solution>&,...). */
static void
attribute_loader(const esio_handle &h,
                 const char *location,
                 const std::string &name,
                 real_t &value)
{
    esio_attribute_read(h, location, name.c_str(), &value);
}

/** Helper for NaNing values within a \ref manufactured_solution. */
static void NaNer(const std::string&, real_t& value)
{
    value = std::numeric_limits<real_t>::quiet_NaN();
}

void load(const esio_handle h,
          shared_ptr<manufactured_solution>& msoln,
          const antioch_constitutive& cmods,
          const specification_grid& grid,
          const char *location)
{
    SUZERAIN_UNUSED(cmods);  // TODO?

    // Only proceed if a manufactured solution is active in the restart
    int in_use = 0;
    esio_line_establish(h, 1, 0, 1); // All ranks load any data
    if (ESIO_NOTFOUND != esio_line_size(h, location, NULL)) {
        esio_line_read(h, location, &in_use, 0);
    }
    if (!in_use) {
        msoln.reset();
        return;
    }

    DEBUG0("Loading manufactured_solution parameters");

    // Allocate storage and defensively NaN every parameter not explicitly
    // loaded below.  Protects us against accidentally missing new solution
    // parameters.
    msoln.reset(new manufactured_solution());
    msoln->foreach_parameter(&NaNer);

    // Non-scenario solution parameters are stored as attributes under location
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->u  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->v  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->w  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->T  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));

    // Scenario params
    esio_line_read(h, name_T_r,  &msoln->T_r,  0);
    esio_line_read(h, name_mu_r, &msoln->mu_r, 0);
    esio_line_read(h, name_beta, &msoln->beta, 0);

    // Parameters set to match supplied arguments
    // Match must happen *after* call to antioch_constitutive::init_antioch
    //msoln->match(cmods);
    msoln->match(grid);
}

void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        contiguous_state<4,complex_t> &swave,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        const real_t simulation_time)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                         swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(boost::numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(boost::numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(boost::numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Initialize operator_base to access decomposition-ready utilities
    operator_base obase(grid, dgrid, cop, b);

    // Allocate one field of temporary storage for scratch purposes
    scoped_ptr<contiguous_state<4,complex_t> > _scratch_ptr(        // RAII
        support::allocate_padded_state<contiguous_state<4,complex_t> >(1,dgrid));
    contiguous_state<4,complex_t> &scratch = *_scratch_ptr;         // Shorthand
    multi_array::fill(scratch, 0);                                  // Defensive

    // Prepare physical-space view of the wave-space scratch storage
    physical_view<1> phys(dgrid, scratch);

    // Prepare factored mass matrix for repeated use
    bsplineop_luz massluz(cop);
    const complex_t scale_factor = 1 / dgrid.chi();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();

    // For each scalar field...
    for (size_t f = 0; f < state_count; ++f) {

        // ...compute the manufactured solution in physical space...
        for (int offset = 0, j = dgrid.local_physical_start.y();
            j < dgrid.local_physical_end.y();
            ++j) {

            const real_t y = obase.y(j);

            for (int k = dgrid.local_physical_start.z();
                k < dgrid.local_physical_end.z();
                ++k) {

                const real_t z = obase.z(k);

                for (int i = dgrid.local_physical_start.x();
                    i < dgrid.local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    const real_t x = obase.x(i);

                    // Ugly, slow switch but performance irrelevant here
                    switch (f) {
                    case ndx::e:
                        phys(0, offset) = msoln.rhoe(x, y, z, simulation_time);
                        break;
                    case ndx::mx:
                        phys(0, offset) = msoln.rhou(x, y, z, simulation_time);
                        break;
                    case ndx::my:
                        phys(0, offset) = msoln.rhov(x, y, z, simulation_time);
                        break;
                    case ndx::mz:
                        phys(0, offset) = msoln.rhow(x, y, z, simulation_time);
                        break;
                    case ndx::rho:
                        phys(0, offset) = msoln.rho (x, y, z, simulation_time);
                        break;
                    default:
                        SUZERAIN_ERROR_REPORT("unknown field",
                                              SUZERAIN_ESANITY);
                    }

                } // end X

            } // end Z

        } // end Y

        // ...convert collocation values to non-dealiased coefficients...
        dgrid.transform_physical_to_wave(&phys.coeffRef(0, 0));  // X, Z
        obase.zero_dealiasing_modes(scratch, 0);
        obase.bop_solve(massluz, scratch, 0);                    // Y

        // ...and accumulate into the corresponding scalar field of swave
        SUZERAIN_ENSURE(std::equal(scratch.shape() + 1, scratch.shape() + 4,
                                   swave.shape() + 1));
        if (SUZERAIN_UNLIKELY(0U == scratch.shape()[1])) {
            continue;  // Sidestep assertions on trivial data
        }
        typedef contiguous_state<4,complex_t>::index index;
        const index ku = boost::numeric_cast<index>(
                                scratch.index_bases()[2] + scratch.shape()[2]);
        const index lu = boost::numeric_cast<index>(
                                scratch.index_bases()[3] + scratch.shape()[3]);
        for (index lx = scratch.index_bases()[3], ly = swave.index_bases()[3];
            lx < lu;
            ++lx, ++ly) {

            for (index kx = scratch.index_bases()[2],
                       ky = swave.index_bases()[2];
                kx < ku;
                ++kx, ++ky) {

                blas::axpby(scratch.shape()[1], alpha,
                        &scratch[0][scratch.index_bases()[1]][kx][lx],
                        scratch.strides()[1], beta,
                        &swave[f][swave.index_bases()[1]][ky][ly],
                        swave.strides()[1]);
            }
        }

    }
}

} // end namespace reacting

} // end namespace suzerain
