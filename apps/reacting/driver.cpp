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
 * @copydoc driver.hpp
 */

#include "driver.hpp"

#include <esio/esio.h>

#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/profile.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/samples.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "reacting_ndx.hpp"

namespace suzerain {

namespace reacting {

driver::driver(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : driver_base(application_synopsis,
                  argument_synopsis,
                  description,
                  revstr)
    , isothermal(make_shared<support::definition_isothermal>())
    , chdef(make_shared<definition_channel>())
    , cmods(make_shared<antioch_constitutive>())
    , fsdef(make_shared<definition_filter>())
    , sgdef(make_shared<support::definition_largo>())
    , who("reacting")
{
    // Sets up usual 5 fields.  If necessary, species are added later.
    this->fields = default_fields();
}

std::vector<std::string>
driver::initialize(
        int argc,
        char **argv)
{
    // msoln is not used by all binaries and is therefore not added below
    options.add_definition(*chdef);
    options.add_definition(*cmods);
    options.add_definition(*fsdef);
    options.add_definition(*isothermal);
    options.add_definition(*sgdef);

    // Delegate to superclass initialization
    std::vector<std::string> positional = super::initialize(argc, argv);

    return positional;
}

void
driver::reset()
{
    isothermal.reset();
    chdef     .reset();
    msoln     .reset();
    cmods     .reset();
    fsdef     .reset();
    sgdef     .reset();
    return super::reset();
}

void
driver::log_manufactured_solution_absolute_error(
        const std::string& timeprefix,
        const time_type t,
        const step_type nt)
{
    SUZERAIN_UNUSED(nt);

    // NOP whenever a manufactured_solution is not in use
    if (!msoln) return;

    // Avoid computational cost when logging is disabled
    support::logging::logger_type mms_abserr
            = support::logging::get_logger("mms.abserr");
    if (!INFO0_ENABLED(mms_abserr)) return;

    // Compute $L^2_{xyz}$ of error of state against manufactured solution
    state_nonlinear->assign_from(*state_linear);
    accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            *grid, *dgrid, *cop, *b, t);
    const std::vector<field_L2xyz> L2
        = compute_field_L2xyz(*state_nonlinear, *grid, *dgrid, *gop);

    // Output absolute global errors for each field
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(mms_abserr, msg.str());
}

void
driver::compute_statistics(
        driver::time_type t)
{
    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    state_nonlinear->assign_from(*state_linear);
    mean = reacting::sample_quantities(*cmods, *grid, *dgrid,
                                       *cop, *state_nonlinear, t);

    // Obtain mean quantities computed via implicit forcing (when possible)
    if (common_block.means.rows() == mean.storage.rows()) {
        mean.f().col(0)      = common_block.fx();
        mean.f().col(1)      = common_block.fy();
        mean.f().col(2)      = common_block.fz();
        mean.f_dot_u()       = common_block.f_dot_u();
        mean.qb()            = common_block.qb();
        mean.CrhoE()         = common_block.CrhoE();
        mean.C2rhoE()        = common_block.C2rhoE();
        mean.Crhou() .col(0) = common_block.Crhou();
        mean.C2rhou().col(0) = common_block.C2rhou();
        mean.Crhou() .col(1) = common_block.Crhov();
        mean.C2rhou().col(1) = common_block.C2rhov();
        mean.Crhou() .col(2) = common_block.Crhow();
        mean.C2rhou().col(2) = common_block.C2rhow();
        mean.Crho()          = common_block.Crho();
        mean.C2rho()         = common_block.C2rho();
        mean.Crhou_dot_u()   = common_block.Crhou_dot_u();
        mean.C2rhou_dot_u()  = common_block.C2rhou_dot_u();

        // Convert values at collocation points to expansion coefficients
        bsplineop_lu mass(*cop);
        mass.opform_mass(*cop);
        mass.factor();
        mass.solve(quantities::nscalars::implicit,
                   mean.storage.middleCols<quantities::nscalars::implicit>(
                       quantities::start::implicit).data(),
                   mean.storage.innerStride(), mean.storage.outerStride());
    } else {
        WARN0(who, "Could not obtain mean quantities set by implicit forcing");
    }
}

void
driver::save_spectra_primitive(
        const esio_handle esioh)
{
    SUZERAIN_TIMER_SCOPED("driver::save_spectra_primitive");

    const size_t Ns = cmods->Ns();
    const size_t state_count = Ns+4;

    // Indices for auxiliary storage to allow the diluter to
    // be included in the two-point corralation computation
    struct aux { enum {
            T       = 0,
            u       = 1,
            v       = 2,
            w       = 3,
            rho     = 4,
            species = 5
        }; };

    // Get total number of fields in aux storage
    const size_t aux_count = (aux::species + Ns);

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    // Maybe FIXME: In this way, instead of reusing the memory of state_nonlinear,
    // we are allocating new memory for all the variables to include the diluter.
    // There might be a more memory-efficient way to do this...
    typedef contiguous_state<4,complex_t> state_type;
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux_count, *dgrid));                               // RAII
    state_type &auxw = *_auxw_ptr;                                 // Shorthand

    // Ensure state and auxw storage meets this routine's assumptions
    const unsigned int swave_count = state_count;
    assert(static_cast<int>(ndx::e      ) < swave_count);
    assert(static_cast<int>(ndx::mx     ) < swave_count);
    assert(static_cast<int>(ndx::my     ) < swave_count);
    assert(static_cast<int>(ndx::mz     ) < swave_count);
    assert(static_cast<int>(ndx::rho    ) < swave_count);
    if (Ns > 1) {
        assert(static_cast<int>(ndx::species) < swave_count);
    }
    SUZERAIN_ENSURE(state_nonlinear->shape()[0] == swave_count);
    SUZERAIN_ENSURE(std::equal(state_nonlinear->shape() + 1,
                               state_nonlinear->shape() + 4,
                               auxw.shape() + 1));
    SUZERAIN_ENSURE(std::equal(state_nonlinear->strides() + 1,
                               state_nonlinear->strides() + 4,
                               auxw.strides() + 1));

    // Obtain information from instantaneous fields in *state_linear
    state_nonlinear->assign_from(*state_linear);

    // Convert to pointwise conserved state in physical space
    physical_view<> sphys(*dgrid, *state_nonlinear);
    shared_ptr<operator_tools> otool = obtain_operator_tools();
    for (size_t f = 0; f < swave_count; ++f) {
        otool->zero_dealiasing_modes(  *state_nonlinear, f);
        otool->bop_apply     (0,    1, *state_nonlinear, f);
        dgrid->transform_wave_to_physical(&sphys.coeffRef(f,0));
    }

    // Pointwise conversion from e, mx, my, mz, rho, rho_s
    // to T, u, v, w, rho, c_s
    physical_view<> auxp(*dgrid, auxw);
    real_t Tguess = -1;
    VectorXr species(Ns); // species densities
    VectorXr cs     (Ns); // species mass fractions
    for (int offset = 0; offset < sphys.cols(); ++offset) {

        // Unpack conserved state
        const real_t   e  (sphys(ndx::e,   offset));
        const Vector3r m  (sphys(ndx::mx,  offset),
                           sphys(ndx::my,  offset),
                           sphys(ndx::mz,  offset));
        const real_t   rho(sphys(ndx::rho, offset));

        // Compute desired primitive state
        const Vector3r u = m / rho;
        real_t T = 0;

        // Unpack species variables
        // NOTE: In species vector, idx 0 is the dilluter (the species
        // that is not explicitly part of the state vector)
        // Compute mass fractions
        const real_t irho = 1.0/rho;
        species(0) = rho;
        for (unsigned int s=1; s<Ns; ++s) {
            species(s) = sphys(ndx::species + s - 1, offset);
            cs(s) = irho * species(s);

            // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
            species(0)        -= species(s);
        }
        cs(0) = irho * species(0);

        // Compute temperature from internal energy
        // (assuming thermal equilibrium)
        const real_t re_internal = e -
            0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        T = cmods->sm_thermo->T_from_e_tot(irho*re_internal, cs, Tguess);
        Tguess = T;

        // Pack desired primitive state
        auxp(aux::T  , offset) = T;
        auxp(aux::u  , offset) = u.x();
        auxp(aux::v  , offset) = u.y();
        auxp(aux::w  , offset) = u.z();
        auxp(aux::rho, offset) = rho;
        for (unsigned int s=0; s<Ns; ++s) {
            auxp(aux::species + s, offset) = cs(s);
        }
    }

    // TODO Zero dealiasing modes?  If so, can't recover m = rho * u.
    // Convert to normalized Fourier coefficients in XZ but points in Y
    // Keep means, compute raw statistics
    auxp *= dgrid->chi();
    for (size_t f = 0; f < aux_count; ++f) {
        dgrid->transform_physical_to_wave(&auxp.coeffRef(f,0));
    }

    // Only rank zero will write the results though all ranks must participate
    const int npairs = (aux_count*(aux_count+1))/2;
    int procid;
    esio_handle_comm_rank(esioh, &procid);

    // Compute and save the two-point correlation as (y_j, k_x, ndxpair)
    // Ordering arises from packing of primitive state in physical space
    {
        shared_array<complex_t> twopoint_x = compute_twopoint_x(
                auxw, aux_count, *grid, *dgrid);
        esio_field_establish(esioh,
                npairs,          0, procid == 0 ? npairs          : 0,
                grid->N.x()/2+1, 0, procid == 0 ? grid->N.x()/2+1 : 0,
                grid->N.y(),     0, procid == 0 ? grid->N.y()     : 0);
        support::complex_field_write(esioh,
                "twopoint_kx", twopoint_x.get(), 0, 0, 0,
                "Streamwise two-point correlations stored row-major"
                " (/collocation_points_y, /kx, scalarpair) for scalarpair"
                " in { T*T, T*u, T*v, T*w, T*rho, T*{c_j}, u*u, u*v, u*w,"
    	        " u*rho, u*{c_j}, v*v, v*w, v*rho, v*{c_j}, w*w, w*rho,"
                " w*{c_j}, rho*rho, rho*{c_j}, {c_i}*{c_j}}");
    }

    // Compute and save the two-point correlation as (y_j, k_z, ndxpair)
    // Ordering arises from packing of primitive state in physical space
    {
        shared_array<complex_t> twopoint_z = compute_twopoint_z(
                auxw, aux_count, *grid, *dgrid);
        esio_field_establish(esioh,
                npairs,          0, procid == 0 ? npairs          : 0,
                grid->N.z(),     0, procid == 0 ? grid->N.z()     : 0,
                grid->N.y(),     0, procid == 0 ? grid->N.y()     : 0);
        support::complex_field_write(esioh,
                "twopoint_kz", twopoint_z.get(), 0, 0, 0,
                "Spanwise two-point correlations stored row-major"
                " (/collocation_points_y, /kz, scalarpair) for scalarpair"
                " in { T*T, T*u, T*v, T*w, T*rho, T*{c_j}, u*u, u*v, u*w,"
                " u*rho, u*{c_j}, v*v, v*w, v*rho, v*{c_j}, w*w, w*rho,"
                " w*{c_j}, rho*rho, rho*{c_j}, {c_i}*{c_j}}");
    }
}

bool
driver::log_status_hook(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    const bool retval = super::log_status_hook(timeprefix, t, nt);
    log_manufactured_solution_absolute_error(timeprefix, t, nt);
    return retval;
}

void
driver::save_metadata_hook(
        const esio_handle esioh)
{
    super::save_metadata_hook(esioh);
    chdef->save(esioh);
    cmods->save(esioh);
    fsdef->save(esioh);
    fsdef->save_filteropz(esioh, dgrid->global_wave_extent.y());
    isothermal->save(esioh);
    sgdef->save(esioh);
    save(esioh, msoln, *cmods, *grid);
    return;
}

void
driver::load_metadata_hook(
        const esio_handle esioh)
{
    super::load_metadata_hook(esioh);
    chdef->load(esioh);
    cmods->load(esioh);
    fsdef->load(esioh);
    isothermal->load(esioh);
    sgdef->load(esioh);

    // After cmods->load, have valid species names, so we can update
    // fields to include species
    add_species_fields(cmods->species_names, this->fields);

    load(esioh, msoln, *cmods, *grid);
    return;
}

bool
driver::save_state_hook(
        const esio_handle esioh)
{
    // Get the state to disk as quickly as possible
    const bool success = super::save_state_hook(esioh);

    // If that went well, go for the spectra
    if (success) {
        save_spectra_primitive(esioh);
    }

    return success;
}

bool
driver::save_statistics_hook(
        const esio_handle esioh,
        const driver::time_type t)
{
    DEBUG0(who, "In driver::save_statistics_hook");

    // Should we compute fresh statistics or re-use cached values?
#pragma warning(push,disable:1572)
    const bool use_cached = controller && (t == mean.t);
#pragma warning(pop)

    // Compute a prefix for logging purposes, including a trailing space
    std::string prefix;
    if (controller) {
        prefix = build_timeprefix(t, controller->current_nt());
        prefix.append(1, ' ');
    }

    // Compute statistics whenever necessary
    if (use_cached) {
        DEBUG0(who, prefix << "Cowardly refusing to re-compute statistics");
    } else {
        const double starttime = MPI_Wtime();
        compute_statistics(t);
        const double elapsed = MPI_Wtime() - starttime;
        INFO0(who, prefix << "Computed statistics in "
              << elapsed << " seconds");
    }

    // Save statistics and invoke superclass hook
    if (!mean.save(esioh)) {
        WARN0(who, "Incomplete statistics saved at time " << t);
    }
    return super::save_statistics_hook(esioh, t);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    support::load_time(esioh, mean.t);

    if (statsdef && statsdef->stale) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Ignoring incoming statistics explicitly marked as stale");
        mean.t = std::numeric_limits<real_t>::quiet_NaN();
    } else if (!mean.load(esioh)) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Incomplete statistics loaded from time " << mean.t);
        mean.t = std::numeric_limits<real_t>::quiet_NaN();
    }
    return super::load_statistics_hook(esioh);
}

// Heavy-handed, but I seemingly cannot get Intel to cooperator otherwise.
// I want non-NaNs to be preferred to NaNs within this max operation.
static
real_t maxmf(const real_t a, const real_t b)
{
    return (boost::math::isnan)(a) ? b
         : (boost::math::isnan)(b) ? a
         : std::max(a, b);
}

// Please review and possibly synchronize apps/perfect/driver.cpp when changing
void
driver::default_restart_interval(
        time_type& t,
        step_type&)
{
    using std::abs;

    // Look for the largest magnitude, problem-dependent, macro velocity scale
    real_t velocity = std::numeric_limits<real_t>::quiet_NaN();

    // In a channel...
    if (grid && grid->two_sided()) {

        // ...any approximate bulk velocity from a driving force...
        if (chdef) {

            if ((boost::math::isnan)(chdef->bulk_rho)) {
                WARN0(who, "No bulk density scale available so assuming 1");
                velocity = maxmf(abs(chdef->bulk_rho_u) / /*rho*/ 1,
                                 velocity);
            } else {
                velocity = maxmf(abs(chdef->bulk_rho_u) / chdef->bulk_rho,
                                 velocity);
            }

        }

        // ...may be trumped by driving the upper and lower walls.
        if (isothermal) {
            velocity = maxmf(abs(isothermal->upper_u), velocity);
            velocity = maxmf(abs(isothermal->lower_u), velocity);
        }

    }

    // On a plate...
    if (grid && grid->one_sided()) {

        if (sgdef->formulation.enabled() && sgdef->baseflow && fields.size()) {

            // ...prefer wall velocity from inviscid baseflow as freestream
            // surrogate as it is independent of the wall-normal domain size...
            const size_t nvalues = fields.size() + /*pressure*/1;
            std::vector<real_t> buf(3 * nvalues);
            sgdef->baseflow->conserved(0.0,               // wall
                                       &buf[0*nvalues],   // base
                                       &buf[1*nvalues],   // dybase
                                       &buf[2*nvalues]);  // dxbase
            const real_t base_rho   = buf[0];
            const real_t base_rho_u = buf[1];
            if (base_rho > 0) {
                velocity = maxmf(abs(base_rho_u) / base_rho, velocity);
            } else {
                TRACE0(who, "Trivial baseflow detected; using upper velocity");
                velocity = maxmf(abs(isothermal->upper_u), velocity);
            }

        } else if (isothermal) {
            // ...taking freestream reference if-and-only-if no baseflow...
            velocity = maxmf(abs(isothermal->upper_u), velocity);
        }

        // ...and permit a driven lower wall velocity to trump.
        if (isothermal) {
            velocity = maxmf(abs(isothermal->lower_u), velocity);
        }

    }

    // Have a domain size and a velocity?  Let's use the matching time scale.
    if (boost::math::isnormal(velocity)) {
        if (grid) {
            t = grid->L.x() / velocity;
            if        (grid->two_sided()) {
                // Eight restarts per flow through in a channel giving 80
                // restarts in 10 flow throughs, a reasonable channel simulation
                // duration.
                t /= 8;
            } else if (grid->one_sided()) {
                // Two restarts per flow through in a boundary layer giving 20
                // restarts in 10 flow throughs and, across O(50) flow throughs
                // for a simulation campaign, O(100) restart files.
                t /= 2;
            } else {
                SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
            }
        } else {
            DEBUG0(who, "No grid details for default_restart_interval");
        }
    } else {
        DEBUG0(who, "No macro scale velocity for default_restart_interval");
    }
}

} // end namespace reacting

} // end namespace suzerain
