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
#include <largo/largo.h>

#include <suzerain/bl.h>
#include <suzerain/channel.h>
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

#include "instantaneous.hpp"
#include "perfect.hpp"

namespace suzerain {

namespace perfect {

driver::driver(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : driver_base(application_synopsis,
                  argument_synopsis,
                  description,
                  revstr)
    , scenario(make_shared<definition_scenario>())
    , isothermal(make_shared<support::definition_isothermal>())
    , sg(make_shared<support::definition_largo>())
    , rad(make_shared<support::definition_radialflow>(/*deltae*/ 1.0))
    , who("perfect")
{
    this->fields = default_fields();
}

driver::~driver()
{
    // Deallocate any lazily allocated slow growth forcing workspace
    if (sg.unique()) {
        DEBUG0(who, "Deallocating largo_workspace");
        if (sg->workspace) {
            largo_deallocate(&sg->workspace);
            sg->workspace = NULL;
        }
    } else {
        DEBUG0(who, "Cowardly refusing to deallocate in-use largo_workspace");
    }
}

std::vector<std::string>
driver::initialize(
        int argc,
        char **argv)
{
    // Only add groups of options when non-trivial at initialization
    if (scenario  ) options.add_definition(*scenario  );
    if (isothermal) options.add_definition(*isothermal);
    if (sg        ) options.add_definition(*sg        );
    if (rad       ) options.add_definition(*rad       );
    if (helm      ) options.add_definition(*helm      );
    if (msoln     ) options.add_definition(*msoln     );

    // Delegate to superclass initialization
    std::vector<std::string> positional = super::initialize(argc, argv);

    // If msoln was provided, match its contents to other members
    if (msoln) {
        if (scenario) msoln->match(*scenario);
        if (grid)     msoln->match(*grid);
    }

    return positional;
}

void
driver::reset()
{
    scenario.reset();
    isothermal.reset();
    sg.reset();
    rad.reset();
    helm.reset();
    msoln.reset();
    mean.reset();
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
driver::log_linearization_error(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    SUZERAIN_UNUSED(nt);

    // Avoid computational cost when logging is disabled
    support::logging::logger_type lin_abserr
            = support::logging::get_logger("lin.abserr");
    if (!INFO0_ENABLED(lin_abserr)) return;

    SUZERAIN_TIMER_SCOPED("driver::log_linearization_error");

    // Algorithm proceeds directly from necessary interface guarantees.
    // Rewriting it symbolically from the code and grokking it is a good
    // exercise towards understanding the low storage time interfaces.
    //
    // Success requires that the following preconditions all hold:
    //    i) At substep zero, this->N computes N(u) - Lu including
    //       the necessary reference quantities.
    //   ii) The reference quantities are stored in common_block and may be
    //       expressly zeroed.  Zeroing them nukes linearized contributions
    //       from I{Nonlinear,Linear}Operator but not linear contributions
    //       from this->L.
    //  iii) On non-zero substeps for zero reference values, this->N computes
    //       only N(u) and this->L computes only
    //       \left(M+\varphi\mathscr{L}\right)u.
    //   iv) At the beginning of this process, state_linear contains
    //       valid information and adheres to boundary conditions.
    //    v) The explicit_operator_nonlinear application requires an auxiliary
    //       scaling factor "chi" to account for Fourier transform
    //       normalization needs.
    //   vi) Using diffwave::apply zeros wavenumbers used only for dealiasing
    //       purposes-- this data is unimportant from the perspective of
    //       measuring actual linearization error.
    state_nonlinear->assign_from(*state_linear);
    common_block.set_zero(grid->dN.y());  // Defensive
    N->apply_operator(t, *state_nonlinear, *method, /*substep*/0);
    L->accumulate_mass_plus_scaled_operator(
            1, *state_linear, dgrid->chi(), *state_nonlinear, /*substep*/0);
    common_block.set_zero(grid->dN.y());  // Zero reference quantities
    L->accumulate_mass_plus_scaled_operator(
            1, *state_linear, -1, *state_nonlinear, /*substep*/0);
    for (size_t k = 0; k < fields.size(); ++k) {
        diffwave::apply(0, 0, 1, (*state_nonlinear)[k].origin(),
            grid->L.x(), grid->L.z(), dgrid->global_wave_extent.y(),
            grid->N.x(), grid->dN.x(),
            dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
            grid->N.z(), grid->dN.z(),
            dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
    }
    state_nonlinear->exchange(*state_linear);
    N->apply_operator(t, *state_nonlinear, *method, /*substep*/1);
    for (size_t k = 0; k < fields.size(); ++k) {
        diffwave::apply(0, 0, 1, (*state_nonlinear)[k].origin(),
            grid->L.x(), grid->L.z(), dgrid->global_wave_extent.y(),
            grid->N.x(), grid->dN.x(),
            dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
            grid->N.z(), grid->dN.z(),
            dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
    }
    state_nonlinear->add_scaled(1/dgrid->chi(), *state_linear);

    // When this->L and this->N match perfectly, state_nonlinear should be
    // identically zero at each collocation point.  Seeing more than acceptable
    // (e.g. 1e-10) floating point error indicates something is amiss.
    shared_ptr<operator_tools> otool = obtain_operator_tools();
    for (size_t f = 0; f < state_nonlinear->shape()[0]; ++f) {
        otool->bop_solve(*otool->massluz(), *state_nonlinear, f);
    }
    const std::vector<field_L2xyz> L2
        = compute_field_L2xyz(*state_nonlinear, *grid, *dgrid, *gop);
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(lin_abserr, msg.str());
}

void driver::log_quantities_of_interest(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    SUZERAIN_UNUSED(nt);

    if (!grid) return;

    SUZERAIN_TIMER_SCOPED("driver::log_quantities_of_interest");

    // If possible, use existing information from mean quantities
    // Otherwise compute from instantaneous fields stored in state_linear
    shared_ptr<operator_tools> otool = obtain_operator_tools();
    profile prof;
    instantaneous inst;
    if (controller && mean && mean->t == t) {
        prof = *mean;
        inst.copy_from(*cop, *mean);
    } else {
        state_nonlinear->assign_from(*state_linear);
        prof = *take_profile(*scenario, *otool, *state_nonlinear, inst);
    }

    this->log_quantities_of_interest(timeprefix, prof, inst);
}

void
driver::log_quantities_of_interest(
        const std::string& prefix,
        const profile& prof,
        const instantaneous& inst)
{
    // A horrible way to cheaply grab a shared pointer to a masslu
    // Uses internal knowledge of obtain_operator_tools requiring a dgrid
    shared_ptr<const bsplineop_lu> masslu;
    if (dgrid) {
        masslu = obtain_operator_tools()->masslu();
    } else {
        shared_ptr<bsplineop_lu> t(new bsplineop_lu(*cop));
        t->factor_mass(*cop);
        masslu = t;
    }

    // Fluctuations only make sense to log when either X or Z is nontrivial
    const bool nontrivial_fluct_possible = grid->N.x() * grid->N.z() > 1;
    if (nontrivial_fluct_possible) {

        // Compute pointwise Reynolds-averaged and Favre-averaged fluctuation
        // profiles as three contributions to the TKE production term.
        ArrayX6r reynolds;
        ArrayX6r favre;
        ArrayX3r prodterms;
        compute_fluctuations(*cop, *masslu, inst, reynolds, favre, prodterms);

        // Find coefficient representations from collocation point values
        masslu->solve(reynolds.cols(), reynolds.data(),
                      reynolds.innerStride(), reynolds.outerStride());
        masslu->solve(favre.cols(), favre.data(),
                      favre.innerStride(), favre.outerStride());
        masslu->solve(prodterms.cols(), prodterms.data(),
                      prodterms.innerStride(), prodterms.outerStride());

        // Reduce coefficients into bulk values
        const VectorXr bulk = support::compute_bulk_weights(*b, *masslu);
        {
            // Broken? reynolds .matrix().applyOnTheLeft(bulk.transpose());
            ArrayXXr tmp;
            tmp = bulk.transpose() * reynolds .matrix();  reynolds  = tmp;
            tmp = bulk.transpose() * favre    .matrix();  favre     = tmp;
            tmp = bulk.transpose() * prodterms.matrix();  prodterms = tmp;
        }
        const real_t production = prodterms.sum();

    }

    // Prepare and log geometry-specific quantities of interest
    if (grid->one_sided()) {

        // Compute many details about the boundary layer for logging
        suzerain_bl_local       wall;
        suzerain_bl_viscous     viscous;
        suzerain_bl_local       edge;
        suzerain_bl_local       edge99;
        suzerain_bl_thicknesses thick;
        suzerain_bl_reynolds    reynolds;
        suzerain_bl_qoi         qoi;
        suzerain_bl_pg          pg;
        summarize_boundary_layer_nature(prof, *scenario, sg, *masslu, *b,
                                        wall, viscous, thick,
                                        edge, edge99,
                                        reynolds, qoi, pg);

        // Log messages using application-agnostic superclass functionality
        this->log_quantities_boundary_layer(prefix,
                                            &wall, &viscous, &thick,
                                            &edge, &edge99,
                                            &reynolds, &qoi, &pg);

    } else if (grid->two_sided()) {

        // Compute many details about the boundary layer for logging
        suzerain_channel_local   wall;
        suzerain_channel_viscous viscous;
        suzerain_channel_local   center;
        suzerain_channel_qoi     qoi;
        summarize_channel_nature(prof, *scenario, *b,
                                 wall, viscous, center, qoi);

        // Log messages using application-agnostic superclass functionality
        this->log_quantities_channel(prefix,
                                     &wall, &viscous, &center, &qoi);

    } else {

        SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();

    }
}

void
driver::compute_statistics(
        driver::time_type t)
{
    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    state_nonlinear->assign_from(*state_linear);
    shared_ptr<operator_tools> otool = obtain_operator_tools();
    mean = take_samples(*scenario, *sg, *otool, *b, *state_nonlinear, t);

    // Obtain mean quantities computed via implicit forcing (when possible)
    if (common_block.imp.rows() == mean->storage.rows()) {
        mean->f().col(0)      = common_block.imp.fx();
        mean->f().col(1)      = common_block.imp.fy();
        mean->f().col(2)      = common_block.imp.fz();
        mean->f_dot_u()       = common_block.imp.f_dot_u();
        mean->qb()            = common_block.imp.qb();
        mean->CrhoE()         = common_block.imp.CrhoE();
        mean->C2rhoE()        = common_block.imp.C2rhoE();
        mean->Crhou() .col(0) = common_block.imp.Crhou();
        mean->C2rhou().col(0) = common_block.imp.C2rhou();
        mean->Crhou() .col(1) = common_block.imp.Crhov();
        mean->C2rhou().col(1) = common_block.imp.C2rhov();
        mean->Crhou() .col(2) = common_block.imp.Crhow();
        mean->C2rhou().col(2) = common_block.imp.C2rhow();
        mean->Crho()          = common_block.imp.Crho();
        mean->C2rho()         = common_block.imp.C2rho();
        mean->Crhou_dot_u()   = common_block.imp.Crhou_dot_u();
        mean->C2rhou_dot_u()  = common_block.imp.C2rhou_dot_u();

        // Convert values from collocation points to B-spline coefficients
        bsplineop_lu mass(*cop);
        mass.opform_mass(*cop);
        mass.factor();
        mass.solve(mean->implicit().cols(),
                   mean->implicit().data(),
                   mean->implicit().innerStride(),
                   mean->implicit().outerStride());
    } else {
        WARN0(who, "Could not obtain mean quantities set by implicit forcing");
    }
}

void
driver::save_spectra_primitive(
        const esio_handle esioh)
{
    SUZERAIN_TIMER_SCOPED("driver::save_spectra_primitive");

    enum { swave_count = 5 };
    assert(static_cast<int>(ndx::e  ) < swave_count);
    assert(static_cast<int>(ndx::mx ) < swave_count);
    assert(static_cast<int>(ndx::my ) < swave_count);
    assert(static_cast<int>(ndx::mz ) < swave_count);
    assert(static_cast<int>(ndx::rho) < swave_count);
    SUZERAIN_ENSURE((int) state_nonlinear->shape()[0] == swave_count);

    // Obtain information from instantaneous fields in *state_linear
    state_nonlinear->assign_from(*state_linear);

    // Convert to pointwise conserved state in physical space
    physical_view<swave_count> sphys(*dgrid, *state_nonlinear);
    shared_ptr<operator_tools> otool = obtain_operator_tools();
    for (size_t f = 0; f < swave_count; ++f) {
        otool->zero_dealiasing_modes(  *state_nonlinear, f);
        otool->bop_apply     (0,    1, *state_nonlinear, f);
        dgrid->transform_wave_to_physical(&sphys.coeffRef(f,0));
    }

    // Pointwise conversion from e, mx, my, mz, rho to T, u, v, w, rho
    const real_t alpha = scenario->alpha;
    const real_t beta  = scenario->beta;
    const real_t gamma = scenario->gamma;
    const real_t Ma    = scenario->Ma;
    for (int offset = 0; offset < sphys.cols(); ++offset) {

        // Unpack conserved state
        const real_t   e  (sphys(ndx::e,   offset));
        const Vector3r m  (sphys(ndx::mx,  offset),
                           sphys(ndx::my,  offset),
                           sphys(ndx::mz,  offset));
        const real_t   rho(sphys(ndx::rho, offset));

        // Compute desired primitive state
        const Vector3r u = m / rho;
        real_t p, T;
        rholut::p_T(alpha, beta, gamma, Ma, rho, m, e, p, T);

        // Pack desired primitive state
        sphys(0, offset) = T;
        sphys(1, offset) = u.x();
        sphys(2, offset) = u.y();
        sphys(3, offset) = u.z();
        sphys(4, offset) = rho;
    }

    // Dealiasing modes not zeroed else can't recover m = rho * u.
    // Convert to normalized Fourier coefficients in XZ but points in Y
    // Leave mean else ensemble of two-point in physical space is wrong!
    sphys *= dgrid->chi();
    for (size_t f = 0; f < swave_count; ++f) {
        dgrid->transform_physical_to_wave(&sphys.coeffRef(f,0));
    }

    // Only rank zero will write the results though all ranks must participate
    const int npairs = (swave_count*(swave_count+1))/2;
    int procid;
    esio_handle_comm_rank(esioh, &procid);

    // Compute and save the two-point correlation as (y_j, k_x, ndxpair)
    // Ordering arises from packing of primitive state in physical space
    {
        shared_array<complex_t> twopoint_x = compute_twopoint_x(
                *state_nonlinear, swave_count, *grid, *dgrid);
        esio_field_establish(esioh,
                npairs,          0, procid == 0 ? npairs          : 0,
                grid->N.x()/2+1, 0, procid == 0 ? grid->N.x()/2+1 : 0,
                grid->N.y(),     0, procid == 0 ? grid->N.y()     : 0);
        support::complex_field_write(esioh,
                "twopoint_kx", twopoint_x.get(), 0, 0, 0,
                "Streamwise two-point correlations stored row-major"
                " (/collocation_points_y, /kx, scalarpair) for scalarpair"
                " in { T*T, T*u, T*v, T*w, T*rho, u*u, u*v, u*w, u*rho,"
                " v*v, v*w, v*rho, w*w, w*rho, rho*rho }");
    }

    // Compute and save the two-point correlation as (y_j, k_z, ndxpair)
    // Ordering arises from packing of primitive state in physical space
    {
        shared_array<complex_t> twopoint_z = compute_twopoint_z(
                *state_nonlinear, swave_count, *grid, *dgrid);
        esio_field_establish(esioh,
                npairs,          0, procid == 0 ? npairs          : 0,
                grid->N.z(),     0, procid == 0 ? grid->N.z()     : 0,
                grid->N.y(),     0, procid == 0 ? grid->N.y()     : 0);
        support::complex_field_write(esioh,
                "twopoint_kz", twopoint_z.get(), 0, 0, 0,
                "Spanwise two-point correlations stored row-major"
                " (/collocation_points_y, /kz, scalarpair) for scalarpair"
                " in { T*T, T*u, T*v, T*w, T*rho, u*u, u*v, u*w, u*rho,"
                " v*v, v*w, v*rho, w*w, w*rho, rho*rho }");
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
    log_quantities_of_interest(timeprefix, t, nt);
    return retval;
}

void
driver::save_metadata_hook(
        const esio_handle esioh)
{
    super::save_metadata_hook(esioh);
    if (scenario)   scenario->save(esioh);
    if (isothermal) isothermal->save(esioh);
    if (sg)         sg->save(esioh);
    if (rad)        rad->save(esioh);
    if (helm)       helm->save(esioh);
    if (msoln)      save(esioh, msoln, *scenario, *grid);
    return;
}

void
driver::load_metadata_hook(
        const esio_handle esioh)
{
    super::load_metadata_hook(esioh);

    if (!scenario) {
        scenario = make_shared<definition_scenario>();
    }
    scenario->load(esioh);

    if (!isothermal) {
        isothermal = make_shared<support::definition_isothermal>();
    }
    isothermal->load(esioh);

    if (!sg) {
        sg = make_shared<support::definition_largo>();
    }
    sg->load(esioh);

    if (!rad) {
        rad = make_shared<support::definition_radialflow>();
    }
    rad->load(esioh);

    if (!helm) {
        helm = make_shared<support::definition_helm>();
    }
    helm->load(esioh);

    load(esioh, msoln, *scenario, *grid);
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
    // Should we compute fresh statistics or re-use cached values?
#pragma warning(push,disable:1572)
    const bool use_cached = controller && mean && (mean->t == t);
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
    // FIXME Also save boundary layer quantities of interest when one-sided?
    if (!support::save_samples(esioh, *mean)) {
        WARN0(who, "Incomplete statistics saved at time " << t);
    }

    return super::save_statistics_hook(esioh, t);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    if (!mean) {
        mean.reset(new samples());
    }
    support::load_time(esioh, mean->t);

    if (statsdef && statsdef->stale) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Incoming statistics explicitly marked stale; ignoring");
        mean->t = std::numeric_limits<real_t>::quiet_NaN();
    } else if (!support::load_samples(esioh, b, *mean)) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Incomplete statistics loaded from time " << mean->t);
        mean->t = std::numeric_limits<real_t>::quiet_NaN();
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
        if (scenario) {

            if ((boost::math::isnan)(scenario->bulk_rho)) {
                TRACE0(who, "No bulk density scale available so assuming 1");
                velocity = maxmf(abs(scenario->bulk_rho_u) / /*rho*/ 1,
                                 velocity);
            } else {
                velocity = maxmf(abs(scenario->bulk_rho_u) / scenario->bulk_rho,
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

        if (sg->formulation.enabled() && sg->baseflow) {
            // ...prefer wall velocity from inviscid baseflow as freestream
            // surrogate as it is independent of the wall-normal domain size...
            largo_state freestream;
            sg->baseflow->conserved(0.0, freestream.as_is());
            velocity = maxmf(abs(freestream.u()), velocity);
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
                // Four restarts per flow through in a boundary layer giving 40
                // restarts in 10 flow throughs and, across O(25) flow throughs
                // for a simulation campaign, O(100) restart files.
                t /= 4;
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

} // end namespace perfect

} // end namespace suzerain
