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

/** @file
 * Application executing \ref suzerain::perfect::driver_advance::run.
 */

#include <esio/esio.h>
#include <largo/largo.h>
#include <gsl/gsl_spline.h>

#include <suzerain/baseflow.hpp>
#include <suzerain/common.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/constraint_treatment.hpp>
#include <suzerain/format.hpp>
#include <suzerain/largo_state.hpp>
#include <suzerain/radialflow.h>
#include <suzerain/rholut.hpp>
#include <suzerain/specification_isothermal.hpp>
#include <suzerain/specification_zgbsv.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/definition_noise.hpp>
#include <suzerain/support/logging.hpp>

#include "driver.hpp"
#include "hybrid_operator.hpp"
#include "isothermal_mass_operator.hpp"
#include "nonlinear_operator.hpp"
#include "nonreflecting_treatment.hpp"
#include "perfect.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

// TODO Additional QoI that would steer the scenario, e.g. wall blowing.
static
int cev_baseflow_laminar(
    const double dstag,
    double& gammae,
    double& Mae,
    double& pexi,
    double& T_ratio);

/** Application for initializing new restart files. */
struct driver_advance : public driver
{
    driver_advance(const std::string& revstr)
        : driver("Compressible, perfect gas simulation time advancement",
                 "RESTART-FILE",
                 "",
                 revstr)
        , who("advance")
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // namespace perfect

} // namespace suzerain

// Provided by main_advance_svnrev.{c,h} so revstr updates are merely relinking
extern "C" const char revstr[];

/** Instantiate and invoke the application */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_advance app(revstr);
    return app.run(argc, argv);
}

int
suzerain::perfect::driver_advance::run(int argc, char **argv)
{
    using namespace std;
    using boost::math::isnan;

    // Storage for binary-specific options
    noz = make_shared<support::radialflow_definition>(/*deltae*/ 1.0);
    const support::noise_definition noisedef;
    string solver_spec(static_cast<string>(suzerain::zgbsv_specification()));
    string implicit("rhome_xyz");
    real_t cevisslam = numeric_limits<real_t>::quiet_NaN();
    string undriven;

    // Register binary-specific options
    options.add_definition(const_cast<support::noise_definition&>(noisedef));
    options.add_options()
        ("explicit",  boost::program_options::bool_switch(),
                      "Use purely explicit operators")
        ("implicit",  boost::program_options::value(&implicit)
                          ->implicit_value(implicit),
                      "Use hybrid implicit/explicit operators, optionally"
                      " choosing 'rhome_xyz' or 'rhome_y' linearized treatment")
        ("solver",    boost::program_options::value(&solver_spec)
                          ->default_value(solver_spec),
                      "Use the specified algorithm for any --implicit solves")
        ("cevisslam", boost::program_options::value(&cevisslam),
                      "Take homogenized boundary layer baseflow conditions from"
                      " given leeward arc length from the laminar CEV"
                      " international space station return peak heating regime"
                      " stagnation point as measured in meters.")
        ("undriven",  boost::program_options::value(&undriven)
                          ->implicit_value("all"),
                      "Disable all or some physics-related driving forces"
                      " by specifying 'all' or 'rho', 'rho_u', etc.")
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    vector<string> positional = initialize(argc, argv);

    // Select type of timestepping operators to use (default implicit)
    options.conflicting_options("explicit", "implicit");
    options.conflicting_options("explicit", "solver");
    const bool use_explicit =  options.variables()["explicit"].as<bool>();
    const bool use_implicit =  options.variables().count("implicit")
                            || !use_explicit;

    // Validate implicit option; implicit_value "captures" following argument
    // Establish common_block.linearization based on the result
    if (use_implicit) {
        boost::algorithm::trim(implicit);
        if (implicit == "rhome_xyz") {
            common_block.linearization = linearize::rhome_xyz;
        } else if (implicit == "rhome_y") {
            common_block.linearization = linearize::rhome_y;
        } else {
            FATAL0(who, "Unknown --implicit argument:  " << implicit);
            return EXIT_FAILURE;
        }
    } else {
        common_block.linearization = linearize::none;
        implicit = "none";
    }

    if (options.variables().count("cevisslam")) {
        INFO0("Mimicking scenario " << cevisslam
              << " meters leeward of the laminar CEV stagnation point");
        double T_ratio;
        cev_baseflow_laminar(cevisslam,
                             scenario->gamma,  // Notice modification before...
                             scenario->Ma,     // ...call to adjust_scenario().
                             noz->pexi,
                             T_ratio);
        INFO0("Adjusting lower temperature to reproduce CEV T_e / T_w  = "
              << T_ratio);
        isothermal->lower_T = 1 / T_ratio;     // Edge temperature is 1.
    }

    if (positional.size() != 1) {
        FATAL0(who, "Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const string restart_file = positional[0];

    INFO0(who, "Loading restart file: " << restart_file);
    real_t initial_t = numeric_limits<real_t>::quiet_NaN();
    {
        // Preserve exact restart file details via Push/Pop/Merge below
        shared_ptr<scenario_definition> restart_scenario
                = make_shared<scenario_definition>();

        // Load the restart details with state going into state_linear
        shared_ptr<boost::remove_pointer<esio_handle>::type> h( // RAII
                esio_handle_initialize(MPI_COMM_WORLD),
                esio_handle_finalize);
        esio_file_open(h.get(), restart_file.c_str(), 0);
        restart_scenario.swap(scenario);                        // Push
        load_restart(h.get(), initial_t);
        restart_scenario.swap(scenario);                        // Pop
        scenario->populate(*restart_scenario, true);            // Merge

        // Adjust total energy as necessary to account for any scenario change
        state_nonlinear->assign_from(*state_linear);
        adjust_scenario(*state_nonlinear, *scenario, *grid, *dgrid, *cop,
                        restart_scenario->Ma, restart_scenario->gamma);
        state_linear->assign_from(*state_nonlinear);
    }

    if (msoln) {
        INFO0(who, "Restart file prescribes a manufactured solution");
        if (!(isnan)(scenario->bulk_rho)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho = "
                  << scenario->bulk_rho);
        }
        if (!(isnan)(scenario->bulk_rho_u)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_u = "
                  << scenario->bulk_rho_u);
        }
        if (!(isnan)(scenario->bulk_rho_E)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_E = "
                  << scenario->bulk_rho_E);
        }
    }

    // If requested, add noise to the momentum fields at startup (expensive).
    if (noisedef.percent > 0) {
        state_nonlinear->assign_from(*state_linear);
        add_noise(*state_nonlinear, noisedef, *scenario,
                  *grid, *dgrid, *cop, *b);
        state_linear->assign_from(*state_nonlinear);
    }

    // Initialize any requested slow growth forcing workspace
    sg->workspace = NULL; // Defensive
    common_block.slow_treatment = slowgrowth::none;
    if (sg->formulation.enabled()) {
        common_block.slow_treatment = slowgrowth::largo;
        const std::string& model = sg->formulation.name();
        enum  { neq = 5, ns = 0, ntvar = 0 };
        static const char ransmodel[] = "dns";
        INFO0(who, "Allocating Largo model \"" << model << "\" with neq=" << neq
              << ", ns=" << ns << ", ntvar=" << ntvar
              << ", ransmodel=" << ransmodel);
        largo_allocate(&sg->workspace, model.c_str(), neq,
                       ns, ntvar, ransmodel);
        if (!sg->workspace) {
            FATAL0(who, "Largo could not allocate requested model");
            return EXIT_FAILURE;
        }
        if ((isnan)(sg->grdelta)) {
            WARN0(who, "Slow growth rate grdelta is NaN");
        }
    }

    // Prepare any necessary, problem-specific constraints
    shared_ptr<constraint::treatment<operator_common_block> > constrainer(
            new constraint::treatment<operator_common_block>(
                    scenario->Ma, *dgrid, *b, common_block));
    if        (grid->two_sided()) { // Channel per channel_treatment.tex

        INFO0(who, "Establishing driving, channel-like state constraints");
        constrainer->physical[ndx::rho].reset(
                new constraint::reference_bulk(scenario->bulk_rho  , *b));
        constrainer->physical[ndx::mx ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_u, *b));
        constrainer->physical[ndx::e  ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_E, *b));

    } else if (grid->one_sided()) { // Flat plate per plate_treatment.tex
                                    // unless complicated baseflow prescribed

        if (!(isnan)(scenario->bulk_rho)) {
            WARN0(who, "Removing channel-like bulk_rho setting");
            scenario->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        }
        if (!(isnan)(scenario->bulk_rho_u)) {
            WARN0(who, "Removing channel-like bulk_rho_u setting");
            scenario->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
        }
        if (!(isnan)(scenario->bulk_rho_E)) {
            WARN0(who, "Removing channel-like bulk_rho_E setting");
            scenario->bulk_rho_E = numeric_limits<real_t>::quiet_NaN();
        }

        DEBUG0(who, "Computing freestream per vanilla flat plate scenario");
        largo_state freestream;
        {
            real_t T_inf   = isothermal->upper_T;                 // Brevity
            real_t mx_inf  = pow(T_inf, scenario->beta);          // Eqn ( 7)
            real_t rho_inf = pow(T_inf, scenario->beta - 0.5);    // Eqn ( 8)
            real_t e_inf   = pow(T_inf, scenario->beta - 0.5) * ( // Eqn (11)
                                 T_inf/(scenario->gamma*(scenario->gamma-1))
                               + pow(scenario->Ma, 2) / 2 * (
                                     T_inf
                                   + pow(isothermal->upper_v, 2)
                                   + pow(isothermal->upper_w, 2)
                                 )
                             );
            real_t p_inf;
            rholut::p(scenario->gamma, rho_inf, T_inf, p_inf);

            // Pack the into Largo-ready storage choosing lower_v for v
            // This matters for establishing isothermal->upper_v below
            freestream.rho = rho_inf;
            freestream.mx  = mx_inf;
            freestream.my  = rho_inf * isothermal->lower_v;
            freestream.mz  = rho_inf * isothermal->upper_w;
            freestream.e   = e_inf;
            freestream.p   = p_inf;
        }

        // If the baseflow is specified by a radial nozzle problem...
        if (sg->formulation.enabled() && !noz->trivial()) {

            INFO0(who, "Finding baseflow per suzerain_radialflow_qoi_match...");
            // ...unpack all information from *noz and *scenario...
            const double deltae = noz->deltae;
            const double gamma  = !(isnan)(noz->gamma) && (noz->gamma != 0)
                                ? noz->gamma : scenario->gamma;
            const double Mae    = !(isnan)(noz->Mae ) && (noz->Mae  != 0)
                                ? noz->Mae  : scenario->Ma;
            const double pexi   = noz->pexi;

#           define POSSIBLYWARN0(who, pre, val)                           \
                if ((isnan)(val)) { WARN0(who, pre << fullprec<>(val)); } \
                else              { INFO0(who, pre << fullprec<>(val)); }

            POSSIBLYWARN0(who, "Baseflow requested deltae:     ", deltae);
            POSSIBLYWARN0(who, "Baseflow requested gamma:      ", gamma );
            POSSIBLYWARN0(who, "Baseflow requested Mae:        ", Mae   );
            POSSIBLYWARN0(who, "Baseflow requested pexi:       ", pexi  );

            // ...compute edge state matching requested targets...
            double Ma0, R[2], uR, rhoR, pR;
            suzerain_radialflow_qoi_match(
                    deltae, gamma, Mae, pexi,
                    &Ma0, R+1, R+0, &uR, &rhoR, &pR);

            POSSIBLYWARN0(who, "Matching radial flow has Ma0:  ", Ma0 );
            POSSIBLYWARN0(who, "Matching radial flow has R:    ", R[0]);
            POSSIBLYWARN0(who, "Matching radial flow has uR:   ", uR  );
            POSSIBLYWARN0(who, "Matching radial flow has rhoR: ", rhoR);
            POSSIBLYWARN0(who, "Matching radial flow has pR:   ", pR  );
            POSSIBLYWARN0(who, "Matching radial flow has R0:   ", R[1]);

            // ...solve for wall state given edge state...
            shared_ptr<suzerain_radialflow_solution> r(
                    suzerain_radialflow_solver(
                        Ma0, gamma, uR, rhoR, pR, R, sizeof(R)/sizeof(R[0])),
                    free);

            const double u1   = r->state[1].u;
            const double rho1 = r->state[1].rho;
            const double p1   = r->state[1].p;
            POSSIBLYWARN0(who, "Matching radial flow has u1:   ", u1  );
            POSSIBLYWARN0(who, "Matching radial flow has rho1: ", rho1);
            POSSIBLYWARN0(who, "Matching radial flow has p1:   ", p1  );

#           undef POSSIBLYWARN0

            // ...solve problem at radii fixed by B-spline collocation points
            INFO0(who, "Preparing baseflow with suzerain_radialflow_solver");
            ArrayXr y(b->n());
            for (int i = 0; i < y.size(); ++i) {
                y[i] = b->collocation_point(i);
            }
            ArrayXr S = (y.abs2() + pow(r->state[1].R, 2)).sqrt();
            shared_ptr<suzerain_radialflow_solution> s(
                    suzerain_radialflow_solver(
                        Ma0, gamma, u1, rho1, p1, S.data(), S.size()),
                    free);

            // ...and tuck solution into a baseflow_map as a function of y
            shared_ptr<baseflow_map> bm(new baseflow_map());
            for (int i = 0; i < y.size(); ++i) {
                baseflow_map::row& row = bm->table[y[i]];
                suzerain_radialflow_cartesian_conserved(
                        s.get(), i, scenario->Ma,
                        &row.  base.rho,
                        &row.  base.mx ,
                        &row.  base.my ,
                        &row.  base.e  ,
                        &row.  base.p  ,
                        &row.dxbase.rho,
                        &row.dxbase.mx ,
                        &row.dxbase.my ,
                        &row.dxbase.e  ,
                        &row.dxbase.p  ,
                        &row.dybase.rho,
                        &row.dybase.mx ,
                        &row.dybase.my ,
                        &row.dybase.e  ,
                        &row.dybase.p);
                assert(row.  base.mz == 0);
                assert(row.dxbase.mz == 0);
                assert(row.dybase.mz == 0);
            }
            sg->baseflow = bm;
        }

        // If required by the slow growth model, add freestream baseflow
        if (     sg->formulation.enabled()
             && !sg->formulation.is_strictly_temporal()
             && (   !sg->baseflow
                 || dynamic_cast<baseflow_uniform*>(sg->baseflow.get()))) {
            INFO0(who, "Adding uniform baseflow at freestream conditions");
            shared_ptr<baseflow_uniform> bu(new baseflow_uniform());
            bu->x.resize(state_linear->shape()[0] + /* pressure */ 1);
            memcpy(bu->x.data(), freestream.as_is(), sizeof(freestream));
            sg->baseflow = bu;
        }

        // If slow growth baseflow in use, obtain freestream from it.  Permits
        // obtaining freestream from, e.g., baseflow_polynomial in a manner
        // consistent with slow growth forcing.  Notice when no slow growth we
        // still have vanilla flat plate scenario values after this block.
        if (sg->formulation.enabled() && sg->baseflow) {
            {
                largo_state dontcare;
                sg->baseflow->conserved(grid->L.y(), freestream.as_is(),
                                        dontcare.as_is(), dontcare.as_is());
            }
            {
                real_t dontcare;
                sg->baseflow->pressure (grid->L.y(), freestream.p,
                                        dontcare, dontcare);
            }
        }

        INFO0(who, "Setting freestream reference state on upper boundary");
        isothermal->upper_rho = freestream.rho;
        isothermal->upper_u   = freestream.u();
        if (sg->formulation.enabled()) {
            // Homogenization modifies the inviscid characteristics.
            // See discussion (or placeholder thereof) at Redmine #2982
            isothermal->upper_v = freestream.v()
                                - grid->L.y()*sg->grdelta;
            INFO0(who, "Matched upper_v with prescribed"
                       " grdelta, freestream.v(), and Ly: "
                       << isothermal->upper_v);

            isothermal->upper_w = freestream.w();
            INFO0(who, "Matched upper_w with prescribed freestream.w(): "
                       << isothermal->upper_w);

            // Employs perfect gas EOS with constant gamma.
            isothermal->upper_T = scenario->gamma*freestream.p/freestream.rho;
            INFO0(who, "Matched upper_T with prescribed freestream T: "
                       << isothermal->upper_T);
        } else {
            DEBUG0(who, "Keeping reference upper_v = " << isothermal->upper_v);
            DEBUG0(who, "Keeping reference upper_w = " << isothermal->upper_w);
            DEBUG0(who, "Keeping reference upper_T = " << isothermal->upper_T);
        }

        INFO0(who, "Establishing driving, freestream-like state constraints");
        constrainer->physical[ndx::e  ].reset(
                new constraint::constant_upper(freestream.e,   *cop, 0));
        constrainer->physical[ndx::mx ].reset(
                new constraint::constant_upper(freestream.mx,  *cop, 0));
        constrainer->physical[ndx::rho].reset(
                new constraint::constant_upper(freestream.rho, *cop, 0));
        if (sg->formulation.enabled() && sg->baseflow) {
            INFO0(who, "Matching freestream constraint enforcement"
                       " profile to baseflow");
            largo_state state, dontcare;
            for (int i = 0; i < b->n(); ++i) {
                sg->baseflow->conserved(b->collocation_point(i),
                                        state.as_is(),
                                        dontcare.as_is(),
                                        dontcare.as_is());
                constrainer->physical[ndx::e  ]->shape[i] = state.e  ;
                constrainer->physical[ndx::mx ]->shape[i] = state.mx ;
                constrainer->physical[ndx::rho]->shape[i] = state.rho;
            }
            for (size_t i = 0; i < constrainer->physical.size(); ++i) {
                if (!constrainer->physical[i]->enabled()) {
                    continue;
                }
                real_t minval = constrainer->physical[i]->shape.minCoeff();
                real_t maxval = constrainer->physical[i]->shape.maxCoeff();
                INFO0(who, "Constraint profile for " << ndx::identifier[i]
                            << " has range ["
                            << fullprec<>(minval) << ", "
                            << fullprec<>(maxval) << ']');
            }
        }

        if (isothermal->upper_v > 0 && (boost::math::isinf)(scenario->Re)) {
            WARN0(who, "Nonreflecting viscous outflow boundary problematic"
                       << " (Redmine #2983)");
        }

    } else {

        FATAL0(who, "Sanity error in constraint selection");
        return EXIT_FAILURE;

    }

    // Prepare spatial operators depending on requested advancement type.
    // Notice constrainer always wraps the implicit operator
    // and that the same nonlinear_operator is used pervasively.
    L = constrainer;
    N.reset(new nonlinear_operator(
                *scenario, *grid, *dgrid, *cop, *b, common_block, *sg, msoln));
    if (use_explicit) {

        INFO0(who, "Initializing fully explicit spatial operators");
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<nonreflecting_treatment> nonreflecting(
                    new nonreflecting_treatment(
                        *scenario, *isothermal, *grid, *dgrid, *cop, *b));
            nonreflecting->N = N;
            N = nonreflecting;
        }
        constrainer->L.reset(new isothermal_mass_operator(
                    *scenario, *isothermal,
                    *grid, *dgrid, *cop, *b, common_block));

    } else if (use_implicit) {

        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");
        INFO0(who, "Implicit linearization employed: " << implicit);
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            FATAL0(who, "Nonreflecting upper boundary treatment"
                        " not currently usable with implicit advance");
            return EXIT_FAILURE;
        }
        constrainer->L.reset(new isothermal_hybrid_linear_operator(
                    solver_spec, *scenario, *isothermal,
                    *grid, *dgrid, *cop, *b, common_block));

    } else {

        FATAL0(who, "Sanity error in operator selection");
        return EXIT_FAILURE;

    }

    // Use --undriven as a testing- and debugging-related tool.
    // For example, to investigate nonreflecting boundary condition behavior.
    if (options.variables().count("undriven")) {
        boost::algorithm::trim(undriven);
        const size_t j = distance(
                ndx::identifier.begin(),
                find(ndx::identifier.begin(), ndx::identifier.end(), undriven));
        if (j < ndx::identifier.static_size) {
            INFO0(who, "Disabling any " << ndx::description[j]
                       << " physical constraint per --undriven=" << undriven);
            constrainer->physical[j] = constrainer->none;
        } else if (undriven == "all") {
            INFO0(who, "Disabling all physical constraints per --undriven="
                       << undriven);
            fill(constrainer->physical.begin(), constrainer->physical.end(),
                 constrainer->none);
        } else {
            FATAL0("Unknown --undriven argument:  " << undriven);
            return EXIT_FAILURE;
        }
    }

    // Perform final housekeeping and then advance time as requested
    establish_ieee_mode();
    log_discretization_quality();
    prepare_controller(initial_t, dgrid->chi());
    save_metadata();
    const real_t elapsed_wall_time // Negative on error
        = advance_controller(timedef->status_final,
                             statsdef->final,
                             restartdef->final);

    // If we advanced by any time steps, log the observed linearization error
    if (   !use_explicit
        && elapsed_wall_time >= 0
        && controller->current_nt() > 0) {
        log_linearization_error(build_timeprefix(controller->current_t(),
                                                 controller->current_nt()),
                                controller->current_t(),
                                controller->current_nt());
    }

    // Deallocate any slow growth forcing workspace
    if (sg->formulation.enabled()) {
        if (sg->workspace) {
            largo_deallocate(&sg->workspace);
        }
        sg->workspace = NULL; // Defensive
    }

    // Report error to the OS iff advance_control reported an error
    return elapsed_wall_time < 0;
}

static
int suzerain::perfect::cev_baseflow_laminar(
    const double dstag,
    double& gammae,
    double& Mae,
    double& pexi,
    double& T_ratio)
{
    // Data extracted from notebooks/cev_laminar.in which was originally is based upon
    // https://svn.ices.utexas.edu/repos/pecos/turbulence/heatshield_bl/trunk/laminar/scenario.dat
    enum { N = 71 };
    static const double data_dstag[N] = {
        -0.6143492210355781,    -0.59390333664534101,  -0.57179300778174724,
        -0.54789092718131283,   -0.52206122007465172,  -0.49415938719521613,
        -0.46403232768676173,   -0.43151718546930518,  -0.3964421451945217,
        -0.35862567789533895,   -0.31787735472377765,  -0.27399655136109757,
        -0.2267747407134193,    -0.17599395375299576,  -0.12142773350753355,
        -0.062843090903680121,  -0,                    0.067346704177323247,
        0.13944580580883326,    0.21654744952740934,   0.29890258933871827,
        0.38675857037813088,    0.48035747367467119,   0.57993312475275527,
        0.68103062660465374,    0.78212886925651937,   0.88322732787220115,
        0.9843255166409417,     1.0854233790807597,    1.186520826236158,
        1.2876186095690088,     1.3887167424250881,    1.4898149297790302,
        1.5909130106137055,     1.6920108449633222,    1.7931083547689379,
        1.8942058645744817,     1.9953036989240984,    2.0964017797587737,
        2.1974999671127158,     2.2985980999687952,    2.399695883301646,
        2.5007933304570438,     2.6018911928968622,    2.7029893816656028,
        2.8040878402812845,     2.9051860829331497,    3.0062835847850486,
        3.1058592358631327,     3.199458139159673,     3.2873141201990856,
        3.3696692600103946,     3.4467709037289707,    3.5188700053604807,
        3.5862167095378044,     3.6490598004414849,    3.7076444430453384,
        3.7622106632908006,     3.8129914502512245,    3.8602132608989028,
        3.9040940642615829,     3.9448423874331442,    3.9826588547323265,
        4.01773389500711,       4.0502490372245665,    4.0803760967330209,
        4.1082779296124556,     4.1341076367191167,    4.1580097173195512,
        4.1801200461831449,     4.2005659305733829
    };
    static const double data_gammae[N] = {
        1.4099815429390452,  1.4098104536887972,  1.4096908057455511,
        1.4096108296401819,  1.4095556146292538,  1.4095189225317255,
        1.4094954447208019,  1.4094812278433841,  1.4094741105376656,
        1.4094728832908761,  1.4094750964867506,  1.409479332470402,
        1.4094845509118477,  1.4094896108287147,  1.409493739436265,
        1.409496087805638,   1.4094965858307458,  1.4094949070017295,
        1.4094907750622041,  1.4094836020905861,  1.4094729365300411,
        1.4094581868477518,  1.409438670326961,   1.4094139594507757,
        1.4093838947568316,  1.409349694990228,   1.4093108365063949,
        1.4092681785505701,  1.4092214742059146,  1.4091712712490843,
        1.4091178420751405,  1.4090617036849444,  1.4090029840768459,
        1.4089424630791352,  1.4088808377294937,  1.4088184299707727,
        1.4087553064055358,  1.4086914743946093,  1.4086264139081992,
        1.4085601149390066,  1.4084929150788372,  1.4084256297807161,
        1.4083595370046604,  1.4082956118387255,  1.4082345113997261,
        1.4081769445147463,  1.4081242617376561,  1.4080790926320617,
        1.4080536636759231,  1.4080294579155366,  1.4080078090065777,
        1.4079951744388131,  1.407991112426707,   1.4079948516205061,
        1.4080058983225989,  1.4080234319050933,  1.4080461078617779,
        1.4080733203251381,  1.4081050590057502,  1.4081403498844409,
        1.4081787684846814,  1.4082204150965079,  1.408265818053513,
        1.4083167289094889,  1.4083781550996435,  1.4084445755187867,
        1.4085010844389505,  1.4085542572945138,  1.4086045117223955,
        1.4086554301854781,  1.4087083877256554
    };
    static const double data_Mae[N] = {
        0.42190694027972758,   0.38261248638894363,    0.3477875752073547,
        0.31640737083790416,   0.2876995764882363,     0.25985405007086548,
        0.23489666492710973,   0.21197219841991385,    0.18901267536866218,
        0.16560121027857475,   0.14356407281720557,    0.12173776923151167,
        0.099465572885036532,  0.076881172916196644,   0.054136052992762411,
        0.030073499145841599,  0.0054429205954211607,  0.020625533435922171,
        0.04710689731459318,   0.074540336218711409,   0.10229648360013355,
        0.13005033212091255,   0.15869200176946574,    0.18735787088001371,
        0.21665239131864536,   0.24514188073129567,    0.27350049613561794,
        0.30155122032277992,   0.32914547794748539,    0.35667346504868797,
        0.38391132699969333,   0.41120414466739352,    0.43862203023951135,
        0.46611985089770513,   0.49373079729648561,    0.52144481034749257,
        0.54927217716869747,   0.57704944407135428,    0.60480443623717972,
        0.63230228827088442,   0.65970364956076455,    0.68695701225763917,
        0.71394455827465408,   0.7409908147908203,     0.76817051568966199,
        0.79517940351262706,   0.82219640239697489,    0.84900156497529833,
        0.87490773103516306,   0.89858115487825818,    0.92143086386963002,
        0.94344288692974954,   0.9637540241995538,     0.98245494114359588,
        1.000847126522,        1.0177358093237081,     1.0332195780131908,
        1.0481561893730225,    1.0623982037546464,     1.075287642592152,
        1.0874605127903292,    1.0985922749827015,     1.1093731426053473,
        1.118775131627562,     1.1272786975850249,     1.1348406174944963,
        1.1433080447195736,    1.1521944684829317,     1.1621859200981024,
        1.1744408936560817,    1.1905934935259677
    };
    static const double data_pexi[N] = {
        -0.10101086946499484,    -0.095023446053762964,   -0.089182343095845121,
        -0.083447045126352787,   -0.079212639987524411,   -0.077883638614280201,
        -0.076400383502228372,   -0.075017607543867923,   -0.076596710412302649,
        -0.081926940030949899,   -0.087216018089626909,   -0.095805208115930338,
        -0.11274893621616011,    -0.14376217440601388,    -0.2019320381281717,
        -0.39128487154052993,    -3.5106047841435903,     -0.32211414005226363,
        -0.17260017496292279,    -0.10940150511034639,    -0.077877994719461052,
        -0.059829843422088434,   -0.047537526508757005,   -0.038927185529409138,
        -0.033207371957848009,   -0.029286208180788198,   -0.026241315356759633,
        -0.023850353644973017,   -0.021862794866939709,   -0.0202783088695449,
        -0.018972065957037879,   -0.017931635525217017,   -0.017105486287050446,
        -0.016420139593348454,   -0.01582059702106817,    -0.01527493730587331,
        -0.014755313619465062,   -0.014217427339774808,   -0.013690304564595569,
        -0.013159861199299431,   -0.012694892400477699,   -0.012285393056196152,
        -0.011892537649319856,   -0.011571158658678039,   -0.011314885044750676,
        -0.01104037110704819,    -0.010840033105601166,   -0.010830296722694014,
        -0.010573127574118263,   -0.010191783208898577,   -0.010111390042733736,
        -0.010087225149769299,   -0.010004726557591544,   -0.0098711124687766924,
        -0.0098858336057264413,  -0.0098399105377029502,  -0.009736713144315564,
        -0.0097314199000317486,  -0.009786887911991721,   -0.0097198961918436976,
        -0.0096911626691443126,  -0.0096502759056394972,  -0.0097415965750651881,
        -0.0097736338699843522,  -0.010045362685481096,   -0.010248414045932849,
        -0.011066019395791476,   -0.012350269196758374,   -0.014683079037586257,
        -0.018803670726498776,   -0.025438605531063266
    };
    static const double data_T_ratio[N] = {
        3.4872999452040814,  3.6479923945555432,  3.7199291691044016,
        3.7644092644152076,  3.7999017268666924,  3.8313715997967619,
        3.8602181515453204,  3.8857189776984429,  3.9097235134249479,
        3.9313475230762478,  3.9505823320748172,  3.9669976109161964,
        3.9807965955454621,  3.9918326739643462,  4.0007546597793091,
        4.0081656310024787,  4.015032259169498,   4.0221316905250157,
        4.0308078701497205,  4.0407977228460457,  4.0515352658071633,
        4.0623404117234481,  4.0725715903486446,  4.0830141612741606,
        4.0933729464407564,  4.1010991108409218,  4.1072300067106804,
        4.1125933915687005,  4.1174438736460042,  4.1217382515014434,
        4.1254436508821914,  4.1291175669260829,  4.1334874468323957,
        4.1379994390070687,  4.142982497780304,   4.1482700949021858,
        4.1541354238844006,  4.1611072377993317,  4.1677384625970166,
        4.1748141000911261,  4.181692718061174,   4.1885944038008942,
        4.1960728537236802,  4.2035455920447395,  4.2112397768750833,
        4.2197680032275438,  4.2295387519018597,  4.2400864198458468,
        4.2510425457759693,  4.2619128750651454,  4.2714687975204821,
        4.2805456956631689,  4.2880414432341372,  4.2952355514753942,
        4.3010590003939635,  4.3057514460810973,  4.3102253879254713,
        4.3134371768417648,  4.3161547307280461,  4.3181456381433216,
        4.3196888662022648,  4.3203828218121583,  4.3205184682086744,
        4.3194370500077701,  4.3168567736825629,  4.3119612190476326,
        4.302792717238483,   4.2853609220662268,  4.2545325614550382,
        4.1864251808993052,  4.0039925486698449
    };

    if (dstag < data_dstag[0]) {
        std::ostringstream oss;
        oss << "Requested dstag " << dstag
            << " is less than minimum " << data_dstag[0];
        throw std::domain_error(oss.str());
    }

    if (dstag > data_dstag[N-1]) {
        std::ostringstream oss;
        oss << "Requested dstag " << dstag
            << " is greater than maximum " << data_dstag[N-1];
        throw std::domain_error(oss.str());
    }

    shared_ptr<gsl_interp_accel> a(gsl_interp_accel_alloc(),
                                   gsl_interp_accel_free);
    shared_ptr<gsl_interp>       f(gsl_interp_alloc(gsl_interp_cspline, N),
                                   gsl_interp_free);

    gsl_interp_init(f.get(), data_dstag, data_gammae, N);
    gammae = gsl_interp_eval(f.get(), data_dstag, data_gammae, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_Mae, N);
    Mae = gsl_interp_eval(f.get(), data_dstag, data_Mae, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_pexi, N);
    pexi = gsl_interp_eval(f.get(), data_dstag, data_pexi, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_T_ratio, N);
    T_ratio = gsl_interp_eval(f.get(), data_dstag, data_T_ratio, dstag, a.get());

    return 0;
}
