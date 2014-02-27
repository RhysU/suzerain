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

#include <suzerain/common.hpp>
#include <suzerain/baseflow.hpp>
#include <suzerain/cev.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/treatment_constraint.hpp>
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
#include "operator_hybrid_isothermal.hpp"
#include "operator_mass_isothermal.hpp"
#include "operator_nonlinear.hpp"
#include "treatment_nonreflecting.hpp"
#include "perfect.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

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
    using boost::math::isfinite;

    // Storage for binary-specific options
    const support::definition_noise noisedef;
    string solver_spec(static_cast<string>(suzerain::specification_zgbsv()));
    string implicit("rhome_xyz");
    real_t cevisslam = numeric_limits<real_t>::quiet_NaN();
    string undriven;

    // Register binary-specific options
    options.add_definition(const_cast<support::definition_noise&>(noisedef));
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
        cev::iss_laminar(cevisslam,
                         scenario->gamma,  // Notice modification before...
                         scenario->Ma,     // ...call to adjust_scenario().
                         rad->pexi,
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
        shared_ptr<definition_scenario> restart_scenario
                = make_shared<definition_scenario>();

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
    if        (grid->two_sided()) { // Channel per treatment_channel.tex

        INFO0(who, "Establishing driving, channel-like state constraints");
        constrainer->physical[ndx::rho].reset(
                new constraint::reference_bulk(scenario->bulk_rho  , *b));
        constrainer->physical[ndx::mx ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_u, *b));
        constrainer->physical[ndx::e  ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_E, *b));

    } else if (grid->one_sided()) { // Flat plate per treatment_plate.tex
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

        // If the baseflow is specified by a radial flow problem...
        if (sg->formulation.enabled() && !rad->trivial()) {

            INFO0(who, "Finding baseflow per suzerain_radialflow_qoi_match...");
            // ...unpack all information from *rad and *scenario...
            const double deltae = rad->deltae;
            const double gamma  = !(isnan)(rad->gamma) && (rad->gamma != 0)
                                ? rad->gamma : scenario->gamma;
            const double Mae    = !(isnan)(rad->Mae ) && (rad->Mae  != 0)
                                ? rad->Mae  : scenario->Ma;
            const double pexi   = rad->pexi;

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

        // If a slow growth model is in use, compute defect growth rates.
        // The gramp_{mean,rms} vectors should already have been zero-filled.
        // All nondimensional expressions below in model document section 5.4.
        // Beware v_w used in wall-normal direction but used in energy result.
        if (sg->formulation.enabled()) {

            if (sg->gramp_mean.size()) {
                INFO0(who, "Ignoring incoming, non-normative defect"
                           " mean slow growth rates");
            }
            sg->gramp_mean.assign(fields.size() + /*pressure*/1, 0.0);

            INFO0(who, "Preparing to compute defect slow growth"
                       " rates from wall state");
            const real_t Tw = isothermal->lower_T;      // Brevity
            const real_t Ew = Tw / (scenario->gamma*(scenario->gamma - 1))
                            + pow(scenario->Ma * isothermal->lower_v, 2) / 2;
            largo_state iw, dxiw;
            if (sg->baseflow) {
                largo_state dontcare;
                sg->baseflow->conserved(
                        0.0, iw.as_is(), dontcare.as_is(), dxiw.as_is());
                sg->baseflow->pressure (
                        0.0, iw.p, dontcare.p, dxiw.p); // as_is()
            }

            // Compute density, model-specific rates, and last pressure
            largo_state grDA;
            grDA.rho = (Tw * dxiw.rho - scenario->gamma * dxiw.p)
                     / (Tw *   iw.rho - scenario->gamma *   iw.p);
            if (sg->formulation.expects_conserved_growth_rates()) {
                INFO0(who, "Calculating conserved defect mean growth rates");
                grDA.mx = dxiw.mx / iw.mx;
                grDA.my = dxiw.my / iw.my;
                grDA.mz = dxiw.mz / iw.mz;
                grDA.e  = (Tw * dxiw.e  - scenario->gamma * Ew * dxiw.p)
                        / (Tw *   iw.e  - scenario->gamma * Ew *   iw.p);
            } else if (sg->formulation.expects_specific_growth_rates()) {
                INFO0(who, "Calculating primitive defect mean growth rates");
                // Primitive order matches conserved, hence weird assignments
                // First lines compute derivative from conserved base flow
                grDA.mx = (dxiw.mx - iw.u()*dxiw.rho) / iw.rho / (iw.u()     );
                grDA.my = (dxiw.my - iw.v()*dxiw.rho) / iw.rho / (iw.v()     );
                grDA.mz = (dxiw.mz - iw.w()*dxiw.rho) / iw.rho / (iw.w()     );
                grDA.e  = (dxiw.e  - iw.E()*dxiw.rho) / iw.rho / (iw.E() - Ew);
            } else {
                FATAL0(who, "Sanity error in growth rate computations");
                return EXIT_FAILURE;
            }
            grDA.p = 0;

            // Coerce any non-finite rates to be zero per goofy convention
            for (size_t i = 0; i < sg->gramp_mean.size(); ++i) {
                if (!(isfinite)(grDA.as_is()[i])) {
                    DEBUG0(who, "Treating non-finite mean growth rate index "
                               << i << " as zero");
                    grDA.as_is()[i] = 0;
                }
            }

            // Report any non-zero growth rates to the user
#           define MAYBE_MENTION(var, name)                   \
                do if (var != 0) {                            \
                    INFO0(who, "Mean defect growth rate for " \
                               << name << fullprec<>(var));   \
                } while (0)
            MAYBE_MENTION(grDA.rho, "density:              ");
            if (sg->formulation.expects_conserved_growth_rates()) {
                MAYBE_MENTION(grDA.mx,  "streamwise momentum:  ");
                MAYBE_MENTION(grDA.my,  "wall-normal momentum: ");
                MAYBE_MENTION(grDA.mz,  "spanwise momentum:    ");
                MAYBE_MENTION(grDA.e,   "total energy:         ");
            } else if (sg->formulation.expects_specific_growth_rates()) {
                MAYBE_MENTION(grDA.mx,  "streamwise velocity:  ");
                MAYBE_MENTION(grDA.my,  "wall-normal velocity: ");
                MAYBE_MENTION(grDA.mz,  "spanwise velocity:    ");
                MAYBE_MENTION(grDA.e,   "specific energy:      ");
            } else {
                FATAL0(who, "Sanity error in growth rate reporting");
                return EXIT_FAILURE;
            }
            MAYBE_MENTION(grDA.p,   "pressure:             ");
#undef      MAYBE_MENTION

            // Finally set the rates into the slow growth specification
            for (size_t i = 0; i < sg->gramp_mean.size(); ++i) {
                sg->gramp_mean[i] = grDA.as_is()[i];
            }

            if (sg->gramp_rms.size()) {
                INFO0(who, "Ignoring incoming, non-normative defect"
                           " root-mean-square slow growth rates");
            }
            INFO0(who, "Disabling root-mean-square defect slow growth");
            sg->gramp_rms.assign(sg->gramp_mean.size(), 0.0);
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
                INFO0(who, "Constraint profile for "
                            << setw(5) << left << ndx::identifier[i]
                            << " has range ["
                            << fullprec<>(minval) << ", "
                            << fullprec<>(maxval) << ']');
            }
        }

        if (isothermal->upper_v > 0 && !(boost::math::isinf)(scenario->Re)) {
            WARN0(who, "Nonreflecting viscous outflow boundary problematic"
                       << " (Redmine #2983)");
        }

    } else {

        FATAL0(who, "Sanity error in constraint selection");
        return EXIT_FAILURE;

    }

    // Prepare spatial operators depending on requested advancement type.
    // Notice constrainer always wraps the implicit operator
    // and that the same operator_nonlinear is used pervasively.
    L = constrainer;
    N.reset(new operator_nonlinear(
                *scenario, *grid, *dgrid, *cop, *b, common_block, *sg, msoln));
    if (use_explicit) {

        INFO0(who, "Initializing fully explicit spatial operators");
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<treatment_nonreflecting> nonreflecting(
                    new treatment_nonreflecting(
                        *scenario, *isothermal, *grid, *dgrid, *cop, *b));
            nonreflecting->N = N;
            N = nonreflecting;
        }
        constrainer->L.reset(new operator_mass_isothermal(
                    *scenario, *isothermal,
                    *grid, *dgrid, *cop, *b, common_block));

    } else if (use_implicit) {

        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");
        INFO0(who, "Implicit linearization employed: " << implicit);
        shared_ptr<operator_hybrid_isothermal> hybrid(
                    new operator_hybrid_isothermal(
                        solver_spec, *scenario, *isothermal,
                        *grid, *dgrid, *cop, *b, common_block));
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            hybrid->N = N;
            N = hybrid;
        }
        constrainer->L = hybrid;  // Constrainer invokes hybrid linear operator

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
        } else if (undriven == "none") { // NOP, specifies the inverse of "all"
            INFO0(who, "Enabling all physical constraints per --undriven="
                       << undriven);
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
