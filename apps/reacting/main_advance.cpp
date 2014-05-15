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
#include <suzerain/constraint.hpp>
#include <suzerain/treatment_constraint.hpp>
#include <suzerain/error.h>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_hybrid_residual.hpp>
#include <suzerain/specification_zgbsv.hpp>
#include <suzerain/support/definition_noise.hpp>
#include <suzerain/support/logging.hpp>

#include "driver.hpp"
#include "operator_explicit.hpp"

#include "operator_hybrid.hpp"
#include "treatment_nonreflecting.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace reacting {

/** Application for initializing new restart files. */
struct driver_advance : public driver
{
    driver_advance()
        : driver("Compressible, reacting flow simulation time advancement",
                 "RESTART-FILE",
                 "",
                 REVISIONSTR)
        , who("advance")
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // namespace reacting

} // namespace suzerain

/** Instantiate and invoke the application */
int main(int argc, char **argv)
{
    suzerain::reacting::driver_advance app;
    return app.run(argc, argv);
}

int
suzerain::reacting::driver_advance::run(int argc, char **argv)
{
    using namespace std;
    using boost::math::isnan;

    // Storage for binary-specific options
    const support::definition_noise noisedef;
    string solver_spec(static_cast<string>(suzerain::specification_zgbsv()));
    string filter_spec("viscous");
    string undriven;

    // Register binary-specific options
    options.add_definition(const_cast<support::definition_noise&>(noisedef));
    options.add_options()
        ("explicit", boost::program_options::bool_switch(),
                     "Use purely explicit operators")
        ("implicit", boost::program_options::bool_switch(),
                     "Use hybrid implicit/explicit operators")
        ("solver",   boost::program_options::value(&solver_spec)
                         ->default_value(solver_spec),
                     "Use the specified algorithm for any implicit solves")
        ("filter",   boost::program_options::value(&filter_spec)
                         ->default_value(filter_spec),
                     "Use the specified type to construct filter source")
        ("undriven", boost::program_options::value(&undriven)
                         ->implicit_value("all"),
                     "Disable all or some physics-related driving forces"
                     " by specifying 'all' or 'rho', 'rho_u', etc.")
        ("preserve_upper_input", boost::program_options::bool_switch(),
                     "Preserve upper input or file values (applies to"
                     " slow growth baseflow scenarios)")
        ("preserve_largo_gramp_input", boost::program_options::bool_switch(),
                     "Preserve gramp file values (applies to"
                     " slow growth bl_spatiotemporal_consistent model)")
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    vector<string> positional = initialize(argc, argv);

    // Select type of timestepping operators to use (default implicit)
    options.conflicting_options("explicit", "implicit");
    const bool use_explicit =  options.variables()["explicit"].as<bool>();
    const bool use_implicit =  options.variables()["implicit"].as<bool>()
                            || !use_explicit;

    // Only one implicit option now, so just use if we got '--implicit'
    if (use_implicit) {
        common_block.linearization = linearize::rhome_y;
    } else {
        common_block.linearization = linearize::none;
    }

    // Select type of filtering to use (default viscous)
    const bool use_filter =  options.variables().count("filter");

    if (use_filter) {
        boost::algorithm::trim(filter_spec);
        if (filter_spec == "none") {
            common_block.filter_treatment = filter::none;
        } else if (filter_spec == "cook") {
            common_block.filter_treatment = filter::cook;
        } else if (filter_spec == "viscous") {
            common_block.filter_treatment = filter::viscous;
        } else {
            FATAL0("Unknown --filter argument:  " << filter_spec);
            return EXIT_FAILURE;
        }
    } else {
        common_block.filter_treatment = filter::viscous;
    }

    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const string restart_file = positional[0];

    // Reset or preserve the file or input upper values
    const bool preserve_upper_input =
        options.variables()["preserve_upper_input"].as<bool>();

    // Reset or preserve the file amplitude growth rate values
    const bool preserve_largo_gramp_input =
        options.variables()["preserve_largo_gramp_input"].as<bool>();

    INFO0(who, "Loading restart file: " << restart_file);
    real_t initial_t = numeric_limits<real_t>::quiet_NaN();
    {
        // Preserve exact restart file details via Push/Pop/Merge below
        shared_ptr<definition_channel> restart_chdef
                = make_shared<definition_channel>();

        // Load the restart details with state going into state_linear
        shared_ptr<boost::remove_pointer<esio_handle>::type> h( // RAII
                esio_handle_initialize(MPI_COMM_WORLD),
                esio_handle_finalize);
        esio_file_open(h.get(), restart_file.c_str(), 0);
        restart_chdef.swap(chdef);
        load_restart(h.get(), initial_t);
        restart_chdef.swap(chdef);
        chdef->populate(*restart_chdef, true);

        // FIXME: filter source
        //        initialize with info from restart here (?)

        // Adjust total energy as necessary to account for any scenario change
        state_nonlinear->assign_from(*state_linear);
        // FIXME: Currently no functionality to adjust scenario
        //adjust_scenario(*state_nonlinear, *scenario, *grid, *dgrid, *cop, restart_scenario->Ma, restart_scenario->gamma);
        state_linear->assign_from(*state_nonlinear);
    }

    INFO0(who, "Initializing antioch_constitutive");
    cmods->init_antioch();

    common_block.Ns = cmods->Ns();

    // However, if msoln was provided, match its contents to other members
    // Do here b/c cmods doesn't make sense until after init_antioch
    if (msoln) {
        if (cmods) msoln->match(*cmods);
        if (grid)  msoln->match(*grid);
    }

    if (msoln) {
        INFO0(who, "Restart file prescribes a manufactured solution");
        if (!(isnan)(chdef->bulk_rho)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho = "
                  << chdef->bulk_rho);
        }
        if (!(isnan)(chdef->bulk_rho_u)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_u = "
                  << chdef->bulk_rho_u);
        }
    }

    // If requested, add noise to the momentum fields at startup (expensive).
    if (noisedef.percent > 0) {
        state_nonlinear->assign_from(*state_linear);
        add_noise(*state_nonlinear, noisedef,
                  *grid, *dgrid, *cop, *b);
        state_linear->assign_from(*state_nonlinear);
    }


    // Prepare any necessary, problem-specific constraints
    shared_ptr<constraint::treatment> constrainer(new constraint::treatment(
                    1.0, *dgrid, *b, common_block, common_block));
    if        (grid->two_sided()) { // Channel per treatment_channel.tex

        INFO0(who, "Establishing driving, channel-like state constraints");
        constrainer->physical[ndx::rho].reset(
                new constraint::reference_bulk(chdef->bulk_rho  , *b));
        constrainer->physical[ndx::mx ].reset(
                new constraint::reference_bulk(chdef->bulk_rho_u, *b));
        constrainer->physical[ndx::e  ].reset(
                new constraint::reference_bulk(chdef->bulk_rho_E, *b));

    } else if (grid->one_sided()) { // Flat plate

        real_t rho_inf;
        real_t mx_inf;
        real_t e_inf;

        // Compute baseflow from coefficients
        if (sgdef->baseflow && !preserve_upper_input) {
            INFO0(who, "Setting mean freestream reference state per baseflow");
            WARN0(who, "... there is no correction to freestream values for"
                       " displacement effects ");
            if (isothermal->lower_v != 0.) {
                WARN0(who, "... there is no correction to freestream energy"
                           " for wall-normal velocity");
            }

            const size_t Ns = cmods->Ns();
            if (Ns > 1) {
                WARN0(who, "... setting mean freestream reference state"
                      " per baseflow NOT tested for multispecies,"
                      " proceed at your own risk");
            }
            real_t   base [Ns+4+1];
            real_t dybase [Ns+4+1];
            real_t dxbase [Ns+4+1];
            sgdef->baseflow->conserved(
                grid->L.y(), &base[0], &dybase[0], &dxbase[0]);
            sgdef->baseflow->pressure(
                grid->L.y(), base[Ns+4], dybase[Ns+4], dxbase[Ns+4]);

            // Clobber upper values except upper verocity
            // FIXME: Account for the correction to upper_v for the 
            //        modified characteristics from slow growth.
            //  NOTE: Maybe the proper place to implement the upper_v 
            //        correction is in the nonreflecting implementation,
            //        instead of setting a corrected upper_v here.
            //        Working with a modified upper velocity here might be 
            //        confusing and lead to wrong results, when computing
            //        derived quantities (eg, temperature), that depend on
            //        the physical upper_v (not the effective one from 
            //        characteristics)
            isothermal->upper_rho = base[0];
            isothermal->upper_u   = base[1]/base[0];
            isothermal->upper_w   = base[3]/base[0];

            for (unsigned int s=1; s<Ns; ++s) {
              isothermal->upper_cs[s] = base[5+s]/base[0];
            }

            // Get upper temperature from antioch
            const real_t irho = 1.0/isothermal->upper_rho ;

            // ... momentum
            // ... use the actual freestream rho_v value at the upper boundary
            //     to compute temperature properly
            const real_t   upper_rho_v  (base[2]);
            const Vector3r m (isothermal->upper_rho * isothermal->upper_u,
                              upper_rho_v,
                              isothermal->upper_rho * isothermal->upper_w);

            // ... total energy
            const real_t e   (base[4]);

            // ... species densities
            VectorXr species(Ns); // species densities
            VectorXr cs     (Ns); // species mass fractions

            // NOTE: In species vector, idx 0 is the dilluter (the species
            // that is not explicitly part of the state vector)
            species(0) = isothermal->upper_rho;
            for (unsigned int s=1; s<Ns; ++s) {
                species(s)  = isothermal->upper_cs[s];

                // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
                species(0) -= species(s);
            }

            // ... compute mass fractions
            for (unsigned int s=0; s<Ns; ++s) {
                cs(s) = irho * species(s);
            }

            // ... temperature and pressure-related quantities
            real_t T, pr;
            cmods->evaluate_pressure(
                e, m, isothermal->upper_rho, species, cs, isothermal->lower_T,
                T, pr);

            // ... finally we got what we wanted...
            isothermal->upper_T = T;

            INFO0(who, "... upper_rho set to " << isothermal->upper_rho );
            INFO0(who, "... upper_u   set to " << isothermal->upper_u   );
            INFO0(who, "... upper_v   set to " << isothermal->upper_v   );
            INFO0(who, "... upper_w   set to " << isothermal->upper_w   );
            INFO0(who, "... upper_T   set to " << isothermal->upper_T   );

            // ... set inf values for freestream constraints
            rho_inf = isothermal->upper_rho;
            mx_inf  = isothermal->upper_u   * rho_inf;
            e_inf   = e;

        } else {
            INFO0(who, "Computing mean freestream behavior per plate scenario");
            if (sgdef->baseflow && preserve_upper_input) {
                WARN0(who, "... preserving upper input or file values as per"
                       " 'preserve_upper_input' option");
            }
            const real_t T_inf   = isothermal->upper_T;
            const real_t u_inf   = isothermal->upper_u;
            rho_inf = isothermal->upper_rho;
            mx_inf  = u_inf * rho_inf;

            // Compute total energy at infinity based on input upper values
            // Neglect contribution of wall-normal velocity to kinetic energy,
            // and assume that species concentrations are at equilibrium
            // FIXME: adjust the value of target rho_E dynamically
            //        to account for wall-normal velocity and variations in
            //        species mass fractions
            const size_t Ns = cmods->Ns();
            vector<real_t> mass_fractions(Ns);
            for (unsigned int s=0; s<Ns; ++s) {
              mass_fractions[s]      = isothermal->upper_cs[s];
            }
            e_inf   = (cmods->e_from_T(T_inf, mass_fractions)
                + 0.5 * (u_inf * u_inf)) * rho_inf;
       }

        INFO0(who, "Setting constraints using freestream reference state"
                   " on upper boundary");

        if (chdef->bulk_rho) {
            WARN0(who, "Removing channel-like bulk_rho setting");
            chdef->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        }
        if (chdef->bulk_rho_u) {
            WARN0(who, "Removing channel-like bulk_rho_u setting");
            chdef->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
        }
        if (chdef->bulk_rho_E) {
            WARN0(who, "Removing channel-like bulk_rho_E setting");
            chdef->bulk_rho_E = numeric_limits<real_t>::quiet_NaN();
        }

        INFO0(who, "Establishing driving, freestream-like state constraints");
        constrainer->physical[ndx::e  ].reset(
                new constraint::constant_upper(e_inf,   *cop, 0));
        constrainer->physical[ndx::mx ].reset(
                new constraint::constant_upper(mx_inf,  *cop, 0));
        constrainer->physical[ndx::rho].reset(
                new constraint::constant_upper(rho_inf, *cop, 0));

    } else {

        FATAL0(who, "Sanity error in constraint selection");
        return EXIT_FAILURE;

    }


    // Nonreflecting must mutate chdef/isothermal before instantiating
    // linear operator
    if (grid->one_sided()) {
        // FIXME: Plate implementation in progress.
        // Need to set reference value for rho (rho_ref)

        INFO0(who, "Preparing nonreflecting upper boundary treatment");

        WARN0(who, "Non-reflecting boundary treatment enabling boundary"
                   " layer simulations is still under development.");


        // Assign some reference values
        const size_t Ns = cmods->Ns();

        // SUZERAIN_ENSURE(!(isnan)(isothermal->upper_rho));
        if ((isnan)(isothermal->upper_rho)) {
            FATAL0(who, "upper_rho for nonreflecting boundary"
                        " not specified and not in restart file");
            return EXIT_FAILURE;
        }

        common_block.rho_ref = isothermal->upper_rho;
        common_block.T_ref   = isothermal->upper_T;
        common_block.u_ref   = isothermal->upper_u;
        common_block.v_ref   = isothermal->upper_v;
        common_block.w_ref   = isothermal->upper_w;

        common_block.cs_ref.resize(Ns);
        vector<real_t> mass_fractions(Ns);
        for (unsigned int s=0; s<Ns; ++s) {
            common_block.cs_ref(s) = isothermal->upper_cs[s];
            mass_fractions[s]      = isothermal->upper_cs[s];
        }

        common_block.etots_ref.resize(Ns);
        cmods->evaluate_for_nonreflecting(common_block.T_ref,
                                          common_block.cs_ref,
                                          common_block.a_ref,
                                          common_block.gamma_ref,
                                          common_block.R_ref,
                                          common_block.etots_ref);

        common_block.Cv_ref =  common_block.R_ref
                            / (common_block.gamma_ref - 1);

        // FIXME: Remove this info or make it a debug output option
        INFO0(who, "rho_ref   = " << common_block.rho_ref  );
        INFO0(who, "T_ref     = " << common_block.T_ref    );
        INFO0(who, "a_ref     = " << common_block.a_ref    );
        INFO0(who, "gamma_ref = " << common_block.gamma_ref);
        INFO0(who, "R_ref     = " << common_block.R_ref    );
        INFO0(who, "E_ref     = " << cmods->e_from_T(common_block.T_ref,
                                                     mass_fractions));
    }

    // Reset amplitude growth rate parameters for spatiotemporal formulation
    if (sgdef->formulation == largo_formulation::spatiotemporal_consistent
        && !preserve_largo_gramp_input) {
        INFO0(who, "Setting largo growth rate amplitude values");

        const size_t Ns = cmods->Ns();

        // resize to 5 flow vars + (Ns-1) species + pressure
        sgdef->gramp_mean.resize(5+(Ns-1)+1);

        // set all values to zero
        for (unsigned int i=0; i<sgdef->gramp_mean.size(); i++) {
            sgdef->gramp_mean[i] = 0.;
        }

        // pre-computations to get wall baseflow values
        const size_t iu  = 1;
        const size_t ie  = 4;
        const size_t is0 = 5;
        const size_t ip  = Ns + 4;
        real_t   wall [Ns+4+1];
        real_t dywall [Ns+4+1];
        real_t dxwall [Ns+4+1];
        sgdef->baseflow->conserved(
            0., &wall[0], &dywall[0], &dxwall[0]);
        sgdef->baseflow->pressure(
            0., wall[ip], dywall[ip], dxwall[ip]);


        // thermo quantities to estimate wall flow values
        VectorXr wall_cs(Ns);
        for (unsigned int s=0; s<Ns; ++s) {
            wall_cs(s) = isothermal->lower_cs[s];
        }

        real_t wall_a;            // not needed
        real_t wall_gamma;
        real_t wall_Rmix;
        VectorXr wall_etots(Ns);  // not needed
        cmods->evaluate_for_nonreflecting(isothermal->lower_T,
                                          wall_cs,
                                          wall_a,
                                          wall_gamma,
                                          wall_Rmix,
                                          wall_etots);

        // pre-compute these for convenience
        const real_t wall_Cv       = wall_Rmix / (wall_gamma - 1.);
        const real_t wall_P        =   wall[ip];
        const real_t wall_ddx_P    = dxwall[ip];
        const real_t wall_rhoI     =   wall[0];
        const real_t wall_ddx_rhoI = dxwall[0];
        const real_t wall_rho      = wall_P     / wall_Rmix / isothermal->lower_T;
        const real_t wall_ddx_rho  = wall_ddx_P / wall_Rmix / isothermal->lower_T;
        const real_t wall_ruI      =   wall[iu];
        const real_t wall_ddx_ruI  = dxwall[iu];
        const real_t wall_EI       =   wall[ie] / wall[0];
        const real_t wall_rEI      =   wall[ie];
        const real_t wall_ddx_rEI  = dxwall[ie];
        const real_t wall_E        = wall_Cv * isothermal->lower_T;

        // growth rate amplitudes for density, streamwise velocity, energy
        sgdef->gramp_mean[0]  = 1./(wall_rho - wall_rhoI) * (wall_ddx_rho - wall_ddx_rhoI);
        sgdef->gramp_mean[1]  = wall_ddx_ruI/wall_ruI - wall_ddx_rhoI/wall_rhoI;
        sgdef->gramp_mean[ie] = - (wall_EI/(wall_E - wall_EI))
                                * (wall_ddx_rEI/wall_rEI - wall_ddx_rhoI/wall_rhoI);

        // growth rate amplitude for species must be set all equal for consistency
        // setting them all equal to zero
        for (size_t s=0; s<Ns-1; s++) {
            sgdef->gramp_mean[is0+s] = sgdef->gramp_mean[0];
        }

        // rms amplitude growth rates
        sgdef->gramp_rms.resize(5+(Ns-1)+1);


        // output computed values
        DEBUG0(who, "... wall rho  is " << wall_rho);
        DEBUG0(who, "... wall Cv   is " << wall_Cv);
        DEBUG0(who, "... wall Rmix is " << wall_Rmix);
        DEBUG0(who, "... wall E    is " << wall_E);

        INFO0(who, "... gramp_mean for rho set to " << sgdef->gramp_mean[0]);
        INFO0(who, "... gramp_mean for u   set to " << sgdef->gramp_mean[iu]);
        INFO0(who, "... gramp_mean for E   set to " << sgdef->gramp_mean[ie]);
        WARN0(who, "... gramp_mean for cs  set all to 0. ");

        // set all values to zero
        for (unsigned int i=0; i<sgdef->gramp_rms.size(); i++) {
            sgdef->gramp_rms[i] = 0.;
        }
    } else if (sgdef->formulation == largo_formulation::spatiotemporal_consistent
              && preserve_largo_gramp_input) {
        WARN0(who, "Preserving largo_gramp values from file"
                   " per 'preserve_largo_gramp_input' option");
    }

    // Prepare spatial operators depending on requested advance type
    // Notice constrainer always wraps the implicit operator
    // and that the same operator_nonlinear is used pervasively.
    L = constrainer;
    if (use_explicit) {
        INFO0(who, "Initializing fully explicit spatial operators");

        N.reset(new explicit_operator_nonlinear(
                    *cmods, *grid, *dgrid, *cop, *b, common_block, *fsdef,
                    *sgdef, msoln));

        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<treatment_nonreflecting> nonreflecting(
                    new treatment_nonreflecting(
                        *grid, *dgrid, *cop, *b, common_block));
            nonreflecting->N = N;
            N = nonreflecting;
        }

        constrainer->L.reset(new operator_mass_isothermal(
                    *cmods, *isothermal,
                    *chdef, *grid, *dgrid, *cop, *b, common_block));

    } else if (use_implicit) {
        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");

        constrainer->L.reset(new isothermal_hybrid_linear_operator(
                    solver_spec, *cmods, *isothermal, *chdef, *grid, *dgrid,
                    *cop, *b, common_block));

        shared_ptr<suzerain::operator_hybrid_residual>
            tmp_hybrid( new operator_hybrid_residual(dgrid->chi()) );

        tmp_hybrid->R.reset(new explicit_operator_nonlinear(
                                *cmods, *grid, *dgrid, *cop, *b,
                                common_block, *fsdef, *sgdef, msoln));

        tmp_hybrid->L = this->L;

        this->N = tmp_hybrid;

        if (grid->one_sided()) {
            WARN0(who, "Non-reflecting boundary treatment with hybrid"
                       " implicit/explicit time marching is experimental.");

            shared_ptr<treatment_nonreflecting> nonreflecting(
                new treatment_nonreflecting(
                    *grid, *dgrid, *cop, *b, common_block));
            nonreflecting->N = this->N;
            this->N = nonreflecting;
        }


    } else {
        FATAL0(who, "Sanity error in operator selection");
        return EXIT_FAILURE;
    }

    // Allocate slow growth workspace
    if (sgdef->formulation.enabled()) {
        const std::string& model = sgdef->formulation.name();
        const int neq            = static_cast<int>(cmods->Ns()) + 4;
        const int ns             = static_cast<int>(cmods->Ns()) - 1;
        enum  { ntvar = 0 };
        static const char ransmodel[] = "dns";
        INFO0(who, "Allocating Largo model \"" << model << "\" with neq="
              << neq << ", ns=" << ns << ", ntvar=" << ntvar
              << ", ransmodel=" << ransmodel);
        largo_allocate(&sgdef->workspace, model.c_str(), neq,
                       ns, ntvar, ransmodel);
        if (!sgdef->workspace) {
            FATAL0(who, "Largo could not allocate requested model");
            return EXIT_FAILURE;
        }
        if ((isnan)(sgdef->grdelta)) {
            WARN0(who, "Slow growth rate grdelta is NaN");
        }
        if (sgdef->ignore_fluctuations) {
            // TODO Permit option by modifying operator_nonlinear.hpp.
            FATAL0(who, "Debugging option 'ignore_fluctuations' not supported");
            return EXIT_FAILURE;
        }
        if (sgdef->ignore_gramp_mean) {
            // TODO Permit option by modifying operator_nonlinear.hpp.
            FATAL0(who, "Debugging option 'ignore_gramp_mean' not supported");
            return EXIT_FAILURE;
        }
        if (sgdef->ignore_gramp_rms) {
            // TODO Permit option by modifying operator_nonlinear.hpp.
            FATAL0(who, "Debugging option 'ignore_gramp_rms' not supported");
            return EXIT_FAILURE;
        }
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
    INFO0(who, "B-spline basis of order " << grid->k
               << " on [" << b->collocation_point(0) << ", "
               << b->collocation_point(b->n() - 1) << "] with "
               << b->n() << " DOF stretched per htdelta " << grid->htdelta);
    log_discretization_quality();
    prepare_controller(initial_t, dgrid->chi());
    save_metadata();
    const real_t elapsed_wall_time // Negative on error
        = advance_controller(timedef->status_final,
                             statsdef->final,
                             restartdef->final);

    // Report error to the OS iff advance_control reported an error
    return elapsed_wall_time < 0;
}
