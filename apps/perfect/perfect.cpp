//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// perfect.cpp: Support logic for the Suzerain perfect gas application
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/rngstream.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/validation.hpp>
#include <sys/file.h>

#include "perfect.hpp"
#include "../logging.hpp"
#include "../support.hpp"

// Manufactured solution classes explicitly instantiated for debugging
template class nsctpl_rholut::manufactured_solution<suzerain::real_t>;

using boost::numeric_cast;
using std::size_t;

namespace suzerain { namespace perfect {

std::vector<support::field> default_fields()
{
    std::vector<support::field> retval(ndx::identifier.static_size);
    for (size_t i = 0; i < ndx::identifier.static_size; ++i) {
        retval[i].identifier   = ndx::identifier[i];
        retval[i].description += "Nondimensional ";
        retval[i].description += ndx::description[i];
        retval[i].location     = ndx::identifier[i];
    }
    return retval;
}

void store(const esio_handle h,
           const ScenarioDefinition& scenario)
{
    DEBUG0("Storing ScenarioDefinition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Re", &scenario.Re, 0,
            scenario.options().find("Re",false).description().c_str());

    esio_line_write(h, "Ma", &scenario.Ma, 0,
            scenario.options().find("Ma",false).description().c_str());

    esio_line_write(h, "Pr", &scenario.Pr, 0,
            scenario.options().find("Pr",false).description().c_str());

    esio_line_write(h, "bulk_rho", &scenario.bulk_rho, 0,
            scenario.options().find("bulk_rho",false).description().c_str());

    esio_line_write(h, "bulk_rho_u", &scenario.bulk_rho_u, 0,
            scenario.options().find("bulk_rho_u",false).description().c_str());

    esio_line_write(h, "alpha", &scenario.alpha, 0,
            scenario.options().find("alpha",false).description().c_str());

    esio_line_write(h, "beta", &scenario.beta, 0,
            scenario.options().find("beta",false).description().c_str());

    esio_line_write(h, "gamma", &scenario.gamma, 0,
            scenario.options().find("gamma",false).description().c_str());
}

void load(const esio_handle h,
          ScenarioDefinition& scenario)
{
    DEBUG0("Loading ScenarioDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (!(boost::math::isnan)(scenario.Re)) {
        INFO0("Overriding scenario using Re = " << scenario.Re);
    } else {
        esio_line_read(h, "Re", &scenario.Re, 0);
    }

    if (!(boost::math::isnan)(scenario.Ma)) {
        INFO0("Overriding scenario using Ma = " << scenario.Ma);
    } else {
        esio_line_read(h, "Ma", &scenario.Ma, 0);
    }

    if (!(boost::math::isnan)(scenario.Pr)) {
        INFO0("Overriding scenario using Pr = " << scenario.Pr);
    } else {
        esio_line_read(h, "Pr", &scenario.Pr, 0);
    }

    if (!(boost::math::isnan)(scenario.bulk_rho)) {
        INFO0("Overriding scenario using bulk_rho = " << scenario.bulk_rho);
    } else {
        esio_line_read(h, "bulk_rho", &scenario.bulk_rho, 0);
    }

    if (!(boost::math::isnan)(scenario.bulk_rho_u)) {
        INFO0("Overriding scenario using bulk_rho_u = " << scenario.bulk_rho_u);
    } else {
        esio_line_read(h, "bulk_rho_u", &scenario.bulk_rho_u, 0);
    }

    if (!(boost::math::isnan)(scenario.alpha)) {
        INFO0("Overriding scenario using alpha = " << scenario.alpha);
    } else {
        esio_line_read(h, "alpha", &scenario.alpha, 0);
    }

    if (!(boost::math::isnan)(scenario.beta)) {
        INFO0("Overriding scenario using beta = " << scenario.beta);
    } else {
        esio_line_read(h, "beta", &scenario.beta, 0);
    }

    if (!(boost::math::isnan)(scenario.gamma)) {
        INFO0("Overriding scenario using gamma = " << scenario.gamma);
    } else {
        esio_line_read(h, "gamma", &scenario.gamma, 0);
    }
}

static void attribute_storer(const esio_handle &h,
                             const char *location,
                             const std::string &name,
                             const real_t &value)
{
    esio_attribute_write(h, location, name.c_str(), &value);
}

void store(const esio_handle h,
           const ScenarioDefinition& scenario,
           const problem::GridDefinition& grid,
           const boost::shared_ptr<manufactured_solution>& msoln)
{
    // Only proceed if a manufactured solution is being provided
    if (!msoln) return;

    static const char location[] = "channel::manufactured_solution";
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
            "Is a channel::manufactured_solution in use?");

    DEBUG0("Storing channel::manufactured_solution parameters");

#pragma warning(push,disable:1572)
    // Check parameters stored with the scenario not the manufactured solution
    // because scenario parameters should be loaded from ScenarioDefinition
    if (msoln->alpha != scenario.alpha)
        WARN0("Manufactured solution alpha mismatches with scenario!");
    if (msoln->beta  != scenario.beta)
        WARN0("Manufactured solution beta mismatches with scenario!");
    if (msoln->gamma != scenario.gamma)
        WARN0("Manufactured solution gamma mismatches with scenario!");
    if (msoln->Ma    != scenario.Ma)
        WARN0("Manufactured solution Ma mismatches with scenario!");
    if (msoln->Re    != scenario.Re)
        WARN0("Manufactured solution Re mismatches with scenario!");
    if (msoln->Pr    != scenario.Pr)
        WARN0("Manufactured solution Pr mismatches with scenario!");

    // Check parameters stored with the grid not the manufactured solution
    // because grid parameters should be loaded from GridDefinition
    if (msoln->Lx    != grid.L.x())
        WARN0("Manufactured solution Lx mismatches with grid!");
    if (msoln->Ly    != grid.L.y())
        WARN0("Manufactured solution Ly mismatches with grid!");
    if (msoln->Lz    != grid.L.z())
        WARN0("Manufactured solution Lz mismatches with grid!");
#pragma warning(pop)

    // Non-scenario solution parameters are stored as attributes under location
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->u  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->v  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->w  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->T  .foreach_parameter(bind(attribute_storer, h, location, _1, _2));
}

static void attribute_loader(const esio_handle &h,
                             const char *location,
                             const std::string &name,
                             real_t &value)
{
    esio_attribute_read(h, location, name.c_str(), &value);
}

/** Helper for NaNing values within a manufactured solution instance */
static void NaNer(const std::string&, real_t& value)
{
    value = std::numeric_limits<real_t>::quiet_NaN();
}

void load(const esio_handle h,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          boost::shared_ptr<manufactured_solution>& msoln)
{
    static const char location[] = "channel::manufactured_solution";

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

    DEBUG0("Loading channel::manufactured_solution parameters");

    // Allocate storage and defensively NaN every parameter not explicitly
    // loaded below.  Protects us against accidentally missing new solution
    // parameters.
    msoln.reset(new manufactured_solution());
    msoln->foreach_parameter(&NaNer);

    // Scenario parameters taken from ScenarioDefinition
    msoln->alpha = scenario.alpha;
    msoln->beta  = scenario.beta;
    msoln->gamma = scenario.gamma;
    msoln->Ma    = scenario.Ma;
    msoln->Re    = scenario.Re;
    msoln->Pr    = scenario.Pr;

    // Grid parameters taken from GridDefinition
    msoln->Lx    = grid.L.x();
    msoln->Ly    = grid.L.y();
    msoln->Lz    = grid.L.z();

    // Non-scenario solution parameters are stored as attributes under location
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->u  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->v  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->w  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->T  .foreach_parameter(bind(attribute_loader, h, location, _1, _2));
}

void store_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& swave,
        const ScenarioDefinition& scenario,
        const problem::GridDefinition& grid,
        const pencil_grid& dgrid,
        bspline& b,
        const bsplineop& bop)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Convert coefficients into collocation point values
    // Transforms state from full-wave coefficients to full-physical points
    operator_base obase(grid, dgrid, b, bop);
    for (size_t i = 0; i < swave.shape()[0]; ++i) {
        obase.bop_apply(0, 1, swave, i);
        obase.zero_dealiasing_modes(swave, i);
        dgrid.transform_wave_to_physical(
                reinterpret_cast<real_t *>(swave[i].origin()));
    }

    // Convert conserved rho{_E,_u,_v,_w,} into u, v, w, p, T
    support::physical_view<state_count>::type sphys
        = support::physical_view<state_count>::create(dgrid, swave);

    const real_t alpha = scenario.alpha;
    const real_t beta  = scenario.beta;
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack conserved quantities from fields
        const real_t     e(sphys(ndx::e,   o));
        Vector3r         m(sphys(ndx::mx,  o),
                           sphys(ndx::my,  o),
                           sphys(ndx::mz,  o));
        const real_t   rho(sphys(ndx::rho, o));

        // Compute primitive quantities to be stored
        real_t p, T;
        rholut::p_T(alpha, beta, gamma, Ma, rho, m, e, p, T);
        m /= rho;

        // Pack primitive quantities back into fields (by position)
        sphys(0, o) = m.x(); // Now just X velocity
        sphys(1, o) = m.y(); // Now just Y velocity
        sphys(2, o) = m.z(); // Now just Z velocity
        sphys(3, o) = p;
        sphys(4, o) = T;
    }

    // HDF5 file storage locations and corresponding descriptions
    const boost::array<const char *,state_count> prim_names = {{
        "u", "v", "w", "p", "T"
    }};
    const boost::array<const char *,state_count> prim_descriptions = {{
        "streamwise velocity",
        "wall-normal velocity",
        "spanwise velocity",
        "pressure",
        "temperature",
    }};
    assert(prim_names.static_size == prim_descriptions.static_size);

    // Establish size of collective writes across all ranks and write fields
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());

    for (size_t i = 0; i < state_count; ++i) {

        std::string comment = "Nondimensional ";
        comment += prim_descriptions[i];
        comment += " stored row-major YZX on the 3D rectilinear grid defined"
                   " by taking the outer product of arrays"
                   " /collocation_points_y, /collocation_points_z, and"
                   " /collocation_points_z";

        esio_field_write(h, prim_names[i],
                reinterpret_cast<real_t *>(swave[i].origin()),
                0, 0, 0, comment.c_str());
    }
}

void load_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& state,
        const ScenarioDefinition& scenario,
        const problem::GridDefinition& grid,
        const pencil_grid& dgrid,
        bspline& b,
        const bsplineop& bop)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // This routine does no grid interpolation.  Yell loudly if necessary
    {
        // Check that restart file size matches runtime dealiased extents
        int cg, bg, ag;
        if (ESIO_SUCCESS != esio_field_size(h, "u", &cg, &bg, &ag)) {
            SUZERAIN_ERROR_VOID("Unable to find /u field size from restart",
                                SUZERAIN_EFAILED);
        }
        if (cg != grid.dN.y() || bg != grid.dN.z() || ag != grid.dN.x()) {
            ERROR0("Physical space restart fields have row-major YZX extents "
                   << "(" << cg << "," << bg << "," << ag << ")" << " but "
                   << "(" << grid.dN.y() << "," << grid.dN.z() << ","
                   << grid.dN.x() << ") are required");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }


        // Check that restart file specifies the same B-spline basis.
        // TODO Too restrictive!  Identical collocation points would be okay.
        // TODO Too restrictive?  Any floating point differences kill us.
        boost::shared_ptr<bspline> Fb;
        boost::shared_ptr<bsplineop> Fbop;
        support::load(h, Fb, Fbop);
        const double bsplines_dist = support::distance(b, *Fb);
        const bool bsplines_same
                = bsplines_dist < support::bsplines_distinct_distance;
        if (!bsplines_same) {
            ERROR0("Physical restart has different wall-normal bases ("
                   << bsplines_dist << ")");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }
    }

    // Establish size of collective reads across all ranks and read data
    support::physical_view<state_count>::type sphys
        = support::physical_view<state_count>::create(dgrid, state);
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());
    esio_field_read(h, "u", &sphys(0,0), 0, 0, 0);
    esio_field_read(h, "v", &sphys(1,0), 0, 0, 0);
    esio_field_read(h, "w", &sphys(2,0), 0, 0, 0);
    esio_field_read(h, "p", &sphys(3,0), 0, 0, 0);
    esio_field_read(h, "T", &sphys(4,0), 0, 0, 0);

    // Convert primitive u, v, w, p, and T into rho{_E,_u,_v,_w,}
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack primitive quantities from fields (by position)
        Vector3r       m(sphys(0, o),  // Now just X velocity
                         sphys(1, o),  // Now just Y velocity
                         sphys(2, o)); // Now just Z velocity
        const real_t   p(sphys(3, o));
        const real_t   T(sphys(4, o));

        // Compute conserved quantities from primitive ones
        const real_t rho = gamma * p / T;   // Assumes EOS
        m               *= rho;             // Now m contains momentum
        const real_t e   = rholut::energy_kinetic(Ma, rho, m)
                         + rholut::energy_internal(gamma, p);

        // Pack conserved quantities into fields (by name)
        sphys(ndx::e,   o) = e;
        sphys(ndx::mx,  o) = m.x();
        sphys(ndx::my,  o) = m.y();
        sphys(ndx::mz,  o) = m.z();
        sphys(ndx::rho, o) = rho;
    }

    // Initialize operator_base to access decomposition-ready utilities
    operator_base obase(grid, dgrid, b, bop);

    // Collectively convert physical state to wave space coefficients
    // Build FFT normalization constant into Y direction's mass matrix
    bsplineop_luz massluz(bop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, bop);
    massluz.factor();

    for (size_t i = 0; i < state.shape()[0]; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0)); // X, Z
        obase.zero_dealiasing_modes(state, i);
        obase.bop_solve(massluz, state, i);                      // Y
    }
}

void load(const esio_handle h,
          contiguous_state<4,complex_t>& state,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          const pencil_grid& dgrid,
          bspline& b,
          const bsplineop& bop)
{
    const std::vector<support::field> fields = default_fields();

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(state.shape()[0]  == fields.size());

    // Check whether load_coefficients(...) should work
    bool trycoeffs = true;
    for (size_t i = 0; i < fields.size(); ++i) {
        int ncomponents = 0;
        switch (esio_field_sizev(h, fields[i].location.c_str(),
                                 0, 0, 0, &ncomponents))
        {
            case ESIO_SUCCESS:
                if (ncomponents != 2) {
                    DEBUG0("Field /" << fields[i].identifier
                                     << " looks fishy...");
                }
                break;
            case ESIO_NOTFOUND:
                DEBUG0("Field /" << fields[i].identifier
                                 << " not found in restart");
                trycoeffs = false;
                break;
            default:
                DEBUG0("Field /" << fields[i].identifier << " looks fishy...");
                break;
        }
    }

    // Dispatch to the appropriate restart loading logic
    DEBUG0("Started loading simulation fields");
    if (trycoeffs) {
        support::load_coefficients(h, fields, state, grid, dgrid, b, bop);
    } else {
        INFO0("Loading collocation-based, physical-space restart data");
        load_collocation_values(h, state, scenario, grid, dgrid, b, bop);
    }
    DEBUG0("Finished loading simulation fields");
}

void
adjust_scenario(contiguous_state<4,complex_t> &swave,
                const ScenarioDefinition& scenario,
                const problem::GridDefinition& grid,
                const pencil_grid& dgrid,
                bspline &b,
                const bsplineop& bop,
                const real_t old_Ma,
                const real_t old_gamma)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    bool quickreturn = true;
#pragma warning(push,disable:1572)
    if (old_Ma != scenario.Ma) {
#pragma warning(pop)
        INFO0("Changing state Mach number from "
              << old_Ma << " to " << scenario.Ma);
        quickreturn = false;
    }
#pragma warning(push,disable:1572)
    if (old_gamma != scenario.gamma) {
#pragma warning(pop)
        INFO0("Changing state ratio of specific heats from "
              << old_gamma << " to " << scenario.gamma);
        quickreturn = false;
    }
    if (quickreturn) {
        DEBUG0("Not rescaling energy as scenario parameters have not changed");
        return;
    }

    // Convert state to physical space collocation points
    operator_base obase(grid, dgrid, b, bop);
    support::physical_view<state_count>::type sphys
            = support::physical_view<state_count>::create(dgrid, swave);
    for (size_t k = 0; k < state_count; ++k) {
        obase.bop_apply(0, 1.0, swave, k);
        dgrid.transform_wave_to_physical(&sphys.coeffRef(k,0));
    }

    // Adjust total energy by the necessary amount at every collocation point
    // This procedure is not the cheapest or least numerically noisy,
    // but it does re-use existing compute kernels in a readable way.
    INFO0("Holding density and temperature constant during changes");
    size_t offset = 0;
    for (int j = dgrid.local_physical_start.y();
        j < dgrid.local_physical_end.y();
        ++j) {
        const size_t last_zxoffset = offset
                                   + dgrid.local_physical_extent.z()
                                   * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            const real_t   e  (sphys(ndx::e,   offset));
            const Vector3r m  (sphys(ndx::mx,  offset),
                               sphys(ndx::my,  offset),
                               sphys(ndx::mz,  offset));
            const real_t   rho(sphys(ndx::rho, offset));

            // Compute temperature using old_gamma, old_Ma
            real_t p, T;
            rholut::p_T(scenario.alpha, scenario.beta, old_gamma, old_Ma,
                        rho, m, e, /*out*/ p, /*out*/ T);
            // Compute total energy from new gamma, Ma, rho, T
            rholut::p(scenario.gamma, rho, T, /*out*/ p);
            sphys(ndx::e, offset) = rholut::energy_internal(scenario.gamma, p)
                                  + rholut::energy_kinetic(scenario.Ma, rho, m);
        }
    }

    // Convert state back to wave space coefficients in X, Y, and Z
    // building FFT normalization constant into the mass matrix
    bsplineop_luz massluz(bop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, bop);
    massluz.factor();
    for (size_t i = 0; i < state_count; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0));
        obase.bop_solve(massluz, swave, i);
    }
}

NoiseDefinition::NoiseDefinition(real_t percent,
                                 unsigned long seed)
    : definition_base("Additive random velocity perturbations on startup"),
      percent(percent),
      kxfrac_min(0),
      kxfrac_max(1),
      kzfrac_min(0),
      kzfrac_max(1),
      seed(seed)
{
    using boost::bind;
    using validation::ensure_positive;
    using validation::ensure_nonnegative;
    std::pointer_to_binary_function<unsigned long,const char*,void>
            ptr_fun_ensure_positive_ulint(ensure_positive<unsigned long>);
    std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative_real(ensure_nonnegative<real_t>);
    this->add_options()
        ("fluct_percent",
         boost::program_options::value(&this->percent)
            ->default_value(this->percent)
            ->notifier(std::bind2nd(ptr_fun_ensure_nonnegative_real,
                                   "fluct_percent")),
         "Maximum fluctuation magnitude to add as a percentage of"
         " centerline mean streamwise velocity")
        ("fluct_kxfrac",
         boost::program_options::value<std::string>(0)
            ->default_value("0:1")
            ->notifier(bind(&support::parse_range<real_t>, _1,
                            &this->kxfrac_min, &this->kxfrac_max,
                            0, 1, 0, 1, "fluct_kxfrac")),
         "Range of X wavenumbers in which to generate fluctuations")
        ("fluct_kzfrac",
         boost::program_options::value<std::string>(0)
            ->default_value("0:1")
            ->notifier(bind(&support::parse_range<real_t>, _1,
                            &this->kzfrac_min, &this->kzfrac_max,
                            0, 1, 0, 1, "fluct_kzfrac")),
         "Range of Z wavenumbers in which to generate fluctuations")
        ("fluct_seed",
         boost::program_options::value(&this->seed)
            ->default_value(this->seed)
            ->notifier(std::bind2nd(ptr_fun_ensure_positive_ulint,
                                    "fluct_seed")),
         "RngStream generator seed (L'Ecuyer et al. 2002)");
}

void
add_noise(contiguous_state<4,complex_t> &state,
          const NoiseDefinition& noisedef,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          const pencil_grid& dgrid,
          bspline &b,
          const bsplineop& bop)
{
    // FIXME Needs to be made aware of channel vs flat plate

    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed collocation-based operator matrices
    SUZERAIN_ENSURE(bop.get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    const int Ny = grid.N.y();

#pragma warning(push,disable:1572)
    if (noisedef.percent == 0) {
#pragma warning(pop)
        DEBUG0("Zero noise added to velocity fields");
        return;
    }

    using inorder::wavenumber_abs;
    using inorder::wavenumber_max;
    using inorder::wavenumber_translatable;
    const real_t twopi = 2 * boost::math::constants::pi<real_t>();

    // Evaluate maximum fluctuation magnitude based on percentage of
    // centerline mean streamwise velocity and broadcast result.
    // Alright, alright... actually an approximate mean velocity.
    real_t maxfluct;
    if (dgrid.has_zero_zero_modes()) {
        complex_t momentum, density;
        const real_t centerline = grid.L.y() / 2;
        b.linear_combination(
                0, &state[ndx::mx][0][0][0], 1, &centerline, &momentum);
        INFO0("Centerline mean streamwise momentum at y = "
              << centerline << " is " << momentum);
        b.linear_combination(
                0, &state[ndx::rho][0][0][0], 1, &centerline, &density);
        INFO0("Centerline mean density at y = "
              << centerline << " is " << density);
        maxfluct = noisedef.percent / 100 * (abs(momentum) / abs(density));
    }
    SUZERAIN_MPICHKR(MPI_Bcast(&maxfluct, 1,
                mpi::datatype_of(maxfluct),
                dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
    INFO0("Adding velocity perturbations with maximum magnitude " << maxfluct);

    // Compute and display kxfrac_min, kxfrac_max constraints
#pragma warning(push,disable:2259)
    const int dkx_max = noisedef.kxfrac_max * wavenumber_max(grid.N.x());
    const int dkx_min = noisedef.kxfrac_min * wavenumber_max(grid.N.x());
#pragma warning(pop)
#pragma warning(push,disable:1572)
    if (noisedef.kxfrac_max != 1 || noisedef.kxfrac_min != 0) {
#pragma warning(pop)
        INFO0("Perturbations added only to absolute X wavenumbers in range ["
                << dkx_min << ":" << dkx_max << "]");
    }

    // Compute and display kzfrac_min, kzfrac_max constraints
#pragma warning(push,disable:2259)
    const int dkz_max = noisedef.kzfrac_max * wavenumber_max(grid.N.z());
    const int dkz_min = noisedef.kzfrac_min * wavenumber_max(grid.N.z());
#pragma warning(pop)
#pragma warning(push,disable:1572)
    if (noisedef.kzfrac_max != 1 || noisedef.kzfrac_min != 0) {
#pragma warning(pop)
    INFO0("Perturbations added only to absolute Z wavenumbers in range ["
            << dkz_min << ":" << dkz_max << "]");
    }

    // Form mass matrix to convert (wave, collocation point values, wave)
    // perturbations to (wave, coefficients, wave)
    bsplineop_luz massluz(bop);
    massluz.factor_mass(bop);

    // Set L'Ecuyer et al.'s RngStream seed.  Use a distinct Substream for each
    // wall-normal pencil to ensure process is a) repeatable despite changes in
    // processor count, b) easy to code, and c) embarrassingly parallel.
    RngStream rng;
    {
        boost::array<unsigned long,6> seed;
        std::fill(seed.begin(), seed.end(), noisedef.seed);
        rng.SetSeed(seed.data());
    }
    rng.IncreasedPrecis(true);  // Use more bits of resolution

    // Approach is the following:
    //  0) Allocate storage for state and three additional scalar fields.
    //  1) Generate a random vector-valued field \tilde{A}.
    //     \tilde{A}'s x- and z-derivatives have zero mean by periodicity;
    //  2) Zero first two B-spline coefficients near walls so partial_y
    //     \tilde{A} vanishes at the wall
    //  3) This step intentionally left blank.
    //  4) Compute curl A in physical space and rescale so maximum
    //     pointwise norm of A is maxfluct.  curl A is solenoidal
    //     and now has the velocity perturbation properties we desire.
    //  5) Store curl A in physical space in the three scalar fields.
    //  6) Copy state into auxiliary state storage and bring to
    //     physical space.
    //  7) At each point compute velocity and internal energy.  Perturb
    //     velocity and compute new total energy using perturbed
    //     velocities.
    //  8) Bring perturbed state information back to wavespace.
    //  9) Overwrite state storage with the new perturbed state.

    //  0) Allocate storage for state and three additional scalar fields.
    boost::scoped_ptr<contiguous_state<4,complex_t> > _s_ptr( // RAII
            support::allocate_padded_state<contiguous_state<4,complex_t> >(
                state_count + 3, dgrid));
    contiguous_state<4,complex_t> &s = *_s_ptr;               // Shorthand
    std::fill(s.range().begin(), s.range().end(), 0);        // Zero memory

    // 1) Generate a random vector-valued field \tilde{A}.
    // For each scalar component of \tilde{A}...
    for (size_t l = 0; l < 3; ++l) {

        for (int k = 0; k < grid.dN.z(); ++k) {
            if (!wavenumber_translatable(grid.N.z(), grid.dN.z(), k)) continue;


            for (int i = 0; i < grid.dN.x(); ++i) {
                if (!wavenumber_translatable(grid.N.x(), grid.dN.x(), i)) continue;

                // ...and advance RngStream to the (i, ., k) substream...
                // ...(necessary for processor-topology independence)...
                rng.ResetNextSubstream();

                // ...but only the rank holding the (i, ., k) pencil continues.
                if (   k <  dgrid.local_wave_start.z()
                    || k >= dgrid.local_wave_end.z()
                    || i <  dgrid.local_wave_start.x()
                    || i >= dgrid.local_wave_end.x()) continue;

                // Satisfy fluct_kzfrac_min, fluct_kzfrac_max constraints
                if (wavenumber_abs(grid.dN.z(), k) < dkz_min) continue;
                if (wavenumber_abs(grid.dN.z(), k) > dkz_max) continue;

                // Satisfy fluct_kxfrac_min, fluct_kxfrac_max constraints
                if (wavenumber_abs(grid.dN.x(), i) < dkx_min) continue;
                if (wavenumber_abs(grid.dN.x(), i) > dkx_max) continue;

                // Compute local indices for global (i, ., k).
                const int local_i = i - dgrid.local_wave_start.x();
                const int local_k = k - dgrid.local_wave_start.z();

                // ...generate coeffs with well-defined pseudorandom order.
                //  2) Zero first two B-spline coefficients near walls
                //     so partial_y \tilde{A} vanishes at the wall
                s[2*l][0][local_i][local_k] = 0;
                s[2*l][1][local_i][local_k] = 0;
                for (int j = 2; j < Ny - 2; ++j) {
                    const real_t magnitude = rng.RandU01();
                    const real_t phase     = rng.RandU01() * twopi;
                    s[2*l][j][local_i][local_k] = std::polar(magnitude, phase);
                }
                s[2*l][Ny - 2][local_i][local_k] = 0;
                s[2*l][Ny - 1][local_i][local_k] = 0;

            } // end X

        } // end Z

    } // end scalar components of A

    //  4) Compute curl A in physical space and rescale so maximum
    //     pointwise norm of A is maxfluct.  curl A is solenoidal
    //     and now has the velocity perturbation properties we desire.

    // Copy s[2l] into s[2l+1] as we need two copies to compute curl A
    for (size_t l = 0; l < 3; ++l) s[2*l+1] = s[2*l];

    // Prepare physical-space view of the wave-space storage
    support::physical_view<state_count+3>::type p
        = support::physical_view<state_count+3>::create(dgrid, s);

    // Initializing operator_base to access decomposition-ready utilities
    operator_base obase(grid, dgrid, b, bop);

    // From Ax in s[0] compute \partial_y Ax
    obase.bop_apply(1, 1.0, s, 0);
    dgrid.transform_wave_to_physical(&p.coeffRef(0,0));

    // From Ax in s[1]  compute \partial_z Ax
    obase.bop_apply(0, 1.0, s, 1);
    obase.diffwave_apply(0, 1, 1.0, s, 1);
    dgrid.transform_wave_to_physical(&p.coeffRef(1,0));

    // From Ay in s[2] compute \partial_x Ay
    obase.bop_apply(0, 1.0, s, 2);
    obase.diffwave_apply(1, 0, 1.0, s, 2);
    dgrid.transform_wave_to_physical(&p.coeffRef(2,0));

    // From Ay in s[3] compute \partial_z Ay
    obase.bop_apply(0, 1.0, s, 3);
    obase.diffwave_apply(0, 1, 1.0, s, 3);
    dgrid.transform_wave_to_physical(&p.coeffRef(3,0));

    // From Az in s[4] compute \partial_x Az
    obase.bop_apply(0, 1.0, s, 4);
    obase.diffwave_apply(1, 0, 1.0, s, 4);
    dgrid.transform_wave_to_physical(&p.coeffRef(4,0));

    // From Az in s[5] compute \partial_y Az
    obase.bop_apply(1, 1.0, s, 5);
    dgrid.transform_wave_to_physical(&p.coeffRef(5,0));

    // Store curl A in s[{5,6,7}] and find global maximum magnitude of curl A
    real_t maxmagsquared = 0;
    size_t offset = 0;
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Assert curl A is identically zero at the walls
                if (j == 0 || j == Ny - 1) {
#pragma warning(push,disable:1572)
                    assert(p(0, offset) == 0.0);
                    assert(p(1, offset) == 0.0);
                    assert(p(2, offset) == 0.0);
                    assert(p(3, offset) == 0.0);
                    assert(p(4, offset) == 0.0);
                    assert(p(5, offset) == 0.0);
#pragma warning(pop)
                }

                // The definition of curl gives the following components
                //   1) \partial_y A_z - \partial_z A_y
                //   2) \partial_z A_x - \partial_x A_z
                //   3) \partial_x A_y - \partial_y A_x
                // where the mean of \partial_x and \partial_z terms must be
                // zero by periodicity.  Components 1 and 3 may have nonzero
                // mean because they wall-normal derivatives contributions.
                const Vector3r curlA(p(5, offset) - p(3, offset),
                                     p(1, offset) - p(4, offset),
                                     p(2, offset) - p(0, offset));

                //  5) Store curl A in physical space in the 3 scalar fields.
                //
                // A nonzero mean in the x, y, and z directions is,
                // respectively, "corrected" by bulk forcing, the wall, and
                // viscous effects.  The spanwise viscous effects presumably
                // have the slowest timescale so rotate the components from 123
                // to 312 to reduce the simulation time before stationarity.
                // This rotation may introduce acoustic noise.
                p(state_count + 0, offset) = curlA.z();
                p(state_count + 1, offset) = curlA.x();
                p(state_count + 2, offset) = curlA.y();

                maxmagsquared = math::maxnan(
                        maxmagsquared, curlA.squaredNorm());

            } // end X

        } // end Z

    } // end Y
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, &maxmagsquared, 1,
                mpi::datatype<real_t>::value,
                MPI_MAX, MPI_COMM_WORLD));

    // Rescale curl A components so max ||curl A|| == maxfluct
    p.row(state_count + 0) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(state_count + 1) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(state_count + 2) *= (maxfluct / std::sqrt(maxmagsquared));

    //  6) Copy state into auxiliary state storage and bring to
    //     physical space.
    for (size_t i = 0; i < state_count; ++i) {
        s[i] = state[i];
        obase.bop_apply(0, 1.0, s, i);
        dgrid.transform_wave_to_physical(&p.coeffRef(i,0));
    }

    //  7) At each point compute velocity and internal energy.  Perturb
    //     velocity and compute new total energy using perturbed
    //     velocities.
    const real_t Ma = scenario.Ma;
    offset = 0;
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Retrieve internal energy
                real_t         e(p(ndx::e,   offset));
                Vector3r       m(p(ndx::mx,  offset),
                                 p(ndx::my,  offset),
                                 p(ndx::mz,  offset));
                const real_t rho(p(ndx::rho, offset));
                const real_t e_int = rholut::energy_internal(Ma, rho, m, e);

                // Perturb momentum and compute updated total energy
                m.x() += rho * p(state_count + 0, offset);
                m.y() += rho * p(state_count + 1, offset);
                m.z() += rho * p(state_count + 2, offset);
                const real_t e_kin = rholut::energy_kinetic(Ma, rho, m);
                e = e_int + e_kin;

                // Store results back to state fields
                p(ndx::e,  offset) = e;
                p(ndx::mx, offset) = m.x();
                p(ndx::my, offset) = m.y();
                p(ndx::mz, offset) = m.z();

            } // end X

        } // end Z

    } // end Y

    //  8) Bring perturbed state information back to wavespace (no rho!)
    // Build FFT normalization constant into Y direction's mass matrix.
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, bop);
    massluz.factor();
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::e , 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::mx, 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::my, 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::mz, 0));  // X, Z
    obase.bop_solve(massluz, s, ndx::e );                       // Y
    obase.bop_solve(massluz, s, ndx::mx);                       // Y
    obase.bop_solve(massluz, s, ndx::my);                       // Y
    obase.bop_solve(massluz, s, ndx::mz);                       // Y

    //  9) Overwrite state storage with the new perturbed state (not rho!)
    state[ndx::e ] = s[ndx::e ];
    state[ndx::mx] = s[ndx::mx];
    state[ndx::my] = s[ndx::my];
    state[ndx::mz] = s[ndx::mz];
}

void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        contiguous_state<4,complex_t> &swave,
        const problem::GridDefinition &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &bop,
        const real_t simulation_time)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Initialize operator_base to access decomposition-ready utilities
    operator_base obase(grid, dgrid, b, bop);

    // Allocate one field of temporary storage for scratch purposes
    boost::scoped_ptr<contiguous_state<4,complex_t> > _scratch_ptr( // RAII
        support::allocate_padded_state<contiguous_state<4,complex_t> >(1,dgrid));
    contiguous_state<4,complex_t> &scratch = *_scratch_ptr;         // Shorthand
    multi_array::fill(scratch, 0);                                 // Defensive

    // Prepare physical-space view of the wave-space scratch storage
    support::physical_view<1>::type phys
            = support::physical_view<1>::create(dgrid, scratch);

    // Prepare factored mass matrix for repeated use
    bsplineop_luz massluz(bop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, bop);
    massluz.factor();

    // For each scalar field...
    for (size_t f = 0; f < state_count; ++f) {

        // ...compute the manufactured solution in physical space...
        size_t offset = 0;
        for (int j = dgrid.local_physical_start.y();
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

// This looks like logic from explicit_op.cpp but does not belong there.
// Reading through that file, especially the NonlinearOperator::apply_operator
// implementation, is recommended before reviewing this logic.  This routine is
// definitely suboptimal but is expected to be invoked very infrequently and
// therefore to not be a prime target for optimization.
mean sample_mean_quantities(
        const ScenarioDefinition &scenario,
        const problem::GridDefinition &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &bop,
        contiguous_state<4,complex_t> &swave,
        const real_t t)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Shorthand
    const size_t Ny = swave.shape()[1];
    namespace acc = boost::accumulators;
    typedef contiguous_state<4,complex_t> state_type;

    // State enters method as coefficients in X, Y, and Z directions

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.
    struct aux { enum {
        e_y,   e_x,   e_z,
        mx_y,  mx_x,  mx_z,
        my_y,  my_x,  my_z,
        mz_y,  mz_x,  mz_z,
        rho_y, rho_x, rho_z,
        count // Sentry
    }; };

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    boost::scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux::count, dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                 // Shorthand

    // Ensure state storage meets this routine's assumptions

    // Sanity check incoming swave's and auxw's shape and contiguity
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    SUZERAIN_ENSURE(std::equal(swave.shape() + 1, swave.shape() + 4,
                               auxw.shape() + 1));
    SUZERAIN_ENSURE(std::equal(swave.strides() + 1, swave.strides() + 4,
                               auxw.strides() + 1));

    // Rank-specific details accumulated in ret to be MPI_Reduce-d later
    mean ret(t, Ny);

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret.rho_E()        = Map<VectorXc>(swave[ndx::e  ].origin(), Ny).real();
        ret.rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret.rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret.rho_u().col(2) = Map<VectorXc>(swave[ndx::mz ].origin(), Ny).real();
        ret.rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Obtain access to helper routines for differentiation
    operator_base obase(grid, dgrid, b, bop);

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    obase.zero_dealiasing_modes(   swave, ndx::e);
    obase.bop_accumulate(1,    1., swave, ndx::e, 0., auxw, aux::e_y);
    obase.bop_apply     (0,    1., swave, ndx::e);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    obase.diffwave_accumulate(1, 0, 1., swave, ndx::e, 0., auxw, aux::e_x);
    obase.diffwave_accumulate(0, 1, 1., swave, ndx::e, 0., auxw, aux::e_z);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    obase.zero_dealiasing_modes(   swave, ndx::mx);
    obase.bop_accumulate(1,    1., swave, ndx::mx, 0., auxw, aux::mx_y);
    obase.bop_apply     (0,    1., swave, ndx::mx);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    obase.diffwave_accumulate(1, 0, 1., swave, ndx::mx, 0., auxw, aux::mx_x);
    obase.diffwave_accumulate(0, 1, 1., swave, ndx::mx, 0., auxw, aux::mx_z);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    obase.zero_dealiasing_modes(   swave, ndx::my);
    obase.bop_accumulate(1,    1., swave, ndx::my, 0., auxw, aux::my_y);
    obase.bop_apply     (0,    1., swave, ndx::my);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    obase.diffwave_accumulate(1, 0, 1., swave, ndx::my, 0., auxw, aux::my_x);
    obase.diffwave_accumulate(0, 1, 1., swave, ndx::my, 0., auxw, aux::my_z);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    obase.zero_dealiasing_modes(   swave, ndx::mz);
    obase.bop_accumulate(1,    1., swave, ndx::mz, 0., auxw, aux::mz_y);
    obase.bop_apply     (0,    1., swave, ndx::mz);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    obase.diffwave_accumulate(1, 0, 1., swave, ndx::mz, 0., auxw, aux::mz_x);
    obase.diffwave_accumulate(0, 1, 1., swave, ndx::mz, 0., auxw, aux::mz_z);

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    obase.zero_dealiasing_modes(   swave, ndx::rho);
    obase.bop_accumulate(1,    1., swave, ndx::rho, 0., auxw, aux::rho_y);
    obase.bop_apply     (0,    1., swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    obase.diffwave_accumulate(1, 0, 1., swave, ndx::rho,  0., auxw, aux::rho_x);
    obase.diffwave_accumulate(0, 1, 1., swave, ndx::rho,  0., auxw, aux::rho_z);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    support::physical_view<aux::count>::type auxp
        = support::physical_view<aux::count>::create(dgrid, auxw);
    support::physical_view<state_count>::type sphys
        = support::physical_view<state_count>::create(dgrid, swave);
    for (size_t i = 0; i < state_count; ++i) {
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    size_t offset = 0;
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Type used to accumulate mean quantities versus wall-normal position
        typedef acc::accumulator_set<
                real_t,
#if BOOST_VERSION >= 104700
                acc::stats<acc::tag::sum_kahan>
#else
                acc::stats<acc::tag::sum>
#endif
            > accumulator_type;

        // Accumulators for each mean quantity computed in physical space.
        // For example, quantity "foo" has accumulator "sum_foo".
        // Declared within Y loop so they are reset on each Y iteration.
#define DECLARE(r, data, tuple)                                              \
        accumulator_type BOOST_PP_CAT(sum_,BOOST_PP_TUPLE_ELEM(2, 0, tuple)) \
                [BOOST_PP_TUPLE_ELEM(2, 1, tuple)];
        BOOST_PP_SEQ_FOR_EACH(DECLARE,,CHANNEL_MEAN_PHYSICAL)
#undef DECLARE

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Unpack total energy-related quantities
                const real_t e(sphys(ndx::e, offset));
                const Vector3r grad_e(auxp(aux::e_x, offset),
                                      auxp(aux::e_y, offset),
                                      auxp(aux::e_z, offset));

                // Unpack momentum-related quantities
                const Vector3r m(sphys(ndx::mx, offset),
                                 sphys(ndx::my, offset),
                                 sphys(ndx::mz, offset));
                const real_t div_m = auxp(aux::mx_x, offset)
                                   + auxp(aux::my_y, offset)
                                   + auxp(aux::mz_z, offset);
                const Matrix3r grad_m;
                const_cast<Matrix3r&>(grad_m) <<
                                     auxp(aux::mx_x,  offset),
                                     auxp(aux::mx_y,  offset),
                                     auxp(aux::mx_z,  offset),
                                     auxp(aux::my_x,  offset),
                                     auxp(aux::my_y,  offset),
                                     auxp(aux::my_z,  offset),
                                     auxp(aux::mz_x,  offset),
                                     auxp(aux::mz_y,  offset),
                                     auxp(aux::mz_z,  offset);

                // Unpack density-related quantities
                const real_t rho(sphys(ndx::rho, offset));
                const Vector3r grad_rho(auxp(aux::rho_x, offset),
                                        auxp(aux::rho_y, offset),
                                        auxp(aux::rho_z, offset));

                // Compute local quantities based upon state.
                const Vector3r u   = rholut::u(
                                        rho, m);
                const real_t div_u = rholut::div_u(
                                        rho, grad_rho, m, div_m);
                const Matrix3r grad_u = rholut::grad_u(
                                        rho, grad_rho, m, grad_m);

                real_t p, T, mu, lambda;
                Vector3r grad_p, grad_T, grad_mu, grad_lambda;
                rholut::p_T_mu_lambda(
                    scenario.alpha, scenario.beta, scenario.gamma, scenario.Ma,
                    rho, grad_rho, m, grad_m, e, grad_e,
                    p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

                const Matrix3r tau = rholut::tau(
                                        mu, lambda, div_u, grad_u);
                const Vector3r tau_u = tau * u;

                // Accumulate quantities into sum_XXX using function syntax.
                sum_E[0](e / rho);

                sum_T[0](T);

                sum_mu[0](mu);

                sum_nu[0](mu / rho);

                sum_u[0](u.x());
                sum_u[1](u.y());
                sum_u[2](u.z());

                sum_sym_grad_u[0]( grad_u(0,0)                   );
                sum_sym_grad_u[1]((grad_u(0,1) + grad_u(1,0)) / 2);
                sum_sym_grad_u[2]((grad_u(0,2) + grad_u(2,0)) / 2);
                sum_sym_grad_u[3]( grad_u(1,1)                   );
                sum_sym_grad_u[4]((grad_u(1,2) + grad_u(2,1)) / 2);
                sum_sym_grad_u[5]( grad_u(2,2)                   );

                sum_sym_rho_grad_u[0](rho *  grad_u(0,0)                   );
                sum_sym_rho_grad_u[1](rho * (grad_u(0,1) + grad_u(1,0)) / 2);
                sum_sym_rho_grad_u[2](rho * (grad_u(0,2) + grad_u(2,0)) / 2);
                sum_sym_rho_grad_u[3](rho *  grad_u(1,1)                   );
                sum_sym_rho_grad_u[4](rho * (grad_u(1,2) + grad_u(2,1)) / 2);
                sum_sym_rho_grad_u[5](rho *  grad_u(2,2)                   );

                sum_grad_T[0](grad_T.x());
                sum_grad_T[1](grad_T.y());
                sum_grad_T[2](grad_T.z());

                sum_rho_grad_T[0](rho * grad_T.x());
                sum_rho_grad_T[1](rho * grad_T.y());
                sum_rho_grad_T[2](rho * grad_T.z());

                sum_tau_colon_grad_u[0]((tau.transpose()*grad_u).trace());

                sum_tau[0](tau(0,0));
                sum_tau[1](tau(0,1));
                sum_tau[2](tau(0,2));
                sum_tau[3](tau(1,1));
                sum_tau[4](tau(1,2));
                sum_tau[5](tau(2,2));

                sum_tau_u[0](tau_u.x());
                sum_tau_u[1](tau_u.y());
                sum_tau_u[2](tau_u.z());

                sum_p_div_u[0](p*div_u);

                sum_rho_u_u[0](rho * u.x() * u.x());
                sum_rho_u_u[1](rho * u.x() * u.y());
                sum_rho_u_u[2](rho * u.x() * u.z());
                sum_rho_u_u[3](rho * u.y() * u.y());
                sum_rho_u_u[4](rho * u.y() * u.z());
                sum_rho_u_u[5](rho * u.z() * u.z());

                sum_rho_u_u_u[0](rho * u.x() * u.x() * u.x());
                sum_rho_u_u_u[1](rho * u.x() * u.x() * u.y());
                sum_rho_u_u_u[2](rho * u.x() * u.x() * u.z());
                sum_rho_u_u_u[3](rho * u.x() * u.y() * u.y());
                sum_rho_u_u_u[4](rho * u.x() * u.y() * u.z());
                sum_rho_u_u_u[5](rho * u.x() * u.z() * u.z());
                sum_rho_u_u_u[6](rho * u.y() * u.y() * u.y());
                sum_rho_u_u_u[7](rho * u.y() * u.y() * u.z());
                sum_rho_u_u_u[8](rho * u.y() * u.z() * u.z());
                sum_rho_u_u_u[9](rho * u.z() * u.z() * u.z());

                sum_rho_T_u[0](rho * T * u.x());
                sum_rho_T_u[1](rho * T * u.y());
                sum_rho_T_u[2](rho * T * u.z());

                sum_rho_mu[0](rho * mu);

                sum_mu_S[0](mu * ( grad_u(0,0)                    - div_u / 3));
                sum_mu_S[1](mu * ((grad_u(0,1) + grad_u(1,0)) / 2            ));
                sum_mu_S[2](mu * ((grad_u(0,2) + grad_u(2,0)) / 2            ));
                sum_mu_S[3](mu * ( grad_u(1,1)                    - div_u / 3));
                sum_mu_S[4](mu * ((grad_u(1,2) + grad_u(2,1)) / 2            ));
                sum_mu_S[5](mu * ( grad_u(2,2)                    - div_u / 3));

                sum_mu_div_u[0](mu * div_u);

                sum_mu_grad_T[0](mu * grad_T.x());
                sum_mu_grad_T[1](mu * grad_T.y());
                sum_mu_grad_T[2](mu * grad_T.z());

                // TODO Sum mean slow growth forcing contributions (Redmine #2496)

                sum_SrhoE[0](0);

                sum_Srhou[0](0);
                sum_Srhou[1](0);
                sum_Srhou[2](0);

                sum_Srho[0](0);

                sum_Srhou_dot_u[0](0);

            } // end X

        } // end Z

        // Move y-specific sums into MPI-reduction-ready storage for y(j) using
        // Eigen comma initialization syntax.  Yes, this is not stride 1.
#define EXTRACT_SUM(z, n, q) acc::sum(sum_##q[n])
#define MOVE_SUM_INTO_TMP(r, data, tuple)                                 \
        ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().row(j) <<                  \
            BOOST_PP_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, tuple),               \
                          EXTRACT_SUM, BOOST_PP_TUPLE_ELEM(2, 0, tuple));
        BOOST_PP_SEQ_FOR_EACH(MOVE_SUM_INTO_TMP,,CHANNEL_MEAN_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // quantities while other ranks have zeros in those locations.

    // Reduce sums onto rank zero and then return garbage from non-zero ranks
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {

        // Reduce operation requires no additional storage on rank-zero
        SUZERAIN_MPICHKR(MPI_Reduce(
                MPI_IN_PLACE, ret.storage.data(), ret.storage.size(),
                mpi::datatype<mean::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

    } else {

        // Reduce operation requires temporary storage on non-zero ranks
        ArrayXXr tmp;
        tmp.resizeLike(ret.storage);
        tmp.setZero();
        SUZERAIN_MPICHKR(MPI_Reduce(ret.storage.data(), tmp.data(), tmp.size(),
                mpi::datatype<mean::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

        // Force non-zero ranks contain all NaNs to help detect usage errors
        ret.storage.fill(
                std::numeric_limits<mean::storage_type::Scalar>::quiet_NaN());

        // Return from all non-zero ranks
        return ret;

    }

    // Only rank zero reaches this logic because of return statement just above
    assert(mpi::comm_rank(MPI_COMM_WORLD) == 0);

    // Physical space sums, which are at collocation points, need to be
    // divided by the dealiased extents and converted to coefficients.
    const real_t scale_factor = dgrid.global_physical_extent.x()
                              * dgrid.global_physical_extent.z();
    bsplineop_lu scaled_mass(bop);
    scaled_mass.opform(1, &scale_factor, bop);
    scaled_mass.factor();
    scaled_mass.solve(mean::nscalars::physical,
            ret.storage.middleCols<mean::nscalars::physical>(
                mean::nscalars::wave).data(),
            ret.storage.innerStride(), ret.storage.outerStride());

    // Fill with NaNs those samples which were not computed by this method
#define FILL(r, data, tuple)                                               \
    ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().fill(                           \
            std::numeric_limits<mean::storage_type::Scalar>::quiet_NaN());
    BOOST_PP_SEQ_FOR_EACH(FILL,,CHANNEL_MEAN_IMPLICIT)
#undef FILL

    return ret;
}

// Helper for the store(...) implementation just below
class mean_storer
{

public:

    mean_storer(const esio_handle esioh, const std::string& prefix)
        : esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    void operator()(const std::string& name, const EigenArray& dat) const
    {
        int procid;
        esio_handle_comm_rank(esioh, &procid);
        esio_field_establish(esioh,
                             1,          0, (procid == 0 ? 1          : 0),
                             dat.cols(), 0, (procid == 0 ? dat.cols() : 0),
                             dat.rows(), 0, (procid == 0 ? dat.rows() : 0));

        std::string key(prefix);
        key.append(name);
        esio_field_write(esioh, key.c_str(), dat.data(), 0, 0, 0,
            "Mean quantity sample stored using row-major indices (B-spline"
            " coefficient, tensor component, sample number) where the"
            " B-spline basis is defined by /Ny, /breakpoints_y, and /knots");
    }

private:

    esio_handle esioh;
    std::string prefix;

};

void store(const esio_handle h, const mean& m)
{
    mean_storer f(h, "bar_");
    m.foreach(f);
}

// Helper for the load(...) implementation just below
class mean_loader
{

public:

    mean_loader(const esio_handle esioh, const std::string& prefix)
        : esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    void operator()(const std::string& name, const EigenArray& dat_) const {

        // http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
        EigenArray& dat = const_cast<EigenArray&>(dat_);

        int procid;
        esio_handle_comm_rank(esioh, &procid);

        std::string key(prefix);
        key.append(name);

        int cglobal, bglobal, aglobal;
        if (ESIO_SUCCESS == esio_field_size(esioh, key.c_str(),
                                            &cglobal, &bglobal, &aglobal)
            && (   EigenArray::ColsAtCompileTime == Eigen::Dynamic
                || dat.cols() == bglobal)) {
            dat.resize(aglobal, bglobal);
            esio_field_establish(esioh,
                                 1,          0, (procid == 0 ? 1          : 0),
                                 dat.cols(), 0, (procid == 0 ? dat.cols() : 0),
                                 dat.rows(), 0, (procid == 0 ? dat.rows() : 0));
            esio_field_read(esioh, key.c_str(), dat.data(), 0, 0, 0);
        } else {
            WARN0("Unable to load " << key << " for nscalar = " << dat.cols());
            dat.fill(std::numeric_limits<
                    typename EigenArray::Scalar>::quiet_NaN());
        }
    }

private:

    esio_handle esioh;
    std::string prefix;

};

void load(const esio_handle h, mean& m)
{
    int cglobal, bglobal, aglobal;
    if (ESIO_SUCCESS == esio_field_size(h, "bar_rho",
                                        &cglobal, &bglobal, &aglobal)) {
        m.storage.resize(aglobal, Eigen::NoChange);
        mean_loader f(h, "bar_");
        m.foreach(f);
        // m.t presumably set externally using load_time
    } else {
        WARN0("No mean samples loaded-- unable to anticipate storage needs");
        m.storage.fill(std::numeric_limits<real_t>::quiet_NaN());
        m.t = std::numeric_limits<real_t>::quiet_NaN();
    }
}

} /* namespace perfect */ } /* namespace suzerain */
