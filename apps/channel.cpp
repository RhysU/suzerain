/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * channel.cpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <suzerain/error.h>
#include <suzerain/htstretch.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/RngStream.hpp>

#include "logger.hpp"
#include "channel.hpp"

// Manufactured solution classes explicitly instantiated for debugging
template class nsctpl_rholut::manufactured_solution<real_t>;

using boost::numeric_cast;

namespace channel {

const boost::array<const char *,field::count> field::name = {{
    "rho", "rhou", "rhov", "rhow", "rhoe"
}};

void mpi_abort_on_error_handler_gsl(const char * reason,
                                    const char * file,
                                    int line,
                                    int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "GSL", gsl_strerror(error_code));
}

void mpi_abort_on_error_handler_suzerain(const char * reason,
                                         const char * file,
                                         int line,
                                         int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "Suzerain", suzerain_strerror(error_code));
}

void mpi_abort_on_error_handler_esio(const char * reason,
                                     const char * file,
                                     int line,
                                     int error_code)
{
    return mpi_abort_on_error_handler(reason, file, line,
            error_code, "ESIO", esio_strerror(error_code));
}

void mpi_abort_on_error_handler(const char * reason,
                                const char * file,
                                int line,
                                int error_code,
                                const char * origin,
                                const char * strerror)
{
    FATAL((origin ? origin : "NULLORIGIN")
          << " reports '"
          << (reason ? reason : "NULLREASON")
          << "' as code #"
          << error_code
          << " ('"
          << (strerror ? strerror : "NULLSTRERROR")
          << "') from "
          << (file ? file : "NULLFILE")
          << ':'
          << line);
    MPI_Abort(MPI_COMM_WORLD, errno ? errno : EXIT_FAILURE);
}

void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario)
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

    esio_line_write(h, "bulk_rhou", &scenario.bulk_rhou, 0,
            scenario.options().find("bulk_rhou",false).description().c_str());

    esio_line_write(h, "alpha", &scenario.alpha, 0,
            scenario.options().find("alpha",false).description().c_str());

    esio_line_write(h, "beta", &scenario.beta, 0,
            scenario.options().find("beta",false).description().c_str());

    esio_line_write(h, "gamma", &scenario.gamma, 0,
            scenario.options().find("gamma",false).description().c_str());

    esio_line_write(h, "Lx", &scenario.Lx, 0,
            scenario.options().find("Lx",false).description().c_str());

    esio_line_write(h, "Ly", &scenario.Ly, 0,
            scenario.options().find("Ly",false).description().c_str());

    esio_line_write(h, "Lz", &scenario.Lz, 0,
            scenario.options().find("Lz",false).description().c_str());
}

void load(const esio_handle h,
          suzerain::problem::ScenarioDefinition<real_t>& scenario)
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

    if (!(boost::math::isnan)(scenario.bulk_rhou)) {
        INFO0("Overriding scenario using bulk_rhou = " << scenario.bulk_rhou);
    } else {
        esio_line_read(h, "bulk_rhou", &scenario.bulk_rhou, 0);
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

    if (!(boost::math::isnan)(scenario.Lx)) {
        INFO0("Overriding scenario using Lx = " << scenario.Lx);
    } else {
        esio_line_read(h, "Lx", &scenario.Lx, 0);
    }

    if (!(boost::math::isnan)(scenario.Ly)) {
        INFO0("Overriding scenario using Ly = " << scenario.Ly);
    } else {
        esio_line_read(h, "Ly", &scenario.Ly, 0);
    }

    if (!(boost::math::isnan)(scenario.Lz)) {
        INFO0("Overriding scenario using Lz = " << scenario.Lz);
    } else {
        esio_line_read(h, "Lz", &scenario.Lz, 0);
    }
}

void store(const esio_handle h,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing GridDefinition parameters");

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Nx", &grid.N.x(), 0,
            grid.options().find("Nx",false).description().c_str());

    esio_line_write(h, "DAFx", &grid.DAF.x(), 0,
            grid.options().find("DAFx",false).description().c_str());

    esio_line_write(h, "Ny", &grid.N.y(), 0,
            grid.options().find("Ny",false).description().c_str());

    esio_line_write(h, "k", &grid.k, 0,
            grid.options().find("k",false).description().c_str());

    esio_line_write(h, "htdelta", &grid.htdelta, 0,
            grid.options().find("htdelta",false).description().c_str());

    esio_line_write(h, "Nz", &grid.N.z(), 0,
            grid.options().find("Nz",false).description().c_str());

    esio_line_write(h, "DAFz", &grid.DAF.z(), 0,
            grid.options().find("DAFz",false).description().c_str());

    DEBUG0("Storing wavenumber vectors for Fourier bases");

    Eigen::ArrayXc buf(std::max(grid.N.x(), grid.N.z()));

    // Obtain wavenumbers via computing 1*(i*kx)/i
    buf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.data(),
            Lx, Lz, 1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(h, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(buf.data()),
            2, "Wavenumbers in streamwise X direction"); // Re(buf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    buf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.data(),
            Lx, Lz, 1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(h, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(buf.data()),
            2, "Wavenumbers in spanwise Z direction"); // Re(buf)
}

void load(const esio_handle h,
          suzerain::problem::GridDefinition& grid)
{
    DEBUG0("Loading GridDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (grid.N.x()) {
        INFO0("Overriding grid using Nx = " << grid.N.x());
    } else {
        int value;
        esio_line_read(h, "Nx", &value, 0);
        grid.Nx(value);
    }

    if (!((boost::math::isnan)(grid.DAF.x()))) {
        INFO0("Overriding grid using DAFx = " << grid.DAF.x());
    } else {
        double factor;
        esio_line_read(h, "DAFx", &factor, 0);
        grid.DAFx(factor);
    }

    if (grid.N.y()) {
        INFO0("Overriding grid using Ny = " << grid.N.y());
    } else {
        int value;
        esio_line_read(h, "Ny", &value, 0);
        grid.Ny(value);
    }

    if (grid.k) {
        INFO0("Overriding grid using k = " << grid.k);
    } else {
        esio_line_read(h, "k", &grid.k, 0);
    }

    if (!((boost::math::isnan)(grid.htdelta))) {
        INFO0("Overriding grid using htdelta = " << grid.htdelta);
    } else {
        esio_line_read(h, "htdelta", &grid.htdelta, 0);
    }

    if (grid.N.z()) {
        INFO0("Overriding grid using Nz = " << grid.N.z());
    } else {
        int value;
        esio_line_read(h, "Nz", &value, 0);
        grid.Nz(value);
    }

    if (!((boost::math::isnan)(grid.DAF.z()))) {
        INFO0("Overriding grid using DAFz = " << grid.DAF.z());
    } else {
        double factor;
        esio_line_read(h, "DAFz", &factor, 0);
        grid.DAFz(factor);
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
           const suzerain::problem::ScenarioDefinition<real_t>& scenario,
           const boost::shared_ptr<manufactured_solution>& msoln)
{
    static const char location[] = "channel::manufactured_solution";

    // Always write a flag to indicate whether or not a MS is in use
    const int in_use = msoln ? 1 : 0;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &in_use, 0,
            "Is a channel::manufactured_solution in use?");

    // Only proceed if an MS is in use
    if (!in_use) return;

    DEBUG0("Storing channel::manufactured_solution parameters");

    // Check parameters stored with the scenario not the manufactured solution
#pragma warning(push,disable:1572)
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
    if (msoln->Lx    != scenario.Lx)
        WARN0("Manufactured solution Lx mismatches with scenario!");
    if (msoln->Ly    != scenario.Ly)
        WARN0("Manufactured solution Ly mismatches with scenario!");
    if (msoln->Lz    != scenario.Lz)
        WARN0("Manufactured solution Lz mismatches with scenario!");
#pragma warning(pop)

    // Non-scenario solution parameters are stored as attributes under location
    // Scenario parameters should be taken from ScenarioDefinition
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_storer, h, location, _1, _2));
    msoln->u.foreach_parameter(  bind(attribute_storer, h, location, _1, _2));
    msoln->v.foreach_parameter(  bind(attribute_storer, h, location, _1, _2));
    msoln->w.foreach_parameter(  bind(attribute_storer, h, location, _1, _2));
    msoln->T.foreach_parameter(  bind(attribute_storer, h, location, _1, _2));
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
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          boost::shared_ptr<manufactured_solution>& msoln)
{
    static const char location[] = "channel::manufactured_solution";


    // Determine if a manufactured solution should be loaded
    int in_use;
    esio_line_establish(h, 1, 0, 1); // All ranks load
    esio_line_read(h, location, &in_use, 0);

    // Only proceed if an MS is in use
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
    msoln->Lx    = scenario.Lx;
    msoln->Ly    = scenario.Ly;
    msoln->Lz    = scenario.Lz;

    // Non-scenario solution parameters are stored as attributes under location
    using boost::bind;
    msoln->rho.foreach_parameter(bind(attribute_loader, h, location, _1, _2));
    msoln->u.foreach_parameter(  bind(attribute_loader, h, location, _1, _2));
    msoln->v.foreach_parameter(  bind(attribute_loader, h, location, _1, _2));
    msoln->w.foreach_parameter(  bind(attribute_loader, h, location, _1, _2));
    msoln->T.foreach_parameter(  bind(attribute_loader, h, location, _1, _2));
}

void create(const int ndof,
            const int k,
            const double left,
            const double right,
            const double htdelta,
            boost::shared_ptr<suzerain::bspline>& b,
            boost::shared_ptr<suzerain::bsplineop>& bop)
{
    // Compute breakpoint locations
    Eigen::ArrayXd breakpoints(ndof + 2 - k);
    suzerain::math::linspace(0.0, 1.0, breakpoints.size(), breakpoints.data());
    for (int i = 0; i < breakpoints.size(); ++i) {
        breakpoints[i] = suzerain_htstretch2(htdelta, 1.0, breakpoints[i]);
    }
    breakpoints = (right - left) * breakpoints + left;

    // Generate the B-spline workspace based on order and breakpoints
    // Maximum non-trivial derivative operators included
    b = boost::make_shared<suzerain::bspline>(
            k, breakpoints.size(), breakpoints.data());
    assert(b->n() == ndof);
    bop.reset(new suzerain::bsplineop(
                *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
    assert(bop->n() == ndof);
}

void store(const esio_handle h,
           const boost::shared_ptr<suzerain::bspline>& b,
           const boost::shared_ptr<suzerain::bsplineop>& bop,
           const boost::shared_ptr<suzerain::bsplineop>& gop)
{
    // Ensure we were handed the appropriate discrete operators
    assert(bop->get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    assert(gop->get()->method == SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing B-spline knot details");

    Eigen::ArrayXr buf(b->nknot());

    for (int i = 0; i < b->nknot(); ++i) buf[i] = b->knot(i);
    esio_line_establish(h, b->nknot(), 0, (procid == 0 ? b->nknot() : 0));
    esio_line_write(h, "knots", buf.data(), 0,
            "Knots used to build B-spline basis");

    for (int i = 0; i < b->nbreak(); ++i) buf[i] = b->breakpoint(i);
    esio_line_establish(h, b->nbreak(), 0, (procid == 0 ? b->nbreak() : 0));
    esio_line_write(h, "breakpoints", buf.data(), 0,
            "Breakpoint locations used to build B-spline basis");

    for (int i = 0; i < b->n(); ++i) buf[i] = b->collocation_point(i);
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "collocation_points", buf.data(), 0,
            "Collocation points used to build discrete operators");

    b->integration_coefficients(0, buf.data());
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "integration_weights", buf.data(), 0,
            "Integrate by dotting B-spline coefficients against weights");

    char name[8]      = {};
    char comment[127] = {};

    for (int k = 0; k <= bop->nderiv(); ++k) {
        snprintf(name, sizeof(name)/sizeof(name[0]), "Dy%d", k);
        snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                "Wall-normal derivative Dy%d(i,j) = D%d[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = bop->ku(k) + 1 + bop->kl(k);
        esio_plane_establish(h,
                bop->n(), 0, (procid == 0 ? bop->n() : 0),
                lda,      0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, bop->D(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", bop->kl(k));
        esio_attribute_write(h, name, "ku", bop->ku(k));
        esio_attribute_write(h, name, "m",  bop->n());
        esio_attribute_write(h, name, "n",  bop->n());
    }

    DEBUG0("Storing B-spline Galerkin L2 derivative operators");

    for (int k = 0; k <= gop->nderiv(); ++k) {
        snprintf(name, sizeof(name)/sizeof(name[0]), "Gy%d", k);
        snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                "Wall-normal Galerkin L2 Gy%d(i,j) = G%d[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = gop->ku(k) + 1 + gop->kl(k);
        esio_plane_establish(h,
                gop->n(), 0, (procid == 0 ? gop->n() : 0),
                lda,      0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, gop->D(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", gop->kl(k));
        esio_attribute_write(h, name, "ku", gop->ku(k));
        esio_attribute_write(h, name, "m",  gop->n());
        esio_attribute_write(h, name, "n",  gop->n());
    }
}

void load(const esio_handle h,
          boost::shared_ptr<suzerain::bspline>& b,
          boost::shared_ptr<suzerain::bsplineop>& bop)
{
    DEBUG0("Loading B-spline workspace based on order and breakpoints");

    // All ranks load B-spline order
    int k;
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, "k", &k, 0);

    // htdelta is ignored

    // All ranks load B-spline breakpoints
    int nbreak;
    esio_line_size(h, "breakpoints", &nbreak);
    esio_line_establish(h, nbreak, 0, nbreak);
    Eigen::ArrayXr buf(nbreak);
    esio_line_read(h, "breakpoints", buf.data(), 0);

    // Collocation points are ignored

    // Construct B-spline workspace
    b = boost::make_shared<suzerain::bspline>(k, nbreak, buf.data());
    bop.reset(new suzerain::bsplineop(
                *b, k-2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE));
}

void store_time(const esio_handle h,
                real_t time)
{
    // Root writes details
    int rank;
    esio_handle_comm_rank(h, &rank);
    esio_line_establish(h, 1, 0, (rank == 0) ? 1 : 0);

    esio_line_write(h, "t", &time, 0, "Simulation physical time");

    DEBUG0("Stored simulation time " << time);
}

void load_time(const esio_handle h,
               real_t &time)
{
    // All ranks read details
    esio_line_establish(h, 1, 0, 1);

    esio_line_read(h, "t", &time, 0);

    DEBUG0("Loaded simulation time " << time);
}

template<>
suzerain::NoninterleavedState<4,complex_t>* allocate_padded_state(
           const std::size_t howmany_fields,
           const suzerain::pencil_grid& dgrid)
{
    // Create instance with appropriate padding to allow P3DFFTification
    suzerain::NoninterleavedState<4,complex_t> * const retval =
        new suzerain::NoninterleavedState<4,complex_t>(
            suzerain::to_yxz(howmany_fields, dgrid.local_wave_extent),
            suzerain::prepend(dgrid.local_wave_storage(), suzerain::strides_cm(
                              suzerain::to_yxz(dgrid.local_wave_extent)))
        );

    // Clear to avoid lingering NaN issues
    suzerain::multi_array::fill(*retval, 0);

    return retval;
}

/** State descriptions used within ESIO-based restart files */
static const boost::array<const char *,5> field_descriptions = {{
    "Nondimensional density coefficients stored row-major ZXY",
    "Nondimensional X momentum coefficients stored row-major ZXY",
    "Nondimensional Y momentum coefficients stored row-major ZXY",
    "Nondimensional Z momentum coefficients stored row-major ZXY",
    "Nondimensional total energy coefficients stored row-major ZXY"
}};

void store(const esio_handle h,
           const suzerain::NoninterleavedState<4,complex_t> &state,
           const suzerain::problem::GridDefinition& grid,
           const suzerain::pencil_grid& dgrid)
{
    typedef suzerain::NoninterleavedState<4,complex_t> store_type;

    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(grid.N.x(),
                                            grid.dN.x(),
                                            dgrid.local_wave_start.x(),
                                            dgrid.local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(grid.N.z(),
                                            grid.dN.z(),
                                            dgrid.local_wave_start.z(),
                                            dgrid.local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // Save each scalar field in turn...
    for (size_t i = 0; i < field::count; ++i) {

        // Create a view of the state for just the i-th scalar
        // Not strictly necessary, but aids greatly in debugging
        boost::multi_array_types::index_range all;
        boost::const_array_view_gen<store_type,3>::type field
            = state[boost::indices[i][all][all][all]];

        // ...which requires two writes per field, once per Z range
        for (int j = 0; j < 2; ++j) {

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &field[0]
                            [mxb[0] - dgrid.local_wave_start.x()]
                            [mzb[j] - dgrid.local_wave_start.z()];
            }

            // Collectively establish size of read across all ranks
            esio_field_establish(h,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Perform collective write operation from state_linear
            complex_field_write(h, field::name[i], src,
                                field.strides()[2],
                                field.strides()[1],
                                field.strides()[0],
                                field_descriptions[i]);
        }
    }
}

void load(const esio_handle h,
          suzerain::NoninterleavedState<4,complex_t> &state,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          const suzerain::bspline& b,
          const suzerain::bsplineop& bop)
{
    typedef suzerain::NoninterleavedState<4,complex_t> load_type;

    // Ensure local state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.global_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Obtain details on the restart field's global sizes
    int Fz, Fx, Fy, ncomponents;
    esio_field_sizev(h, field::name[0], &Fz, &Fx, &Fy, &ncomponents);
    assert(ncomponents == 2);

    // Prepare a file-specific B-spline basis
    boost::shared_ptr<suzerain::bspline> Fb;
    boost::shared_ptr<suzerain::bsplineop> Fbop;
    load(h, Fb, Fbop);
    assert(Fy == Fb->n());

    // Check if the B-spline basis in the file differs from ours.  Use strict
    // equality as minor knot differences magnify once collocation points
    // and operators are computed.
    bool bsplines_same =    b.k()     == Fb->k()
                         && b.n()     == Fb->n()
                         && b.nknot() == Fb->nknot();
    for (int j = 0; bsplines_same && j < b.nknot(); ++j) {
#pragma warning(push,disable:1572)
        bsplines_same = (b.knot(j) == Fb->knot(j));
#pragma warning(pop)
    }

    // Compute wavenumber translation logistics for X direction.
    // Requires turning a C2R FFT complex-valued coefficient count into a
    // real-valued coefficient count.  Further, need to preserve even- or
    // odd-ness of the coefficient count to handle, for example, Fx == 1.
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(2 * (Fx - 1) + (Fx & 1),
                                            grid.dN.x(),
                                            dgrid.local_wave_start.x(),
                                            dgrid.local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(Fz,
                                            grid.dN.z(),
                                            dgrid.local_wave_start.z(),
                                            dgrid.local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // Possibly prepare a temporary buffer into which to read each scalar field
    // and a factorization of b's mass matrix.  Used only when
    // !bsplines_same.
    typedef boost::multi_array<
        complex_t, 3, suzerain::blas::allocator<complex_t>::type
    > tmp_type;
    boost::scoped_ptr<tmp_type> tmp;
    boost::scoped_ptr<suzerain::bsplineop_luz> mass;
    if (!bsplines_same) {
        DEBUG0("Differences in B-spline basis require restart projection");
        const boost::array<tmp_type::index,3> extent = {{
            Fy, state.shape()[2], state.shape()[3]
        }};
        tmp.reset(new tmp_type(extent, boost::fortran_storage_order()));
        mass.reset(new suzerain::bsplineop_luz(bop));
        mass->form_mass(bop);
    }

    DEBUG0("Started loading simulation fields");

    // Load each scalar field in turn
    for (size_t i = 0; i < field::count; ++i) {

        // Create a view of the state for just the i-th scalar
        boost::multi_array_types::index_range all;
        boost::array_view_gen<load_type, 3>::type field
                = state[boost::indices[i][all][all][all]];

        // Clear storage prior to load to zero not-loaded coefficients
        if (bsplines_same) {
            suzerain::multi_array::fill(field, 0);
        } else {
            suzerain::multi_array::fill(*tmp, 0);
        }

        // Two ESIO read operations per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Collectively establish size of read across all ranks
            esio_field_establish(h, Fz, fzb[j], (fze[j] - fzb[j]),
                                    Fx, fxb[0], (fxe[0] - fxb[0]),
                                    Fy,      0, (            Fy));


            // Destination of read is NULL for empty READ operations.
            // Otherwise find starting point for coefficient loading.
            complex_t * dst = NULL;
            boost::array<load_type::index,3> dst_strides = {{ 0, 0, 0 }};
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                const boost::array<load_type::index,3> index_list = {{
                        0,
                        mxb[0] - dgrid.local_wave_start.x(),
                        mzb[j] - dgrid.local_wave_start.z()
                }};
                if (bsplines_same) {
                    dst = &field(index_list);
                    std::copy(field.strides(), field.strides() + 3,
                              dst_strides.begin());
                } else {
                    dst = &(*tmp)(index_list);
                    std::copy(tmp->strides(), tmp->strides() + 3,
                              dst_strides.begin());
                }
            }

            // Perform collective read operation into dst
            complex_field_read(h, field::name[i], dst,
                               dst_strides[2], dst_strides[1], dst_strides[0]);
        }


        // If necessary, interpolate between B-spline bases.
        // Relies heavily on both bases being collocation-based.
        // This will change dramatically if we ever go the L_2 route.
        if (!bsplines_same) {

            // Step 0: Obtain collocation points for new basis
            Eigen::ArrayXd points(b.n());
            for (int k = 0; k < b.n(); ++k) points[k] = b.collocation_point(k);

            // Step 1: Affine transformation of points into old basis domain
            // Allows stretching of old range onto new range.
            const double newmin = b.collocation_point(0);
            const double newmax = b.collocation_point(b.n() - 1);
            const double oldmin = Fb->collocation_point(0);
            const double oldmax = Fb->collocation_point(Fb->n() - 1);
            // Only pay for floating point loss if strictly necessary
#pragma warning(push,disable:1572)
            if (oldmin != newmin || oldmax != newmax) {
#pragma warning(pop)
                points = (points - newmin)
                    * ((oldmax - oldmin) / (newmax - newmin)) + oldmin;
            }

            // Step 2: Evaluate old basis + coefficients at points into new
            const int jmax = numeric_cast<int>(field.shape()[1]);
            const int kmax = numeric_cast<int>(field.shape()[2]);
            for (int k = 0; k < kmax; ++k) {
                for (int j = 0; j < jmax; ++j) {
                    Fb->linear_combination(0, &(*tmp)[0][j][k], b.n(),
                                           points.data(), &field[0][j][k], 0);
                }
            }

            // Step 3: Invert for new coefficients using factored mass matrix
            for (int k = 0; k < kmax; ++k) {
                mass->solve(jmax, &field[0][0][k],
                            numeric_cast<int>(field.strides()[0]),
                            numeric_cast<int>(field.strides()[1]));
            }
        }
    }

    DEBUG0("Finished loading simulation fields");
}

NoiseDefinition::NoiseDefinition(real_t fluctpercent,
                                 unsigned long fluctseed)
    : IDefinition("Additive random velocity perturbations on startup"),
      fluctpercent(fluctpercent),
      fluctseed(fluctseed)
{
    using ::suzerain::validation::ensure_positive;
    using ::suzerain::validation::ensure_nonnegative;
    ::std::pointer_to_binary_function<unsigned long,const char*,void>
        ptr_fun_ensure_positive_ulint(ensure_positive<unsigned long>);
    ::std::pointer_to_binary_function<real_t,const char*,void>
        ptr_fun_ensure_nonnegative_real(ensure_nonnegative<real_t>);
    this->add_options()
        ("fluctpercent",
         boost::program_options::value(&this->fluctpercent)
            ->default_value(this->fluctpercent)
            ->notifier(std::bind2nd(ptr_fun_ensure_nonnegative_real,
                                   "fluctpercent")),
         "Maximum fluctuation magnitude to add as a percentage of"
         " centerline mean streamwise velocity")
        ("fluctseed",
         boost::program_options::value(&this->fluctseed)
            ->default_value(this->fluctseed)
            ->notifier(std::bind2nd(ptr_fun_ensure_positive_ulint,
                                    "fluctseed")),
         "RngStream generator seed (L'Ecuyer et al. 2002)");
}

void
add_noise(suzerain::NoninterleavedState<4,complex_t> &state,
          const NoiseDefinition& noisedef,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline &b,
          const suzerain::bsplineop& bop)
{

#pragma warning(push,disable:1572)
    if (noisedef.fluctpercent == 0) {
#pragma warning(pop)
        DEBUG0("Zero noise added to velocity fields");
        return;
    }

    using suzerain::inorder::wavenumber_translatable;
    using suzerain::NoninterleavedState;
    namespace ndx = channel::field::ndx;
    const real_t twopi = 2 * boost::math::constants::pi<real_t>();
    Eigen::VectorXc scratch;

    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed collocation-based operator matrices
    assert(bop.get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    // Evaluate maximum fluctuation magnitude based on percentage of
    // centerline mean streamwise velocity and broadcast result.
    // Alright, alright... actually an approximate mean velocity.
    real_t maxfluct;
    if (dgrid.local_wave_start.x() == 0 && dgrid.local_wave_start.z() == 0) {
        assert(suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0);
        scratch.setConstant(2, 0);
        const real_t centerline = scenario.Ly / 2;
        b.linear_combination(
                0, &state[ndx::rhou][0][0][0], 1, &centerline, scratch.data());
        INFO0("Centerline mean streamwise momentum at y = "
              << centerline << " is " << scratch[0]);
        b.linear_combination(
                0, &state[ndx::rho][0][0][0], 1, &centerline, scratch.data()+1);
        INFO0("Centerline mean density at y = "
              << centerline << " is " << scratch[1]);
        maxfluct = noisedef.fluctpercent / 100
                 * (abs(scratch[0]) / abs(scratch[1]));
    }
    SUZERAIN_MPICHKR(MPI_Bcast(&maxfluct, 1,
                suzerain::mpi::datatype_of(maxfluct), 0, MPI_COMM_WORLD));
    INFO0("Adding velocity perturbations with maximum magnitude " << maxfluct);

    // Form mass matrix to convert (wave, collocation point values, wave)
    // perturbations to (wave, coefficients, wave)
    suzerain::bsplineop_luz massluz(bop);
    massluz.form_mass(bop);

    // Obtain integration coefficients for computing mean quantities across Y
    Eigen::VectorXr meancoeff(grid.N.y());
    b.integration_coefficients(0, meancoeff.data());
    meancoeff /= scenario.Ly;

    // Set L'Ecuyer et al.'s RngStream seed.  Use a distinct Substream for each
    // wall-normal pencil to ensure process is a) repeatable despite changes in
    // processor count, b) easy to code, and c) embarrassingly parallel.
    suzerain::RngStream rng;
    {
        boost::array<unsigned long,6> seed;
        std::fill(seed.begin(), seed.end(), noisedef.fluctseed);
        rng.SetSeed(seed.data());
    }
    rng.IncreasedPrecis(true);  // Use more bits of resolution

    // Approach is the following:
    //  0) Allocate storage for state and three additional scalar fields.
    //  1) Generate a random vector-valued field \tilde{A}.
    //     \tilde{A}'s x- and z-derivatives have zero mean by periodicity;
    //  2) Zero first two B-spline coefficients near walls so partial_y
    //     \tilde{A} vanishes at the wall
    //  3) Compute mean of \partial_y \tilde{A}(i, ., k) and subtract
    //     y * mean from \tilde{A}.  Again zero first two B-spline
    //     coefficients near walls.  Mean of \partial_y A is now approx zero.
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
    boost::scoped_ptr<NoninterleavedState<4,complex_t> > _s_ptr( // RAII
            allocate_padded_state<NoninterleavedState<4,complex_t> >(
                field::count + 3, dgrid));
    NoninterleavedState<4,complex_t> &s = *_s_ptr;               // Shorthand

    //  1) Generate a random vector-valued field \tilde{A}.
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

            // Compute local indices for global (i, ., k).
            const int local_i = i - dgrid.local_wave_start.x();
            const int local_k = k - dgrid.local_wave_start.z();

            // For each scalar component of \tilde{A}...
            for (std::size_t l = 0; l <  3; ++l) {

                // ...generate coeffs with well-defined pseudorandom order.
                //  2) Zero first two B-spline coefficients near walls
                //     so partial_y \tilde{A} vanishes at the wall
                scratch.setConstant(grid.N.y(), 0);
                for (int j = 2; j < grid.N.y() - 2; ++j) {
                    const real_t magnitude = rng.RandU01();
                    const real_t phase     = rng.RandU01() * twopi;
                    scratch[j] = std::polar(magnitude, phase);
                }

                //  3) Compute mean of \partial_y \tilde{A}(i, ., k) and
                //  subtract y * mean from \tilde{A}.  Again zero first two
                //  B-spline coefficients near walls.  Mean of \partial_y A is
                //  now approximately complex zero.
                bop.accumulate(1, 1,
                               1.0, scratch.data(), 1, grid.N.y(),
                               0.0, &s[2*l][0][local_i][local_k], 1, grid.N.y());
                massluz.solve(1, &s[2*l][0][local_i][local_k], 1, grid.N.y());
                std::complex<real_t> mean = 0;
                for (int m = 0; m < grid.N.y(); ++m) {
                    mean += meancoeff[m] * state[2*l][m][local_i][local_k];
                }
                bop.accumulate(1, 0,
                               1.0, scratch.data(), 1, grid.N.y(),
                               0.0, &s[2*l][0][local_i][local_k], 1, grid.N.y());
                for (int m = 0; m < grid.N.y(); ++m) {
                    s[2*l][m][local_i][local_k] -= b.collocation_point(m) * mean;
                }
                massluz.solve(1, &s[2*l][0][local_i][local_k], 1, grid.N.y());
                s[2*l][0           ][local_i][local_k] = 0;
                s[2*l][1           ][local_i][local_k] = 0;
                s[2*l][grid.N.y()-2][local_i][local_k] = 0;
                s[2*l][grid.N.y()-1][local_i][local_k] = 0;

                // Copy so coefficients are stored in s[2l] and s[2l+1]
                // as we need two copies to compute components of curl A.
                for (int m = 0; m < grid.N.y(); ++m) {
                    s[2*l+1][m][local_i][local_k] = s[2*l][m][local_i][local_k];
                }

            } // end scalar components of A

        } // end X

    } // end Z

    //  4) Compute curl A in physical space and rescale so maximum
    //     pointwise norm of A is maxfluct.  curl A is solenoidal
    //     and now has the velocity perturbation properties we desire.

    // Prepare physical-space view of the wave-space storage
    physical_view<field::count+3>::type p
        = physical_view<field::count+3>::create(dgrid, s);

    // Initializing OperatorBase to access decomposition-ready utilities
    suzerain::OperatorBase<real_t> obase(scenario, grid, dgrid, b, bop);

    // From Ax in s[0] compute \partial_y Ax
    obase.bop_apply(1, 1.0, s, 0);
    dgrid.transform_wave_to_physical(&p(0,0));

    // From Ax in s[1]  compute \partial_z Ax
    obase.bop_apply(0, 1.0, s, 1);
    obase.diffwave_apply(0, 1, 1.0, s, 1);
    dgrid.transform_wave_to_physical(&p(1,0));

    // From Ay in s[2] compute \partial_x Ay
    obase.bop_apply(0, 1.0, s, 2);
    obase.diffwave_apply(1, 0, 1.0, s, 2);
    dgrid.transform_wave_to_physical(&p(2,0));

    // From Ay in s[3] compute \partial_z Ay
    obase.bop_apply(0, 1.0, s, 3);
    obase.diffwave_apply(0, 1, 1.0, s, 3);
    dgrid.transform_wave_to_physical(&p(3,0));

    // From Az in s[4] compute \partial_x Az
    obase.bop_apply(0, 1.0, s, 4);
    obase.diffwave_apply(1, 0, 1.0, s, 4);
    dgrid.transform_wave_to_physical(&p(4,0));

    // From Az in s[5] compute \partial_y Az
    obase.bop_apply(1, 1.0, s, 5);
    dgrid.transform_wave_to_physical(&p(5,0));

    // Store Curl A in s[{5,6,7}] and find global maximum magnitude of curl A
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
                if (j == 0 || j == grid.N.y() - 1) {
#pragma warning(push,disable:1572)
                    assert(p(0, offset) == 0.0);
                    assert(p(1, offset) == 0.0);
                    assert(p(2, offset) == 0.0);
                    assert(p(3, offset) == 0.0);
                    assert(p(4, offset) == 0.0);
                    assert(p(5, offset) == 0.0);
#pragma warning(pop,disable:1572)
                }

                const real_t curlA_x        = p(5, offset) - p(3, offset);
                const real_t curlA_y        = p(1, offset) - p(4, offset);
                const real_t curlA_z        = p(2, offset) - p(0, offset);
                const real_t magsquared     = curlA_x * curlA_x
                                            + curlA_y * curlA_y
                                            + curlA_z * curlA_z;

                //  5) Store curl A in physical space in the 3 scalar fields.
                p(field::count + 0, offset) = curlA_x;
                p(field::count + 1, offset) = curlA_y;
                p(field::count + 2, offset) = curlA_z;

                maxmagsquared = suzerain::math::maxnan(
                        maxmagsquared,magsquared);

            } // end X

        } // end Z

    } // end Y
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, &maxmagsquared, 1,
                suzerain::mpi::datatype<real_t>::value,
                MPI_MAX, MPI_COMM_WORLD));

    // Rescale curl A components so max ||curl A|| == maxfluct
    p.row(field::count + 0) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(field::count + 1) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(field::count + 2) *= (maxfluct / std::sqrt(maxmagsquared));

    //  6) Copy state into auxiliary state storage and bring to
    //     physical space.
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        s[i] = state[i];
        obase.bop_apply(0, 1.0, s, i);
        dgrid.transform_wave_to_physical(&p(i,0));
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

                namespace rholut = suzerain::rholut;

                // Retrieve internal energy
                const real_t rho( p(ndx::rho,  offset));
                Eigen::Vector3r m(p(ndx::rhou, offset),
                                  p(ndx::rhov, offset),
                                  p(ndx::rhow, offset));
                real_t       e  ( p(ndx::rhoe, offset));
                const real_t e_int = rholut::energy_internal(Ma, rho, m, e);

                // Perturb momentum and compute updated total energy
                m.x() += rho * p(field::count + 0, offset);
                m.y() += rho * p(field::count + 1, offset);
                m.z() += rho * p(field::count + 2, offset);
                const real_t e_kin = rholut::energy_kinetic(Ma, rho, m);
                e = e_int + e_kin;

                // Store results back to state fields
                p(ndx::rhou, offset) = m.x();
                p(ndx::rhov, offset) = m.y();
                p(ndx::rhow, offset) = m.z();
                p(ndx::rhoe, offset) = e;

            } // end X

        } // end Z

    } // end Y

    //  8) Bring perturbed state information back to wavespace.
    // Build FFT normalization constant into Y direction's mass matrix.
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.form(1, &scale_factor, bop);
    assert(field::ndx::rho == 0);
    assert(static_cast<int>(field::ndx::rho) + 1 == field::ndx::rhou);
    for (std::size_t i = field::ndx::rhou; i < field::count; ++i) {
        dgrid.transform_physical_to_wave(&p(i, 0));      // X, Z
        obase.bop_solve(massluz, s, i);                  // Y
    }

    //  9) Overwrite state storage with the new perturbed state.
    assert(field::ndx::rho == 0);
    assert(static_cast<int>(field::ndx::rho) + 1 == field::ndx::rhou);
    for (std::size_t i = field::ndx::rhou; i < channel::field::count; ++i) {
        state[i] = s[i];
    }
}

boost::array<L2,field::count>
field_L2(const suzerain::NoninterleavedState<4,complex_t> &state,
         const suzerain::problem::ScenarioDefinition<real_t>& scenario,
         const suzerain::problem::GridDefinition& grid,
         const suzerain::pencil_grid& dgrid,
         const suzerain::bsplineop& gop)
{
    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed Galerkin L2 operator matrices
    assert(gop.get()->method == SUZERAIN_BSPLINEOP_GALERKIN_L2);

    // Only want non-dealiased X-direction modes to contribute to L2
    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(grid.N.x(),
                                            grid.dN.x(),
                                            dgrid.local_wave_start.x(),
                                            dgrid.local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Only want non-dealiased Z-direction modes to contribute to L2
    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(grid.N.z(),
                                            grid.dN.z(),
                                            dgrid.local_wave_start.z(),
                                            dgrid.local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // Hosed if MPI_COMM_WORLD rank 0 does not contain the zero-zero mode!
    if (mxb[0] == 0 && mzb[0] == 0) {
        assert(suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0);
    }

    // Temporary storage for inner product computations
    Eigen::VectorXc tmp;
    tmp.setZero(grid.N.y());

    // Temporary storage for accumulating and broadcasting results
    complex_t buf[2*field::count];
    typedef Eigen::Map<Eigen::Array<complex_t,field::count,1> > results_type;
    results_type total2(buf);
    results_type mean2(buf + field::count);

    // Compute the local L2 contribution towards each L^2 norm squared
    total2.setZero();
    for (size_t k = 0; k < field::count; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int n = mzb[j]; n < mze[j]; ++n) {
                for (int m = mxb[0]; m < mxe[0]; ++m) {
                    const complex_t * u_mn
                        = &state[k][0][m - dgrid.local_wave_start.x()]
                                      [n - dgrid.local_wave_start.z()];
                    gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
                    complex_t dot = suzerain::blas::dot(
                            grid.N.y(), u_mn, 1, tmp.data(), 1);
                    if (m > 0 && m < grid.dN.x()/2) {
                        dot *= 2;
                    }
                    total2[k] += dot;
                }
            }
        }
    }
    total2 *= scenario.Lx * scenario.Lz;

    // Reduce total2 sum onto processor housing the zero-zero mode using
    // mean2 as a scratch buffer to simulate (the disallowed) MPI_IN_PLACE
    SUZERAIN_MPICHKR(MPI_Reduce(total2.data(),
                mean2.data(), field::count * sizeof(complex_t)/sizeof(real_t),
                suzerain::mpi::datatype<real_t>(),
                MPI_SUM, 0, MPI_COMM_WORLD));
    total2 = mean2;

    // Compute the mean-only L^2 squared for each field on the root processor
    if (mzb[0] == 0 && mxb[0] == 0) {
        for (size_t k = 0; k < field::count; ++k) {
            const complex_t * u_mn = &state[k][0][0][0];
            gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
            mean2[k] = suzerain::blas::dot(grid.N.y(), u_mn, 1, tmp.data(), 1);
        }
        mean2 *= scenario.Lx * scenario.Lz;
    }

    // Broadcast total2 and mean2 values to all processors
    SUZERAIN_MPICHKR(MPI_Bcast(
                buf,
                sizeof(buf)/sizeof(buf[0]) * sizeof(complex_t)/sizeof(real_t),
                suzerain::mpi::datatype<real_t>(),
                0, MPI_COMM_WORLD));

    // Obtain fluctuating2 = total2 - mean2 and pack the return structure
    boost::array<L2,field::count> retval;
    for (size_t k = 0; k < field::count; ++k) {
        retval[k].mean2        = std::abs(mean2[k]);
        retval[k].fluctuating2 = std::abs(total2[k] - mean2[k]);
    }

    return retval;
}

void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        suzerain::NoninterleavedState<4,complex_t> &swave,
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        const real_t simulation_time)
{
    // Initializing OperatorBase to access decomposition-ready utilities
    suzerain::OperatorBase<real_t> obase(scenario, grid, dgrid, b, bop);

    // Prepare physical-space view of the wave-space storage
    physical_view<field::count>::type sphys
        = physical_view<field::count>::create(dgrid, swave);

    // Depending on whether or not we need swave's data...
#pragma warning(push,disable:1572)
    if (beta == 0) {
#pragma warning(pop)
        // ...either clear out any lingering NaNs or Infs from storage...
        suzerain::multi_array::fill(swave, 0);
    } else {
        // ...or scale data by beta and transform it to physical space.
        for (std::size_t i = 0; i < field::count; ++i) {
            obase.bop_apply(0, beta, swave, i);
            dgrid.transform_wave_to_physical(&sphys(i,0));
        }
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
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

                // Initialize manufactured solution's conserved state...
                const real_t rho  = msoln.rho (x, y, z, simulation_time);
                const real_t rhou = msoln.rhou(x, y, z, simulation_time);
                const real_t rhov = msoln.rhov(x, y, z, simulation_time);
                const real_t rhow = msoln.rhow(x, y, z, simulation_time);
                const real_t rhoe = msoln.rhoe(x, y, z, simulation_time);

                // ...and accumulate it into the desired storage
                sphys(field::ndx::rho,  offset) += alpha * rho;
                sphys(field::ndx::rhou, offset) += alpha * rhou;
                sphys(field::ndx::rhov, offset) += alpha * rhov;
                sphys(field::ndx::rhow, offset) += alpha * rhow;
                sphys(field::ndx::rhoe, offset) += alpha * rhoe;

            } // end X

        } // end Z

    } // end Y

    // Now convert physical state back to wave space coefficients,
    // building FFT normalization constant into Y direction's mass matrix.
    suzerain::bsplineop_luz massluz(bop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.form(1, &scale_factor, bop);
    for (std::size_t i = 0; i < field::count; ++i) {
        dgrid.transform_physical_to_wave(&sphys(i, 0));      // X, Z
        obase.bop_solve(massluz, swave, i);                  // Y
    }
}

} // end namespace channel
