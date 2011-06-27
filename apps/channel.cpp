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

#include "logger.hpp"
#include "channel.hpp"
#include "nsctpl_rholut.hpp"

// Manufactured solution classes explicitly instantiated for debugging
template class nsctpl_rholut::manufactured_solution<real_t>;
template class nsctpl_rholut::primitive<real_t>;

// TODO Refactor -0.0 handling for scenario.{alpha,beta}, grid.htdelta

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
           const boost::shared_ptr<
                 nsctpl_rholut::manufactured_solution<real_t> >& msoln)
{
    static const char location[] = "nsctpl_rholut::manufactured_solution";

    // Always write a flag to indicate whether or not a MS is in use
    const int in_use = msoln ? 1 : 0;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &in_use, 0,
            "Is a nsctpl_rholut::manufactured_solution in use?");

    // Only proceed if an MS is in use
    if (!in_use) return;

    DEBUG0("Storing nsctpl_rholut::manufactured_solution parameters");

    // Check parameters stored with the scenario not the manufactured solution
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
          boost::shared_ptr<
                 nsctpl_rholut::manufactured_solution<real_t> >& msoln)
{
    static const char location[] = "nsctpl_rholut::manufactured_solution";


    // Determine if a manufactured solution should be loaded
    int in_use;
    esio_line_establish(h, 1, 0, 1); // All ranks load
    esio_line_read(h, location, &in_use, 0);

    // Only proceed if an MS is in use
    if (!in_use) {
        msoln.reset();
        return;
    }

    DEBUG0("Loading nsctpl_rholut::manufactured_solution parameters");

    // Allocate storage and defensively NaN every parameter not explicitly
    // loaded below.  Protects us against accidentally missing new solution
    // parameters.
    msoln.reset(new nsctpl_rholut::manufactured_solution<real_t>());
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

    char name[8];
    char comment[127];

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
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.global_wave_extent.z());

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
        bsplines_same = (b.knot(j) == Fb->knot(j));
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

    // Compute wavenumber translation logistics for Y direction
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

        // Clear storage prior to load to zero not-loaded coefficents
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
            for (int i = 0; i < b.n(); ++i) points[i] = b.collocation_point(i);

            // Step 1: Affine transformation of points into old basis domain
            // Allows stretching of old range onto new range.
            const double newmin = b.collocation_point(0);
            const double newmax = b.collocation_point(b.n() - 1);
            const double oldmin = Fb->collocation_point(0);
            const double oldmax = Fb->collocation_point(Fb->n() - 1);
            // Only pay for floating point loss if strictly necessary
            if (oldmin != newmin || oldmax != newmax) {
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
                    const complex_t * u_mn = &state[k][0][m][n];
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

    // Reduce total2 sum onto processor housing the zero-zero mode
    SUZERAIN_MPICHKR(MPI_Reduce(MPI_IN_PLACE,
                total2.data(), total2.size() * sizeof(complex_t)/sizeof(real_t),
                suzerain::mpi::datatype<real_t>(),
                MPI_SUM, 0, MPI_COMM_WORLD));

    // Compute the mean-only L^2 squared for each field on the root processor
    if (mzb[0] == 0 && mxb[0] == 0) {
        for (size_t k = 0; k < field::count; ++k) {
            const complex_t * u_mn = &state[k][0][0][0];
            gop.accumulate(0, 1.0, u_mn, 1, 0.0, tmp.data(), 1);
            mean2[k] = suzerain::blas::dot(grid.N.y(), u_mn, 1,
                                                       tmp.data(), 1);
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

} // end namespace channel
