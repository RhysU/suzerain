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
#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/RngStream.hpp>

#include "logging.hpp"
#include "channel.hpp"

// Manufactured solution classes explicitly instantiated for debugging
template class nsctpl_rholut::manufactured_solution<real_t>;

using boost::numeric_cast;

namespace channel {

// See "Configuration" at http://logging.apache.org/log4cxx/index.html
const char log4cxx_config[] =
    "log4j.rootLogger=INFO, CONSOLE, LOG\n"
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %-10c %m%n\n"
    "log4j.appender.LOG=org.apache.log4j.FileAppender\n"
    "log4j.appender.LOG.append=true\n"
    "log4j.appender.LOG.filename=log.dat\n"
    "log4j.appender.LOG.layout=${log4j.appender.CONSOLE.layout}\n"
    "log4j.appender.LOG.layout.ConversionPattern=${log4j.appender.CONSOLE.layout.ConversionPattern}\n"
;

const boost::array<const char *,field::count> field::name = {{
    "rho", "rhou", "rhov", "rhow", "rhoe"
}};

const boost::array<const char *,field::count> field::description = {{
    "density",
    "streamwise momentum", "wall-normal momentum", "spanwise momentum",
    "total energy",
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
    Eigen::ArrayXc cbuf(std::max(grid.N.x(), grid.N.z()));

    // Obtain wavenumbers via computing 1*(i*kx)/i
    cbuf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), cbuf.data(),
            Lx, Lz, 1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(h, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in streamwise X direction"); // Re(cbuf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    cbuf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), cbuf.data(),
            Lx, Lz, 1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(h, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in spanwise Z direction"); // Re(cbuf)

    DEBUG0("Storing collocation point vectors for Fourier bases");
    Eigen::ArrayXr rbuf;

    // Obtain collocation points in x using [-Lx/2, Lx/2]) and dN.x()
    if (grid.dN.x() > 1) {
        rbuf = Eigen::ArrayXr::LinSpaced(
                Eigen::Sequential, grid.dN.x(), 0, grid.dN.x() - 1);
        rbuf *= Lx / grid.dN.x();
        rbuf -= Lx/2;
    } else {
        rbuf = Eigen::ArrayXr::Constant(grid.dN.x(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_x", rbuf.data(), 0,
            "Collocation points for the dealiased, streamwise X direction");

    // Obtain collocation points in z using [-Lz/2, Lz/2]) and dN.z()
    if (grid.dN.z() > 1) {
        rbuf = Eigen::ArrayXr::LinSpaced(
                Eigen::Sequential, grid.dN.z(), 0, grid.dN.z() - 1);
        rbuf *= Lz / grid.dN.z();
        rbuf -= Lz/2;
    } else {
        rbuf = Eigen::ArrayXr::Constant(grid.dN.z(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_z", rbuf.data(), 0,
            "Collocation points for the dealiased, spanwise Z direction");
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
            k, suzerain::bspline::from_breakpoints(),
            breakpoints.size(), breakpoints.data());
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
    esio_line_write(h, "breakpoints_y", buf.data(), 0,
            "Breakpoint locations used to build wall-normal B-spline basis");

    for (int i = 0; i < b->n(); ++i) buf[i] = b->collocation_point(i);
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "collocation_points_y", buf.data(), 0,
            "Collocation points used to build wall-normal discrete operators");

    b->integration_coefficients(0, buf.data());
    esio_line_establish(h, b->n(), 0, (procid == 0 ? b->n() : 0));
    esio_line_write(h, "integration_weights", buf.data(), 0,
            "Integrate by dotting B-spline coefficients against weights");

    char name[8]      = {};
    char comment[127] = {};

    for (int k = 0; k <= bop->nderiv(); ++k) {
        snprintf(name, sizeof(name), "Dy%d", k);
        snprintf(comment, sizeof(comment),
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
        snprintf(name, sizeof(name), "Gy%d", k);
        snprintf(comment, sizeof(comment),
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

    // All ranks load B-spline breakpoints_y (with backward compatibility)
    Eigen::ArrayXr breakpoints;
    const char *names[] = { "breakpoints_y", "breakpoints" };
    for (std::size_t i = 0; i < SUZERAIN_COUNTOF(names); ++i) {
        int nbreak;
        if (ESIO_NOTFOUND == esio_line_size(h, names[i], &nbreak)) {
            DEBUG0("Wall-normal breakpoints not found at /" << names[i]);
            continue;
        }
        esio_line_establish(h, nbreak, 0, nbreak);
        breakpoints.resize(nbreak);
        esio_line_read(h, names[i], breakpoints.data(), 0);
        break;
    }
    if (!breakpoints.size()) {
        SUZERAIN_ERROR_VOID(
                "Restart did not contain wall-normal breakpoint locations",
                SUZERAIN_EFAILED);
    }

    // Collocation points are ignored

    // Construct B-spline workspace
    b = boost::make_shared<suzerain::bspline>(
                k, suzerain::bspline::from_breakpoints(),
                breakpoints.size(), breakpoints.data());
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
suzerain::ContiguousState<4,complex_t>* allocate_padded_state(
           const std::size_t howmany_fields,
           const suzerain::pencil_grid& dgrid)
{
    // Create instance with appropriate padding to allow P3DFFTification
    suzerain::ContiguousState<4,complex_t> * const retval =
        new suzerain::ContiguousState<4,complex_t>(
            suzerain::to_yxz(howmany_fields, dgrid.local_wave_extent),
            suzerain::prepend(dgrid.local_wave_storage(), suzerain::strides_cm(
                              suzerain::to_yxz(dgrid.local_wave_extent)))
        );

    // Clear to avoid lingering NaN issues
    suzerain::multi_array::fill(*retval, 0);

    return retval;
}

void store_coefficients(
        const esio_handle h,
        const suzerain::ContiguousState<4,complex_t> &state,
        suzerain::ContiguousState<4,complex_t> &scratch,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid)
{
    // Ensure state and scratch storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());
    assert(scratch.shape()[0]  >= 1);
    assert(scratch.strides()[1] == state.strides()[1]);
    assert(scratch.strides()[2] == state.strides()[2]);
    assert(scratch.strides()[3] == state.strides()[3]);

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

        // ...first generating a metadata comment...
        std::string comment = "Nondimensional ";
        comment += field::description[i];
        comment += " stored row-major ZXY using a Fourier basis in Z stored"
                   " in-order per /kz; a Fourier basis in X stored using"
                   " Hermitian symmetry per /kx; and a B-spline basis in Y"
                   " defined by /k, /breakpoints_y, and /knots";

        // ...next eliminating dealiased content via a copy into scratch[0]...
        suzerain::diffwave::accumulate(
                0, 0,
                complex_t(1), state[i].origin(),
                complex_t(0), scratch[0].origin(),
                scenario.Lx, scenario.Lz, dgrid.global_wave_extent.y(),
                grid.N.x(), grid.dN.x(),
                dgrid.local_wave_start.x(), dgrid.local_wave_end.x(),
                grid.N.z(), grid.dN.z(),
                dgrid.local_wave_start.z(), dgrid.local_wave_end.z());

        // ...followed by two collective writes per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &scratch[0][0][mxb[0] - dgrid.local_wave_start.x()]
                                    [mzb[j] - dgrid.local_wave_start.z()];
            }

            // Collectively establish size of write across all ranks
            esio_field_establish(h,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Perform collective write operation from state_linear
            complex_field_write(h, field::name[i], src,
                                scratch.strides()[3],
                                scratch.strides()[2],
                                scratch.strides()[1],
                                comment.c_str());
        }
    }
}

// Check if two B-spline bases are identical.  Use strict equality as minor
// knot differences magnify once collocation points and operators are computed.
static bool bspline_bases_identical(const suzerain::bspline& a,
                                    const suzerain::bspline& b)
{
    bool retval =    a.k()     == b.k()
                  && a.n()     == b.n()
                  && a.nknot() == b.nknot();
    for (int j = 0; retval && j < b.nknot(); ++j) {
#pragma warning(push,disable:1572)
        retval = (a.knot(j) == b.knot(j));
#pragma warning(pop)
    }

    return retval;
}

void load_coefficients(const esio_handle h,
                       suzerain::ContiguousState<4,complex_t> &state,
                       const suzerain::problem::GridDefinition& grid,
                       const suzerain::pencil_grid& dgrid,
                       const suzerain::bspline& b,
                       const suzerain::bsplineop& bop)
{
    typedef suzerain::ContiguousState<4,complex_t> load_type;

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

    // Check if the B-spline basis in the file differs from ours.
    const bool bsplines_same = bspline_bases_identical(b, *Fb);

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
}

void store_collocation_values(
        const esio_handle h,
        const suzerain::ContiguousState<4,complex_t>& state,
        suzerain::ContiguousState<4,complex_t>& scratch,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid,
        suzerain::bspline& b,
        const suzerain::bsplineop& bop)
{
    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());
    assert(                  scratch.shape()[0]  >= field::count);
    assert(numeric_cast<int>(scratch.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(scratch.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(scratch.shape()[3]) == dgrid.local_wave_extent.z());

    // Copy-and-convert coefficients into collocation point values
    // Transforms from full-wave in state to full-physical in scratch
    suzerain::OperatorBase<real_t> obase(scenario, grid, dgrid, b, bop);
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        obase.diffwave_accumulate(0, 0, 1., state, i, 0., scratch, i);
        obase.bop_apply(0, 1, scratch, i);
        dgrid.transform_wave_to_physical(
                reinterpret_cast<real_t *>(scratch[i].origin()));
    }

    // Convert conserved rho, rhou, rhov, rhow, rhoe into u, v, w, p, T
    physical_view<field::count>::type sphys
        = physical_view<field::count>::create(dgrid, scratch);

    const real_t alpha = scenario.alpha;
    const real_t beta  = scenario.beta;
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack conserved quantities from fields
        real_t rho =      sphys(field::ndx::rho,  o);
        Eigen::Vector3r m(sphys(field::ndx::rhou, o),
                          sphys(field::ndx::rhov, o),
                          sphys(field::ndx::rhow, o));
        real_t e =        sphys(field::ndx::rhoe, o);

        // Compute primitive quantities to be stored
        real_t p, T;
        suzerain::rholut::p_T(alpha, beta, gamma, Ma, rho, m, e, p, T);
        m /= rho;

        // Pack primitive quantities back into fields (by position)
        sphys(0, o) = m.x(); // Now just X velocity
        sphys(1, o) = m.y(); // Now just Y velocity
        sphys(2, o) = m.z(); // Now just Z velocity
        sphys(3, o) = p;
        sphys(4, o) = T;
    }

    // HDF5 file storage locations and corresponding descriptions
    const boost::array<const char *,5> prim_names = {{
        "u", "v", "w", "p", "T"
    }};
    const boost::array<const char *,5> prim_descriptions = {{
        "X velocity", "Y velocity", "Z velocity", "pressure", "temperature",
    }};
    assert(prim_names.static_size == prim_descriptions.static_size);

    // Establish size of collective writes across all ranks and write fields
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());

    for (size_t i = 0; i < prim_names.static_size; ++i) {

        std::string comment = "Nondimensional ";
        comment += prim_descriptions[i];
        comment += " stored row-major YZX on the 3D rectilinear grid defined"
                   " by taking the outer product of arrays"
                   " /collocation_points_y, /collocation_points_z, and"
                   " /collocation_points_z";

        esio_field_write(h, prim_names[i],
                reinterpret_cast<real_t *>(scratch[i].origin()),
                0, 0, 0, comment.c_str());
    }

    // TODO Save mean primitive state at wall-normal collocation points
}

void load_collocation_values(
        const esio_handle h,
        suzerain::ContiguousState<4,complex_t>& state,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid,
        suzerain::bspline& b,
        const suzerain::bsplineop& bop)
{
    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field::count);
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

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

        // Check that restart file specifies the same B-spline basis
        // TODO Too restrictive-- identical collocation points would be okay
        boost::shared_ptr<suzerain::bspline> Fb;
        boost::shared_ptr<suzerain::bsplineop> Fbop;
        load(h, Fb, Fbop);
        if (!bspline_bases_identical(b, *Fb)) {
            ERROR0("Physical restart fields have different wall-normal basis");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }
    }

    // Establish size of collective reads across all ranks and read data
    physical_view<field::count>::type sphys
        = physical_view<field::count>::create(dgrid, state);
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

    // Convert primitive u, v, w, p, and T into rho, rhou, rhov, rhow, rhoe
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack primitive quantities from fields (by position)
        Eigen::Vector3r m(sphys(0, o),  // Now just X velocity
                          sphys(1, o),  // Now just Y velocity
                          sphys(2, o)); // Now just Z velocity
        const real_t p =  sphys(3, o);
        const real_t T =  sphys(4, o);

        // Compute conserved quantities from primitive ones
        const real_t rho = gamma * p / T;   // Assumes EOS
        m               *= rho;             // Not m contains momentum
        const real_t e   = suzerain::rholut::energy_kinetic(Ma, rho, m)
                            + suzerain::rholut::energy_internal(gamma, p);

        // Pack conserved quantities into fields (by name)
        sphys(field::ndx::rho,  o) = rho;
        sphys(field::ndx::rhou, o) = m.x();
        sphys(field::ndx::rhov, o) = m.y();
        sphys(field::ndx::rhow, o) = m.z();
        sphys(field::ndx::rhoe, o) = e;
    }

    // Initialize OperatorBase to access decomposition-ready utilities
    suzerain::OperatorBase<real_t> obase(scenario, grid, dgrid, b, bop);

    // Collectively convert physical state to wave space coefficients
    // Build FFT normalization constant into Y direction's mass matrix
    suzerain::bsplineop_luz massluz(bop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.form(1, &scale_factor, bop);

    for (std::size_t i = 0; i < field::count; ++i) {
        dgrid.transform_physical_to_wave(&sphys(i, 0)); // X, Z
        obase.bop_solve(massluz, state, i);             // Y
    }
}

void load(const esio_handle h,
          suzerain::ContiguousState<4,complex_t>& state,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline& b,
          const suzerain::bsplineop& bop)
{
    // Check whether load_coefficients(...) should work
    bool trycoeffs = true;
    for (std::size_t i = 0; i < field::count; ++i) {
        int ncomponents = 0;
        switch (esio_field_sizev(h, field::name[i], 0, 0, 0, &ncomponents)) {
            case ESIO_SUCCESS:
                if (ncomponents != 2) {
                    DEBUG0("Field /" << field::name[i] << " looks fishy...");
                }
                break;
            case ESIO_NOTFOUND:
                DEBUG0("Field /" << field::name[i] << " not found in restart");
                trycoeffs = false;
                break;
            default:
                DEBUG0("Field /" << field::name[i] << " looks fishy...");
                break;
        }
    }

    // Dispatch to the appropriate restart loading logic
    DEBUG0("Started loading simulation fields");
    if (trycoeffs) {
        load_coefficients(h, state, grid, dgrid, b, bop);
    } else {
        INFO0("Loading collocation-based, physical-space restart data");
        load_collocation_values(h, state, scenario, grid, dgrid, b, bop);
    }
    DEBUG0("Finished loading simulation fields");
}

/**
 * Parses "min:max", "min:[defaultmax]", or "[defaultmin]:max" into valmin, \c
 * valmax where \c absmin <= \c valmin <= \c valmax <= \c absmax is enforced
 * with the outer two inequalities being considered a validation failure.
 */
template<typename T>
static void parse_range(const std::string& s,
                        T *valmin, T *valmax,
                        const T defaultmin, const T defaultmax,
                        const T absmin, const T absmax,
                        const char *name)
{
    assert(absmin <= defaultmin);
    assert(defaultmin <= defaultmax);
    assert(defaultmax <= absmax);

    // Split s on a mandatory colon into whitespace-trimmed s_{min,max}
    const std::size_t colonpos = s.find_first_of(':');
    if (colonpos == std::string::npos) {
        throw std::invalid_argument(std::string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    }
    std::string s_min(s, 0, colonpos);
    std::string s_max(s, colonpos + 1);
    boost::algorithm::trim(s_min);
    boost::algorithm::trim(s_max);

    // Parse recognized formats into valmin and valmax
    if (s_min.length() == 0 && s_max.length() == 0) {
        throw std::invalid_argument(std::string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    } else if (s_min.length() == 0) {
        *valmin = defaultmin;
        *valmax = suzerain::exprparse<T>(s_max, name);
    } else if (s_max.length() == 0) {
        *valmin = suzerain::exprparse<T>(s_min, name);
        *valmax = defaultmax;
    } else {
        *valmin = suzerain::exprparse<T>(s_min, name);
        *valmax = suzerain::exprparse<T>(s_max, name);
    }

    // Ensure valmin <= valmax
    if (*valmin > *valmax) std::swap(*valmin, *valmax);

    // Validate range is within [absmin, absmax]
    if (*valmin < absmin || absmax < *valmax) {
        std::ostringstream oss;
        oss << name << " value [" << *valmin << ":" << *valmax
            << "] is outside valid range [" << absmin << ":" << absmax <<  "]";
        throw std::invalid_argument(oss.str());
    }
}

NoiseDefinition::NoiseDefinition(real_t percent,
                                 unsigned long seed)
    : IDefinition("Additive random velocity perturbations on startup"),
      percent(percent),
      kxfrac_min(0),
      kxfrac_max(1),
      kzfrac_min(0),
      kzfrac_max(1),
      seed(seed)
{
    using ::boost::bind;
    using ::suzerain::validation::ensure_positive;
    using ::suzerain::validation::ensure_nonnegative;
    ::std::pointer_to_binary_function<unsigned long,const char*,void>
        ptr_fun_ensure_positive_ulint(ensure_positive<unsigned long>);
    ::std::pointer_to_binary_function<real_t,const char*,void>
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
            ->notifier(bind(&parse_range<real_t>, _1,
                            &this->kxfrac_min, &this->kxfrac_max,
                            0, 1, 0, 1, "fluct_kxfrac")),
         "Range of X wavenumbers in which to generate fluctuations")
        ("fluct_kzfrac",
         boost::program_options::value<std::string>(0)
            ->default_value("0:1")
            ->notifier(bind(&parse_range<real_t>, _1,
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
add_noise(suzerain::ContiguousState<4,complex_t> &state,
          const NoiseDefinition& noisedef,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline &b,
          const suzerain::bsplineop& bop)
{
    const int Ny = grid.N.y();

#pragma warning(push,disable:1572)
    if (noisedef.percent == 0) {
#pragma warning(pop)
        DEBUG0("Zero noise added to velocity fields");
        return;
    }

    using suzerain::inorder::wavenumber_abs;
    using suzerain::inorder::wavenumber_max;
    using suzerain::inorder::wavenumber_translatable;
    using suzerain::ContiguousState;
    namespace ndx = channel::field::ndx;
    const real_t twopi = 2 * boost::math::constants::pi<real_t>();

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
        complex_t momentum, density;
        const real_t centerline = scenario.Ly / 2;
        b.linear_combination(
                0, &state[ndx::rhou][0][0][0], 1, &centerline, &momentum);
        INFO0("Centerline mean streamwise momentum at y = "
              << centerline << " is " << momentum);
        b.linear_combination(
                0, &state[ndx::rho][0][0][0], 1, &centerline, &density);
        INFO0("Centerline mean density at y = "
              << centerline << " is " << density);
        maxfluct = noisedef.percent / 100 * (abs(momentum) / abs(density));
    }
    SUZERAIN_MPICHKR(MPI_Bcast(&maxfluct, 1,
                suzerain::mpi::datatype_of(maxfluct), 0, MPI_COMM_WORLD));
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
    suzerain::bsplineop_luz massluz(bop);
    massluz.form_mass(bop);

    // Set L'Ecuyer et al.'s RngStream seed.  Use a distinct Substream for each
    // wall-normal pencil to ensure process is a) repeatable despite changes in
    // processor count, b) easy to code, and c) embarrassingly parallel.
    suzerain::RngStream rng;
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
    boost::scoped_ptr<ContiguousState<4,complex_t> > _s_ptr( // RAII
            allocate_padded_state<ContiguousState<4,complex_t> >(
                field::count + 3, dgrid));
    ContiguousState<4,complex_t> &s = *_s_ptr;               // Shorthand

    // 1) Generate a random vector-valued field \tilde{A}.
    // For each scalar component of \tilde{A}...
    for (std::size_t l = 0; l < 3; ++l) {

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
    for (std::size_t l = 0; l < 3; ++l) s[2*l+1] = s[2*l];

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
                Eigen::Vector3r curlA(p(5, offset) - p(3, offset),
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
                p(field::count + 0, offset) = curlA.z();
                p(field::count + 1, offset) = curlA.x();
                p(field::count + 2, offset) = curlA.y();

                maxmagsquared = suzerain::math::maxnan(
                        maxmagsquared, curlA.squaredNorm());

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
field_L2(const suzerain::ContiguousState<4,complex_t> &state,
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
                buf, SUZERAIN_COUNTOF(buf) * sizeof(complex_t)/sizeof(real_t),
                suzerain::mpi::datatype<real_t>(), 0, MPI_COMM_WORLD));

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
        suzerain::ContiguousState<4,complex_t> &swave,
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
