//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// grid_definition.cpp: classes handling grid definitions
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/grid_definition.hpp>

#include <suzerain/diffwave.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

template<typename T>
static void parse_option(const std::string& s,
                         T* value, void (*validator)(T, const char*),
                         const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

static const char grid_definition_description[]
    = "Mixed Fourier/B-spline computational grid definition";

grid_definition::grid_definition()
    : grid_specification()
    , definition_base(grid_definition_description)
{
    this->initialize_options("NaN", "NaN", "NaN");  // Must match L!
}

grid_definition::grid_definition(const char* Lx,
                                 int         Nx,
                                 real_t      DAFx,
                                 const char* Ly,
                                 int         Ny,
                                 int         k,
                                 real_t      htdelta,
                                 const char* Lz,
                                 int         Nz,
                                 real_t      DAFz)
    : grid_specification(exprparse<real_t>(Lx, "grid_definition(..., Lx, ...)"),
                         Nx,
                         DAFx,
                         exprparse<real_t>(Ly, "grid_definition(..., Ly, ...)"),
                         Ny,
                         k,
                         htdelta,
                         exprparse<real_t>(Lz, "grid_definition(..., Lz, ...)"),
                         Nz,
                         DAFz)
    , definition_base("Mixed Fourier/B-spline computational grid definition")
{
    this->initialize_options(Lx, Ly, Lz);
}

void grid_definition::initialize_options(
        const char* default_Lx,
        const char* default_Ly,
        const char* default_Lz)
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::math::isnan;
    using std::auto_ptr;
    using std::bind1st;
    using std::bind2nd;
    using std::mem_fun;
    using std::ptr_fun;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    // Used to help resolve pointers-to-members taking strings
    grid_specification& (grid_specification::*f)(const std::string&) = NULL;

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;
    std::string* nullstr = NULL;

    // Lx
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<real_t>, _1, &L.x(),
                     &ensure_positive<real_t>, "Lx"));
    if (default_Lx)
    {
        p->default_value(default_Lx);
    }
    this->add_options()("Lx", p.release(),
                        "Nondimensional domain length in streamwise X direction");

    // Nx
    p.reset(boost::program_options::value(nullstr));
    f = &grid_specification::Nx;
    p->notifier(bind(f, this, _1));
    if (N.x())
    {
        p->default_value(lexical_cast<string>(N.x()));
    }
    this->add_options()("Nx", p.release(),
                        "Global logical extents in streamwise X direction");

    // DAFx
    p.reset(boost::program_options::value(nullstr));
    f = &grid_specification::DAFx;
    p->notifier(bind(f, this, _1));
    if (!(isnan)(DAF.x()))
    {
        p->default_value(lexical_cast<string>(DAF.x()));
    }
    this->add_options()("DAFx", p.release(),
                        "Dealiasing factor in streamwise X direction");

    // Ly
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<real_t>, _1, &L.y(),
                     &ensure_positive<real_t>, "Ly"));
    if (default_Ly)
    {
        p->default_value(default_Ly);
    }
    this->add_options()("Ly", p.release(),
                        "Nondimensional domain length in wall normal Y direction");

    // Ny
    p.reset(boost::program_options::value(nullstr));
    f = &grid_specification::Ny;
    p->notifier(bind(f, this, _1));
    if (N.y())
    {
        p->default_value(lexical_cast<string>(N.y()));
    }
    this->add_options()("Ny", p.release(),
                        "Global logical extents in wall-normal Y direction");

    // k
    p.reset(boost::program_options::value(nullstr));
    if (k)
    {
        p->notifier(bind(&parse_option<int>, _1, &k,
                         &ensure_positive<int>, "k"));
        p->default_value(lexical_cast<string>(k));
    }
    else
    {
        p->notifier(bind(&parse_option<int>, _1, &k,
                         &ensure_nonnegative<int>, "k"));
    }
    this->add_options()("k", p.release(),
                        "Wall-normal B-spline order (4 is piecewise cubic)");

    // htdelta
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<real_t>, _1, &htdelta,
                     &ensure_nonnegative<real_t>, "htdelta"));
    if (!(isnan)(htdelta))
    {
        p->default_value(lexical_cast<string>(htdelta));
    }
    this->add_options()("htdelta", p.release(),
                        "Wall-normal breakpoint hyperbolic tangent stretching");

    // Lz
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<real_t>, _1, &L.z(),
                     &ensure_positive<real_t>, "Lz"));
    if (default_Lz)
    {
        p->default_value(default_Lz);
    }
    this->add_options()("Lz", p.release(),
                        "Nondimensional domain length in spanwise Z direction");

    // Nz
    p.reset(boost::program_options::value(nullstr));
    f = &grid_specification::Nz;
    p->notifier(bind(f, this, _1));
    if (N.z())
    {
        p->default_value(lexical_cast<string>(N.z()));
    }
    this->add_options()("Nz", p.release(),
                        "Global logical extents in spanwise Z direction");

    // DAFz
    p.reset(boost::program_options::value(nullstr));
    f = &grid_specification::DAFz;
    p->notifier(bind(f, this, _1));
    if (!(isnan)(DAF.z()))
    {
        p->default_value(lexical_cast<string>(DAF.z()));
    }
    this->add_options()("DAFz", p.release(),
                        "Dealiasing factor in spanwise Z direction");

    // Pa
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<int>, _1, &P[0],
                     &ensure_nonnegative<int>, "Pa"));
    if (P[0])
    {
        p->default_value(lexical_cast<string>(P[0]));
    }
    this->add_options()("Pa", p.release(),
                        "Processor count in the Pa decomposition direction");

    // Pb
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<int>, _1, &P[1],
                     &ensure_nonnegative<int>, "Pb"));
    if (P[1])
    {
        p->default_value(lexical_cast<string>(P[1]));
    }
    this->add_options()("Pb", p.release(),
                        "Processor count in the Pb decomposition direction");
}

void save(const esio_handle h,
          const grid_definition& grid)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Saving grid_definition parameters");

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Lx", &grid.L.x(), 0,
            grid.options().find("Lx",false).description().c_str());

    esio_line_write(h, "Nx", &grid.N.x(), 0,
            grid.options().find("Nx",false).description().c_str());

    esio_line_write(h, "DAFx", &grid.DAF.x(), 0,
            grid.options().find("DAFx",false).description().c_str());

    esio_line_write(h, "Ly", &grid.L.y(), 0,
            grid.options().find("Ly",false).description().c_str());

    esio_line_write(h, "Ny", &grid.N.y(), 0,
            grid.options().find("Ny",false).description().c_str());

    esio_line_write(h, "k", &grid.k, 0,
            grid.options().find("k",false).description().c_str());

    esio_line_write(h, "htdelta", &grid.htdelta, 0,
            grid.options().find("htdelta",false).description().c_str());

    esio_line_write(h, "Lz", &grid.L.z(), 0,
            grid.options().find("Lz",false).description().c_str());

    esio_line_write(h, "Nz", &grid.N.z(), 0,
            grid.options().find("Nz",false).description().c_str());

    esio_line_write(h, "DAFz", &grid.DAF.z(), 0,
            grid.options().find("DAFz",false).description().c_str());

    DEBUG0("Storing wavenumber vectors for Fourier bases");
    ArrayXc cbuf(std::max(grid.N.x(), grid.N.z()));

    // Obtain wavenumbers via computing 1*(i*kx)/i
    cbuf.fill(complex_t(1,0));
    diffwave::apply(1, 0, complex_t(0,-1), cbuf.data(),
            grid.L.x(), grid.L.z(),
            1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(h, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in streamwise X direction"); // Re(cbuf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    cbuf.fill(complex_t(1,0));
    diffwave::apply(0, 1, complex_t(0,-1), cbuf.data(),
            grid.L.x(), grid.L.z(),
            1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(h, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(cbuf.data()),
            2, "Wavenumbers in spanwise Z direction"); // Re(cbuf)

    DEBUG0("Storing collocation point vectors for Fourier bases");
    ArrayXr rbuf;

    // Obtain collocation points in x using [-Lx/2, Lx/2]) and dN.x()
    if (grid.dN.x() > 1) {
        rbuf = ArrayXr::LinSpaced(Sequential, grid.dN.x(), 0, grid.dN.x() - 1);
        rbuf *= grid.L.x() / grid.dN.x();
        rbuf -= grid.L.x() / 2;
    } else {
        rbuf = ArrayXr::Constant(grid.dN.x(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_x", rbuf.data(), 0,
            "Collocation points for the dealiased, streamwise X direction");

    // Obtain collocation points in z using [-Lz/2, Lz/2]) and dN.z()
    if (grid.dN.z() > 1) {
        rbuf = ArrayXr::LinSpaced(Sequential, grid.dN.z(), 0, grid.dN.z() - 1);
        rbuf *= grid.L.z() / grid.dN.z();
        rbuf -= grid.L.z() / 2;
    } else {
        rbuf = ArrayXr::Constant(grid.dN.z(), 0);
    }
    esio_line_establish(h, rbuf.size(), 0, (procid == 0 ? rbuf.size() : 0));
    esio_line_write(h, "collocation_points_z", rbuf.data(), 0,
            "Collocation points for the dealiased, spanwise Z direction");
}

void load(const esio_handle h,
          grid_definition& grid)
{
    DEBUG0("Loading grid_definition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (!(boost::math::isnan)(grid.L.x())) {
        INFO0("Overriding grid using Lx = " << grid.L.x());
    } else {
        esio_line_read(h, "Lx", &grid.L.x(), 0);
    }

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

    if (!(boost::math::isnan)(grid.L.y())) {
        INFO0("Overriding grid using Ly = " << grid.L.y());
    } else {
        esio_line_read(h, "Ly", &grid.L.y(), 0);
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

    if (!(boost::math::isnan)(grid.L.z())) {
        INFO0("Overriding grid using Lz = " << grid.L.z());
    } else {
        esio_line_read(h, "Lz", &grid.L.z(), 0);
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

} // end namespace support

} // end namespace suzerain
