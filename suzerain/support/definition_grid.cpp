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
 * @copydoc definition_grid.hpp
 */

#include <suzerain/support/definition_grid.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/diffwave.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

template<typename T>
static void parse_option(
        const std::string& s,
        T* value,
        const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    *value = t;
}

template<typename T>
static void validate_option(
        const std::string& s,
        T* value,
        void (*validator)(T, const char*),
        const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

definition_grid::definition_grid()
    : specification_grid(std::numeric_limits<real_t>::quiet_NaN(),
                         0,
                         std::numeric_limits<real_t>::quiet_NaN(),
                         std::numeric_limits<real_t>::quiet_NaN(),
                         0,
                         0,
                         std::numeric_limits<real_t>::quiet_NaN(),
                         std::numeric_limits<real_t>::quiet_NaN(),
                         0,
                         std::numeric_limits<real_t>::quiet_NaN())
{
}

definition_grid::definition_grid(const real_t Lx,
                                 const int    Nx,
                                 const real_t DAFx,
                                 const real_t Ly,
                                 const int    Ny,
                                 const int    k,
                                 const real_t htdelta,
                                 const real_t Lz,
                                 const int    Nz,
                                 const real_t DAFz)
    : specification_grid(Lx,
                         Nx,
                         DAFx,
                         Ly,
                         Ny,
                         k,
                         htdelta,
                         Lz,
                         Nz,
                         DAFz)
{
}

// Strings used in options_description and populate/override/save/load.
static const char name_Lx[]      = "Lx";
static const char name_Nx[]      = "Nx";
static const char name_DAFx[]    = "DAFx";
static const char name_Ly[]      = "Ly";
static const char name_Ny[]      = "Ny";
static const char name_k[]       = "k";
static const char name_htdelta[] = "htdelta";
static const char name_Lz[]      = "Lz";
static const char name_Nz[]      = "Nz";
static const char name_DAFz[]    = "DAFz";
static const char name_Pa[]      = "Pa";
static const char name_Pb[]      = "Pb";

// Description used in options_description and populate/override/save/load.
static const char desc_Lx[]
        = "Domain length in streamwise X direction";

static const char desc_Nx[]
        = "Global logical extents in streamwise X direction";

static const char desc_DAFx[]
        = "Dealiasing factor in streamwise X direction";

static const char desc_Ly[]
        = "Domain length in wall-normal Y direction";

static const char desc_Ny[]
        = "Global logical extents in wall-normal Y direction";

static const char desc_k[]
        = "Wall-normal B-spline order (4 is piecewise cubic)";

static const char desc_htdelta[]
        = "Wall-normal breakpoint hyperbolic tangent stretching. "
          "Positive and negative values indicate two- and one-sided "
          "stretching, respectively.";

static const char desc_Lz[]
        = "Domain length in spanwise Z direction";

static const char desc_Nz[]
        = "Global logical extents in spanwise Z direction";

static const char desc_DAFz[]
        = "Dealiasing factor in spanwise Z direction";

static const char desc_Pa[]
        = "Processor count in the Pa decomposition direction";

static const char desc_Pb[]
        = "Processor count in the Pb decomposition direction";

boost::program_options::options_description
definition_grid::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    options_description retval(
            "Mixed Fourier/B-spline computational grid definition");

    // Used to help resolve pointers-to-members taking strings
    specification_grid& (specification_grid::*f)(const string&) = NULL;

    auto_ptr<typed_value<string> > p;

    // Lx
    p.reset(value<string>());
    p->notifier(bind(&validate_option<real_t>, _1, &L.x(),
                     &ensure_positive<real_t>, name_Lx));
    if (!(boost::math::isnan)(L.x())) {
        p->default_value(lexical_cast<string>(L.x()));
    }
    retval.add_options()(name_Lx, p.release(), desc_Lx);

    // Nx
    p.reset(value<string>());
    f = &specification_grid::Nx;
    p->notifier(bind(f, this, _1));
    if (N.x()) {
        p->default_value(lexical_cast<string>(N.x()));
    }
    retval.add_options()(name_Nx, p.release(), desc_Nx);

    // DAFx
    p.reset(value<string>());
    f = &specification_grid::DAFx;
    p->notifier(bind(f, this, _1));
    if (!(boost::math::isnan)(DAF.x())) {
        p->default_value(lexical_cast<string>(DAF.x()));
    }
    retval.add_options()(name_DAFx, p.release(), desc_DAFx);

    // Ly
    p.reset(value<string>(NULL));
    p->notifier(bind(&validate_option<real_t>, _1, &L.y(),
                     &ensure_positive<real_t>, name_Ly));
    if (!(boost::math::isnan)(L.y())) {
        p->default_value(lexical_cast<string>(L.y()));
    }
    retval.add_options()(name_Ly, p.release(), desc_Ly);

    // Ny
    p.reset(value<string>());
    f = &specification_grid::Ny;
    p->notifier(bind(f, this, _1));
    if (N.y()) {
        p->default_value(lexical_cast<string>(N.y()));
    }
    retval.add_options()(name_Ny, p.release(), desc_Ny);

    // k
    p.reset(value<string>());
    if (k) {
        p->notifier(bind(&validate_option<int>, _1, &k,
                         &ensure_positive<int>, name_k));
        p->default_value(lexical_cast<string>(k));
    } else {
        p->notifier(bind(&validate_option<int>, _1, &k,
                         &ensure_nonnegative<int>, name_k));
    }
    retval.add_options()(name_k, p.release(), desc_k);

    // htdelta
    p.reset(value<string>());
    p->notifier(bind(&parse_option<real_t>, _1, &htdelta, name_htdelta));
    if (!(boost::math::isnan)(htdelta)) {
        p->default_value(lexical_cast<string>(htdelta));
    }
    retval.add_options()(name_htdelta, p.release(), desc_htdelta);

    // Lz
    p.reset(value<string>());
    p->notifier(bind(&validate_option<real_t>, _1, &L.z(),
                     &ensure_positive<real_t>, name_Lz));
    if (!(boost::math::isnan)(L.z())) {
        p->default_value(lexical_cast<string>(L.z()));
    }
    retval.add_options()(name_Lz, p.release(), desc_Lz);

    // Nz
    p.reset(value<string>());
    f = &specification_grid::Nz;
    p->notifier(bind(f, this, _1));
    if (N.z()) {
        p->default_value(lexical_cast<string>(N.z()));
    }
    retval.add_options()(name_Nz, p.release(), desc_Nz);

    // DAFz
    p.reset(value<string>());
    f = &specification_grid::DAFz;
    p->notifier(bind(f, this, _1));
    if (!(boost::math::isnan)(DAF.z())) {
        p->default_value(lexical_cast<string>(DAF.z()));
    }
    retval.add_options()(name_DAFz, p.release(), desc_DAFz);

    // Pa
    p.reset(value<string>());
    p->notifier(bind(&validate_option<int>, _1, &P[0],
                     &ensure_nonnegative<int>, name_Pa));
    if (P[0]) {
        p->default_value(lexical_cast<string>(P[0]));
    }
    retval.add_options()(name_Pa, p.release(), desc_Pa);

    // Pb
    p.reset(value<string>());
    p->notifier(bind(&validate_option<int>, _1, &P[1],
                     &ensure_nonnegative<int>, name_Pb));
    if (P[1]) {
        p->default_value(lexical_cast<string>(P[1]));
    }
    retval.add_options()(name_Pb, p.release(), desc_Pb);

    return retval;
}

void
definition_grid::populate(
        const specification_grid& that,
        const bool verbose)
{
    maybe_populate(name_Lx, desc_Lx, L.x(), that.L.x(), verbose);
    {
        int value = N.x();
        maybe_populate(name_Nx, desc_Nx, value, that.N.x(), verbose);
        Nx(value);
    }
    {
        double factor = DAF.x();
        maybe_populate(name_DAFx, desc_DAFx, factor, that.DAF.x(), verbose);
        DAFx(factor);
    }
    maybe_populate(name_Ly, desc_Ly, L.y(), that.L.y(), verbose);
    {
        int value = N.y();
        maybe_populate(name_Ny, desc_Ny, value, that.N.y(), verbose);
        Ny(value);
    }
    maybe_populate(name_k, desc_k, k, that.k, verbose);
    maybe_populate(name_htdelta, desc_htdelta, htdelta, that.htdelta, verbose);
    maybe_populate(name_Lz, desc_Lz, L.z(), that.L.z(), verbose);
    {
        int value = N.z();
        maybe_populate(name_Nz, desc_Nz, value, that.N.z(), verbose);
        Nz(value);
    }
    {
        double factor = DAF.z();
        maybe_populate(name_DAFz, desc_DAFz, factor, that.DAF.z(), verbose);
        DAFz(factor);
    }
    maybe_populate(name_Pa, desc_Pa, P[0], that.P[0], verbose);
    maybe_populate(name_Pb, desc_Pb, P[1], that.P[1], verbose);
}

void
definition_grid::override(
        const specification_grid& that,
        const bool verbose)
{
    maybe_override(name_Lx, desc_Lx, L.x(), that.L.x(), verbose);
    {
        int value = N.x();
        maybe_override(name_Nx, desc_Nx, value, that.N.x(), verbose);
        Nx(value);
    }
    {
        double factor = DAF.x();
        maybe_override(name_DAFx, desc_DAFx, factor, that.DAF.x(), verbose);
        DAFx(factor);
    }
    maybe_override(name_Ly, desc_Ly, L.y(), that.L.y(), verbose);
    {
        int value = N.y();
        maybe_override(name_Ny, desc_Ny, value, that.N.y(), verbose);
        Ny(value);
    }
    maybe_override(name_k, desc_k, k, that.k, verbose);
    maybe_override(name_htdelta, desc_htdelta, htdelta, that.htdelta, verbose);
    maybe_override(name_Lz, desc_Lz, L.z(), that.L.z(), verbose);
    {
        int value = N.z();
        maybe_override(name_Nz, desc_Nz, value, that.N.z(), verbose);
        Nz(value);
    }
    {
        double factor = DAF.z();
        maybe_override(name_DAFz, desc_DAFz, factor, that.DAF.z(), verbose);
        DAFz(factor);
    }
    maybe_override(name_Pa, desc_Pa, P[0], that.P[0], verbose);
    maybe_override(name_Pb, desc_Pb, P[1], that.P[1], verbose);
}

void
definition_grid::save(
        const esio_handle h) const
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Saving definition_grid parameters");

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_Lx,      &this->L.x(),   0, desc_Lx     );
    esio_line_write(h, name_Nx,      &this->N.x(),   0, desc_Nx     );
    esio_line_write(h, name_DAFx,    &this->DAF.x(), 0, desc_DAFx   );
    esio_line_write(h, name_Ly,      &this->L.y(),   0, desc_Ly     );
    esio_line_write(h, name_Ny,      &this->N.y(),   0, desc_Ny     );
    esio_line_write(h, name_k,       &this->k,       0, desc_k      );
    esio_line_write(h, name_htdelta, &this->htdelta, 0, desc_htdelta);
    esio_line_write(h, name_Lz,      &this->L.z(),   0, desc_Lz     );
    esio_line_write(h, name_Nz,      &this->N.z(),   0, desc_Nz     );
    esio_line_write(h, name_DAFz,    &this->DAF.z(), 0, desc_DAFz   );

    DEBUG0("Storing wavenumber vectors for Fourier bases");
    ArrayXc cbuf(std::max(N.x()/2+1, N.z()));

    // Obtain Hermitian symmetric X wavenumbers via computing Re(1*(i*kx)/i)
    // Notice N.x() and dN.x() represent logical, not stored, extents.
    cbuf.fill(complex_t(1,0));
    diffwave::apply(1, 0, complex_t(0,-1), cbuf.data(),
            L.x(), L.z(),
            1, N.x(), dN.x(), 0, N.x()/2+1, 1, 1, 0, 1);
    esio_line_establish(h, N.x()/2+1, 0, (procid == 0 ? N.x()/2+1 : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(cbuf.data()), 2,
        "Wavenumbers in streamwise Hermitian-symmetric X direction");

    // Obtain Z wavenumbers via computing Re(1*(i*kz)/i)
    cbuf.fill(complex_t(1,0));
    diffwave::apply(0, 1, complex_t(0,-1), cbuf.data(),
            L.x(), L.z(),
            1, 1, 1, 0, 1, N.z(), N.z(), 0, N.z());
    esio_line_establish(h, N.z(), 0, (procid == 0 ? N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(cbuf.data()), 2,
        "Wavenumbers in spanwise Z direction");

    DEBUG0("Storing collocation point vectors for Fourier bases");
    ArrayXr rbuf(std::max(dN.x(), dN.z()));

    // Save collocation points in x
    for (int i = 0; i < dN.x(); ++i) {
        rbuf[i] = this->x(i);
    }
    esio_line_establish(h, dN.x(), 0, (procid == 0 ? dN.x() : 0));
    esio_line_write(h, "collocation_points_x", rbuf.data(), 0,
            "Collocation points for the dealiased, streamwise X direction");

    // Save collocation points in z
    for (int k = 0; k < dN.z(); ++k) {
        rbuf[k] = this->z(k);
    }
    esio_line_establish(h, dN.z(), 0, (procid == 0 ? dN.z() : 0));
    esio_line_write(h, "collocation_points_z", rbuf.data(), 0,
            "Collocation points for the dealiased, spanwise Z direction");
}

void
definition_grid::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading definition_grid parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    definition_grid t;
    esio_line_read(h, name_Lx, &t.L.x(), 0);
    {
        int value;
        esio_line_read(h, name_Nx, &value, 0);
        t.Nx(value);
    }
    {
        double factor;
        esio_line_read(h, name_DAFx, &factor, 0);
        t.DAFx(factor);
    }
    esio_line_read(h, name_Ly, &t.L.y(), 0);
    {
        int value;
        esio_line_read(h, name_Ny, &value, 0);
        t.Ny(value);
    }
    esio_line_read(h, name_k, &t.k, 0);
    esio_line_read(h, name_htdelta, &t.htdelta, 0);
    esio_line_read(h, name_Lz, &t.L.z(), 0);
    {
        int value;
        esio_line_read(h, name_Nz, &value, 0);
        t.Nz(value);
    }
    {
        double factor;
        esio_line_read(h, name_DAFz, &factor, 0);
        t.DAFz(factor);
    }
    this->populate(t, verbose);  // Prefer incoming to temporary
}

} // end namespace support

} // end namespace suzerain
