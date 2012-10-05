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
#include <suzerain/common.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/validation.hpp>
#include <suzerain/exprparse.hpp>

namespace suzerain {

namespace problem {

template<typename T>
static void parse_option(const std::string &s,
                         T *value, void (*validator)(T, const char *),
                         const char *name)
{
#pragma warning(push,disable:2259)
    const T t = suzerain::exprparse<double>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

static const char grid_definition_description[]
        = "Mixed Fourier/B-spline computational grid definition";

GridDefinition::GridDefinition()
    : IDefinition(grid_definition_description),
      L(std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()),
      N(0, 0, 0),
      DAF(std::numeric_limits<double>::quiet_NaN(),
          1 /* Never dealiased */,
          std::numeric_limits<double>::quiet_NaN()),
      dN(0, 0, 0),
      P(0, 0),
      k(0),
      htdelta(std::numeric_limits<double>::quiet_NaN())
{
    this->initialize_options("NaN", "NaN", "NaN");  // Must match L!
}

GridDefinition::GridDefinition(const char * Lx,
                               int          Nx,
                               double       DAFx,
                               const char * Ly,
                               int          Ny,
                               int          k,
                               double       htdelta,
                               const char * Lz,
                               int          Nz,
                               double       DAFz)
    : IDefinition("Mixed Fourier/B-spline computational grid definition"),
      L(suzerain::exprparse<double>(Lx, "GridDefinition(..., Lx, ...)"),
        suzerain::exprparse<double>(Lx, "GridDefinition(..., Ly, ...)"),
        suzerain::exprparse<double>(Lx, "GridDefinition(..., Lz, ...)")),
      N(Nx, Ny, Nz),
      DAF(DAFx, 1 /* Never dealiased */, DAFz),
      dN(Nx * DAFx, Ny, Nz * DAFz),
      P(0, 0),
      k(k),
      htdelta(htdelta)
{
    this->initialize_options(Lx, Ly, Lz);
}

void GridDefinition::initialize_options(const char * default_Lx,
                                        const char * default_Ly,
                                        const char * default_Lz)
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
    using suzerain::validation::ensure_nonnegative;
    using suzerain::validation::ensure_positive;

    // Used to help resolve pointers-to-members taking strings
    GridDefinition& (GridDefinition::*f)(const std::string&) = NULL;

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;
    std::string *nullstr = NULL;

    // Lx
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<double>, _1, &L.x(),
                     &ensure_positive<double>, "Lx"));
    if (default_Lx) p->default_value(default_Lx);
    this->add_options()("Lx", p.release(),
            "Nondimensional domain length in streamwise X direction");

    // Nx
    p.reset(boost::program_options::value(nullstr));
    f = &GridDefinition::Nx;
    p->notifier(bind(f, this, _1));
    if (N.x()) p->default_value(lexical_cast<string>(N.x()));
    this->add_options()("Nx", p.release(),
            "Global logical extents in streamwise X direction");

    // DAFx
    p.reset(boost::program_options::value(nullstr));
    f = &GridDefinition::DAFx;
    p->notifier(bind(f, this, _1));
    if (!(isnan)(DAF.x())) p->default_value(lexical_cast<string>(DAF.x()));
    this->add_options()("DAFx", p.release(),
            "Dealiasing factor in streamwise X direction");

    // Ly
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<double>, _1, &L.y(),
                     &ensure_positive<double>, "Ly"));
    if (default_Ly) p->default_value(default_Ly);
    this->add_options()("Ly", p.release(),
            "Nondimensional domain length in wall normal Y direction");

    // Ny
    p.reset(boost::program_options::value(nullstr));
    f = &GridDefinition::Ny;
    p->notifier(bind(f, this, _1));
    if (N.y()) p->default_value(lexical_cast<string>(N.y()));
    this->add_options()("Ny", p.release(),
            "Global logical extents in wall-normal Y direction");

    // k
    p.reset(boost::program_options::value(nullstr));
    if (k) {
        p->notifier(bind(&parse_option<int>, _1, &k,
                         &ensure_positive<int>, "k"));
        p->default_value(lexical_cast<string>(k));
    } else {
        p->notifier(bind(&parse_option<int>, _1, &k,
                         &ensure_nonnegative<int>, "k"));
    }
    this->add_options()("k", p.release(),
            "Wall-normal B-spline order (4 indicates piecewise cubics)");

    // htdelta
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<double>, _1, &htdelta,
                        &ensure_nonnegative<double>, "htdelta"));
    if (!(isnan)(htdelta)) p->default_value(lexical_cast<string>(htdelta));
    this->add_options()("htdelta", p.release(),
            "Wall-normal breakpoint hyperbolic tangent stretching");

    // Lz
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(bind(&parse_option<double>, _1, &L.z(),
                     &ensure_positive<double>, "Lz"));
    if (default_Lz) p->default_value(default_Lz);
    this->add_options()("Lz", p.release(),
            "Nondimensional domain length in spanwise Z direction");

    // Nz
    p.reset(boost::program_options::value(nullstr));
    f = &GridDefinition::Nz;
    p->notifier(bind(f, this, _1));
    if (N.z()) p->default_value(lexical_cast<string>(N.z()));
    this->add_options()("Nz", p.release(),
            "Global logical extents in spanwise Z direction");

    // DAFz
    p.reset(boost::program_options::value(nullstr));
    f = &GridDefinition::DAFz;
    p->notifier(bind(f, this, _1));
    if (!(isnan)(DAF.z())) p->default_value(lexical_cast<string>(DAF.z()));
    this->add_options()("DAFz", p.release(),
            "Dealiasing factor in spanwise Z direction");

    // Pa
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<int>, _1, &P[0],
                     &ensure_nonnegative<int>, "Pa"));
    if (P[0]) p->default_value(lexical_cast<string>(P[0]));
    this->add_options()("Pa", p.release(),
            "Processor count in the Pa decomposition direction");

    // Pb
    p.reset(boost::program_options::value(nullstr));
    p->notifier(bind(&parse_option<int>, _1, &P[1],
                     &ensure_nonnegative<int>, "Pb"));
    if (P[1]) p->default_value(lexical_cast<string>(P[1]));
    this->add_options()("Pb", p.release(),
            "Processor count in the Pb decomposition direction");
}

GridDefinition& GridDefinition::Nx(int value)
{
    if (N.x()) {
        suzerain::validation::ensure_positive(value,"Nx");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Nx");
    }
    const_cast<int&>(N.x())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.x()) = N.x() * DAF.x();
#pragma warning(pop)
    return *this;
}

GridDefinition& GridDefinition::Ny(int value)
{
    if (N.y()) {
        suzerain::validation::ensure_positive(value,"Ny");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Ny");
    }
    const_cast<int&>(N.y())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.y()) = N.y() * DAF.y();
#pragma warning(pop)
    return *this;
}

GridDefinition& GridDefinition::Nz(int value)
{
    if (N.z()) {
        suzerain::validation::ensure_positive(value,"Nz");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Nz");
    }
    const_cast<int&>(N.z())  = value;
#pragma warning(push,disable:2259)
    const_cast<int&>(dN.z()) = N.z() * DAF.z();
#pragma warning(pop)
    return *this;
}

GridDefinition& GridDefinition::DAFx(double factor)
{
#pragma warning(push,disable:1572)
    if (DAF.x() != 0) {
#pragma warning(pop)
        suzerain::validation::ensure_positive(factor,"DAFx");
    } else {
        suzerain::validation::ensure_nonnegative(factor,"DAFx");
    }
    const_cast<double&>(DAF.x()) = factor;
#pragma warning(push,disable:2259)
    const_cast<int&   >(dN.x())  = N.x() * DAF.x();
#pragma warning(pop)
    return *this;
}

GridDefinition& GridDefinition::DAFz(double factor)
{
#pragma warning(push,disable:1572)
    if (DAF.z() != 0) {
#pragma warning(pop)
        suzerain::validation::ensure_positive(factor,"DAFz");
    } else {
        suzerain::validation::ensure_nonnegative(factor,"DAFz");
    }
    const_cast<double&>(DAF.z()) = factor;
#pragma warning(push,disable:2259)
    const_cast<int&   >(dN.z())  = N.z() * DAF.z();
#pragma warning(pop)
    return *this;
}

GridDefinition& GridDefinition::Nx(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nx(static_cast<int>(suzerain::exprparse<double>(value, "Nx")));
#pragma warning(pop)
}

GridDefinition& GridDefinition::Ny(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Ny(static_cast<int>(suzerain::exprparse<double>(value, "Ny")));
#pragma warning(pop)
}

GridDefinition& GridDefinition::Nz(const std::string& value)
{
#pragma warning(push,disable:2259)
    return Nz(static_cast<int>(suzerain::exprparse<double>(value, "Nz")));
#pragma warning(pop)
}

GridDefinition& GridDefinition::DAFx(const std::string& value)
{
    return DAFx(suzerain::exprparse<double>(value, "DAFx"));
}

GridDefinition& GridDefinition::DAFz(const std::string& value)
{
    return DAFz(suzerain::exprparse<double>(value, "DAFz"));
}

} // namespace problem

} // namespace suzerain
