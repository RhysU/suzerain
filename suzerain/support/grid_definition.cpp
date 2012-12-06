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

#include <suzerain/exprparse.hpp>
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

} // end namespace support

} // end namespace suzerain
