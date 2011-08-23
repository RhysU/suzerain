/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * grid_definition.cpp: classes handling grid definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/validation.hpp>
#include <suzerain/exprparse.hpp>

namespace suzerain {

namespace problem {

void GridDefinition::initialize_options()
{
    using ::boost::bind;
    using ::boost::lexical_cast;
    using ::boost::math::isnan;
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;
    using ::std::auto_ptr;
    using ::std::bind1st;
    using ::std::bind2nd;
    using ::std::mem_fun;
    using ::std::ptr_fun;
    using ::std::string;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    // Used to resolve pointers to members taking strings
    GridDefinition& (GridDefinition::*f)(const std::string&) = NULL;

    { // Nx
        auto_ptr<typed_value<string> > p(value<string>(NULL));
        f = &GridDefinition::Nx;
        p->notifier(bind(f, this, _1));
        if (N.x()) p->default_value(lexical_cast<string>(N.x()));
        this->add_options()("Nx", p.release(),
                "Global logical extents in streamwise X direction");
    }

    { // DAFx
        auto_ptr<typed_value<string> > p(value<string>(NULL));
        f = &GridDefinition::DAFx;
        p->notifier(bind(f, this, _1));
        if (!(isnan)(DAF.x())) p->default_value(lexical_cast<string>(DAF.x()));
        this->add_options()("DAFx", p.release(),
                "Dealiasing factor in streamwise X direction");
    }

    { // Ny
        auto_ptr<typed_value<string> > p(value<string>(NULL));
        f = &GridDefinition::Ny;
        p->notifier(bind(f, this, _1));
        if (N.y()) p->default_value(lexical_cast<string>(N.y()));
        this->add_options()("Ny", p.release(),
                "Global logical extents in wall-normal Y direction");
    }

    { // k
        auto_ptr<typed_value<int> > p(value(&this->k));
        if (k) {
            p->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "k"));
            p->default_value(k);
        } else {
            p->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"k"));
        }
        this->add_options()("k", p.release(),
                "Wall-normal B-spline order (4 indicates piecewise cubics)");
    }

    { // htdelta
        std::auto_ptr<typed_value<double> > p(value(&this->htdelta));
        p->notifier(bind2nd(ptr_fun(ensure_nonnegative<double>),"htdelta"));
        if (!(isnan)(htdelta)) p->default_value(htdelta);
        this->add_options()("htdelta", p.release(),
                "Wall-normal breakpoint hyperbolic tangent stretching");
    }

    { // Nz
        auto_ptr<typed_value<string> > p(value<string>(NULL));
        f = &GridDefinition::Nz;
        p->notifier(bind(f, this, _1));
        if (N.z()) p->default_value(lexical_cast<string>(N.z()));
        this->add_options()("Nz", p.release(),
                "Global logical extents in spanwise Z direction");
    }

    { // DAFz
        auto_ptr<typed_value<string> > p(value<string>(NULL));
        f = &GridDefinition::DAFz;
        p->notifier(bind(f, this, _1));
        if (!(isnan)(DAF.z())) p->default_value(lexical_cast<string>(DAF.z()));
        this->add_options()("DAFz", p.release(),
                "Dealiasing factor in spanwise Z direction");
    }

    { // Pa
        auto_ptr<typed_value<int> > p(value(&P[0]));
        p->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pa"));
        if (P[0]) p->default_value(P[0]);
        this->add_options()("Pa", p.release(),
                "Processor count in the Pa decomposition direction");
    }

    { // Pb
        auto_ptr<typed_value<int> > p(value(&P[1]));
        p->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pb"));
        if (P[1]) p->default_value(P[1]);
        this->add_options()("Pb", p.release(),
                "Processor count in the Pb decomposition direction");
    }
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
    double d;
    suzerain::exprparse(value, d, "Nx");
    return Nx(static_cast<int>(d));
}

GridDefinition& GridDefinition::Ny(const std::string& value)
{
    double d;
    suzerain::exprparse(value, d, "Ny");
    return Ny(static_cast<int>(d));
}

GridDefinition& GridDefinition::Nz(const std::string& value)
{
    double d;
    suzerain::exprparse(value, d, "Nz");
    return Nz(static_cast<int>(d));
}

GridDefinition& GridDefinition::DAFx(const std::string& value)
{
    double d;
    suzerain::exprparse(value, d, "DAFx");
    return DAFx(d);
}

GridDefinition& GridDefinition::DAFz(const std::string& value)
{
    double d;
    suzerain::exprparse(value, d, "DAFz");
    return DAFz(d);
}

} // namespace problem

} // namespace suzerain
