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
#include <suzerain/grid_definition.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace problem {


void GridDefinition::initialize_options()
{
    using ::boost::math::isnan;
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;
    using ::std::auto_ptr;
    using ::std::bind1st;
    using ::std::bind2nd;
    using ::std::mem_fun;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    { // Nx
        auto_ptr<typed_value<int> > v(value<int>(NULL));
        v->notifier(bind1st(mem_fun(&GridDefinition::Nx),this));
        if (N.x()) v->default_value(N.x());
        this->add_options()("Nx", v.release(),
                "Global logical extents in streamwise X direction");
    }

    { // DAFx
        auto_ptr<typed_value<double> > v(value<double>(NULL));
        v->notifier(bind1st(mem_fun(&GridDefinition::DAFx),this));
        if (!(isnan)(DAF.x())) v->default_value(DAF.x());
        this->add_options()("DAFx", v.release(),
                "Dealiasing factor in streamwise X direction");
    }

    { // Ny
        auto_ptr<typed_value<int> > v(value<int>(NULL));
        v->notifier(bind1st(mem_fun(&GridDefinition::Ny),this));
        if (N.y()) v->default_value(N.y());
        this->add_options()("Ny", v.release(),
                "Global logical extents in wall-normal Y direction");
    }

    { // k
        auto_ptr<typed_value<int> > v(value(&this->k));
        if (k) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "k"));
            v->default_value(k);
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"k"));
        }
        this->add_options()("k", v.release(),
                "Wall-normal B-spline order (4 indicates piecewise cubics)");
    }

    { // htdelta
        std::auto_ptr<typed_value<double> > v(value(&this->htdelta));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<double>),"htdelta"));
        if (!(isnan)(htdelta)) v->default_value(htdelta);
        this->add_options()("htdelta", v.release(),
                "Wall-normal breakpoint hyperbolic tangent stretching");
    }

    { // Nz
        auto_ptr<typed_value<int> > v(value<int>(NULL));
        v->notifier(bind1st(mem_fun(&GridDefinition::Nz),this));
        if (N.z()) v->default_value(N.z());
        this->add_options()("Nz", v.release(),
                "Global logical extents in spanwise Z direction");
    }

    { // DAFz
        auto_ptr<typed_value<double> > v(value<double>(NULL));
        v->notifier(bind1st(mem_fun(&GridDefinition::DAFz),this));
        if (!(isnan)(DAF.z())) v->default_value(DAF.z());
        this->add_options()("DAFz", v.release(),
                "Dealiasing factor in spanwise Z direction");
    }

    { // Pa
        auto_ptr<typed_value<int> > v(value(&P[0]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pa"));
        if (P[0]) v->default_value(P[0]);
        this->add_options()("Pa", v.release(),
                "Processor count in the Pa decomposition direction");
    }

    { // Pb
        auto_ptr<typed_value<int> > v(value(&P[1]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pb"));
        if (P[1]) v->default_value(P[1]);
        this->add_options()("Pb", v.release(),
                "Processor count in the Pb decomposition direction");
    }
}

GridDefinition& GridDefinition::Nx(int value) {
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

GridDefinition& GridDefinition::Ny(int value) {
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

GridDefinition& GridDefinition::Nz(int value) {
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

GridDefinition& GridDefinition::DAFx(double factor) {
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

GridDefinition& GridDefinition::DAFz(double factor) {
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

} // namespace problem

} // namespace suzerain
