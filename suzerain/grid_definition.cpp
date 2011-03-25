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

GridDefinition::GridDefinition(int     default_Nx,
                               double  default_DAFx,
                               int     default_Ny,
                               int     default_k,
                               double  default_htdelta,
                               int     default_Nz,
                               double  default_DAFz)
    : IDefinition("Mixed Fourier/B-spline computational grid definition"),
      N(default_Nx, default_Ny, default_Nz),
      DAF(default_DAFx, 1, default_DAFz),
      dN(default_Nx * default_DAFx, default_Ny, default_Nz * default_DAFz),
      P(0, 0),
      k(default_k),
      htdelta(default_htdelta)
{
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;
    using ::std::bind1st;
    using ::std::bind2nd;
    using ::std::mem_fun;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    this->add_options()
        ("Nx", value<int>()->default_value(N.x())
            ->notifier(bind1st(mem_fun(&GridDefinition::Nx),this)),
         "Global logical extents in streamwise X direction")
        ("DAFx", value<double>()->default_value(DAF.x())
            ->notifier(bind1st(mem_fun(&GridDefinition::DAFx),this)),
         "Dealiasing factor in streamwise X direction")
        ("Ny", value<int>()->default_value(N.y())
            ->notifier(bind1st(mem_fun(&GridDefinition::Ny),this)),
         "Global logical extents in wall-normal Y direction")
        ;

    { // k requires handling to change notifier per default_k
        std::auto_ptr<typed_value<int> > v(value(&this->k));
        if (default_k) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "k"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"k"));
        }
        v->default_value(default_k);
        this->add_options()("k", v.release(),
                "Wall-normal B-spline order (4 indicates piecewise cubics)");
    }

    { // htdelta requires handling to change notifier per default_htdelta
        std::auto_ptr<typed_value<double> > v(value(&this->htdelta));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<double>),"htdelta"));
        v->default_value(default_htdelta);
        this->add_options()("htdelta", v.release(),
                "Wall-normal breakpoint hyperbolic tangent stretching");
    }

    this->add_options()
        ("Nz", value<int>()->default_value(N.z())
            ->notifier(bind1st(mem_fun(&GridDefinition::Nz),this)),
         "Global logical extents in spanwise Z direction")
        ("DAFz", value<double>()->default_value(DAF.z())
            ->notifier(bind1st(mem_fun(&GridDefinition::DAFz),this)),
         "Dealiasing factor in spanwise Z direction")
        ("Pa", value<int>(&P[0])->default_value(P[0])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pa")),
            "Processor count in the Pa decomposition direction")
        ("Pb", value<int>(&P[1])->default_value(P[1])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pb")),
            "Processor count in the Pb decomposition direction")
        ;
}

GridDefinition& GridDefinition::Nx(int value) {
    if (N.x()) {
        suzerain::validation::ensure_positive(value,"Nx");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Nx");
    }
    const_cast<int&>(N.x())  = value;
    const_cast<int&>(dN.x()) = N.x() * DAF.x();
    return *this;
}

GridDefinition& GridDefinition::Ny(int value) {
    if (N.y()) {
        suzerain::validation::ensure_positive(value,"Ny");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Ny");
    }
    const_cast<int&>(N.y())  = value;
    const_cast<int&>(dN.y()) = N.y() * DAF.y();
    return *this;
}

GridDefinition& GridDefinition::Nz(int value) {
    if (N.z()) {
        suzerain::validation::ensure_positive(value,"Nz");
    } else {
        suzerain::validation::ensure_nonnegative(value,"Nz");
    }
    const_cast<int&>(N.z())  = value;
    const_cast<int&>(dN.z()) = N.z() * DAF.z();
    return *this;
}

GridDefinition& GridDefinition::DAFx(double factor) {
    if (DAF.x()) {
        suzerain::validation::ensure_positive(factor,"DAFx");
    } else {
        suzerain::validation::ensure_nonnegative(factor,"DAFx");
    }
    const_cast<double&>(DAF.x()) = factor;
    const_cast<int&   >(dN.x())  = N.x() * DAF.x();
    return *this;
}

GridDefinition& GridDefinition::DAFz(double factor) {
    if (DAF.z()) {
        suzerain::validation::ensure_positive(factor,"DAFz");
    } else {
        suzerain::validation::ensure_nonnegative(factor,"DAFz");
    }
    const_cast<double&>(DAF.z()) = factor;
    const_cast<int&   >(dN.z())  = N.z() * DAF.z();
    return *this;
}

} // namespace problem

} // namespace suzerain
