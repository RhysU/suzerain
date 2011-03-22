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
                               int     default_Nz,
                               double  default_DAFz)
    : IDefinition("Mixed Fourier/B-spline computational grid definition"),
      // global_extents below
      N(default_Nx, default_Ny, default_Nz),
      dN(default_Nx * default_DAFx, default_Ny, default_Nz * default_DAFz),
      k(default_k),
      P(0, 0)
{
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;
    using ::std::auto_ptr;
    using ::std::bind1st;
    using ::std::bind2nd;
    using ::std::mem_fun;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<double,const char*,void>
        ptr_fun_ensure_positive_double(ensure_positive<double>);
    ::std::pointer_to_binary_function<double,const char*,void>
        ptr_fun_ensure_nonnegative_double(ensure_nonnegative<double>);

    // Complicated add_options() calls done to allow changing the validation
    // routine in use when the default provided value is zero.  Zero is
    // generally used a NOP value by some client code.

    { // Nx
        auto_ptr<typed_value<int> > v(value(&this->N.x()));
        if (default_Nx) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "Nx"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Nx"));
        }
        v->default_value(default_Nx);
        this->add_options()("Nx", v.release(),
                "Global logical extents in streamwise X direction");
    }

    { // DAFx
        auto_ptr<typed_value<double> > v(value<double>());
        GridDefinition& (GridDefinition::*callback)(double)
            = &GridDefinition::DAFx;
        v->notifier(bind1st(mem_fun(callback),this));
        v->default_value(default_DAFx);
        this->add_options()("DAFx", v.release(),
            "Dealiasing factor in streamwise X direction");
    }

    { // Ny
        auto_ptr<typed_value<int> > v(value(&this->N.y()));
        if (default_Ny) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "Ny"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Ny"));
        }
        v->default_value(default_Ny);
        this->add_options()("Ny", v.release(),
                "Global logical extents in wall-normal Y direction");
    }

    { // k
        auto_ptr<typed_value<int> > v(value(&this->k));
        if (default_k) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "k"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"k"));
        }
        v->default_value(default_k);
        this->add_options()("k", v.release(),
                "B-spline basis order where k = 4 indicates piecewise cubics");
    }

    { // Nz
        auto_ptr<typed_value<int> > v(value(&this->N.z()));
        if (default_Nz) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<int>),   "Nz"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Nz"));
        }
        v->default_value(default_Nz);
        this->add_options()("Nz", v.release(),
                "Global logical extents in spanwise Z direction");
    }

    { // DAFz
        auto_ptr<typed_value<double> > v(value<double>());
        GridDefinition& (GridDefinition::*callback)(double)
            = &GridDefinition::DAFz;
        v->notifier(bind1st(mem_fun(callback),this));
        v->default_value(default_DAFz);
        this->add_options()("DAFz", v.release(),
            "Dealiasing factor in spanwise Z direction");
    }

    { // Pa
        auto_ptr<typed_value<int> > v(value(&this->P[0]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pa"));
        v->default_value(0);
        this->add_options()("Pa", v.release(),
            "Processor count in the Pa decomposition direction");
    }

    { // Pb
        auto_ptr<typed_value<int> > v(value(&this->P[0]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"Pb"));
        v->default_value(0);
        this->add_options()("Pb", v.release(),
            "Processor count in the Pb decomposition direction");
    }
}

GridDefinition& GridDefinition::DAFx(double factor) {
    suzerain::validation::ensure_nonnegative(factor,"DAFx");
    dN.x() = N.x() * factor;
    return *this;
}

GridDefinition& GridDefinition::DAFz(double factor) {
    suzerain::validation::ensure_nonnegative(factor,"DAFz");
    dN.z() = N.z() * factor;
    return *this;
}

} // namespace problem

} // namespace suzerain
