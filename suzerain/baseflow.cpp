//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc baseflow.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/baseflow.hpp>

#include <gsl/gsl_poly.h>

#include <suzerain/bspline.hpp>
#include <suzerain/radial_nozzle.h>

namespace suzerain {

baseflow_interface::baseflow_interface()
{
}

baseflow_interface::~baseflow_interface()
{
}

baseflow_uniform::baseflow_uniform()
    : x ()
{
}

void
baseflow_uniform::conserved(
        const real_t      y,
        real_t*        base,
        real_t*      dybase,
        real_t*      dxbase) const
{
    SUZERAIN_UNUSED(y);
    const int nstate = x.size() - 1;
    memcpy(base,   x.data(), nstate*sizeof(real_t));
    memset(dybase, 0,        nstate*sizeof(real_t));
    memset(dxbase, 0,        nstate*sizeof(real_t));
}

void
baseflow_uniform::pressure(
        const real_t y,
        real_t&      P,
        real_t&      dyP,
        real_t&      dxP) const
{
    SUZERAIN_UNUSED(y);
    P   = x.tail<1>()[0];
    dyP = 0;
    dxP = 0;
}

baseflow_polynomial::baseflow_polynomial()
    : x ()
    , dx()
{
}

void
baseflow_polynomial::conserved(
        const real_t      y,
        real_t*        base,
        real_t*      dybase,
        real_t*      dxbase) const
{
    // Consistent number of state variables + pressure must be provided.
    assert(x.cols() > 0);
    assert(x.cols() == dx.cols());
    const int nstate = x.cols() - 1;

    // Compute baseflow and y-derivative
    if (x.rows()) {

        double res[2];
        for (int j = 0; j < nstate; ++j) {
            gsl_poly_eval_derivs(x.col(j).data(), x.rows(), y,
                                 &res[0], sizeof(res) / sizeof(res[0]));
            base  [j] = res[0];
            dybase[j] = res[1];
        }

    } else {

        using namespace std;
        memset(base,   0, nstate*sizeof(real_t));
        memset(dybase, 0, nstate*sizeof(real_t));

    }

    // Compute x-derivative of baseflow
    if (dx.rows()) {

        for (int j = 0; j < nstate; j++) {
            dxbase[j] = gsl_poly_eval(dx.col(j).data(), dx.rows(), y);
        }

    } else {

        using namespace std;
        memset(dxbase, 0, nstate*sizeof(real_t));

    }
}

void
baseflow_polynomial::pressure(
        const real_t   y,
        real_t&        P,
        real_t&      dyP,
        real_t&      dxP) const
{
    // Ensure consistent number of state variables + pressure
    assert(x.cols() > 0);
    assert(x.cols() == dx.cols());

    // Compute pressure and y-derivative
    double res[2] = { 0, 0 };
    if (x.rows()) {
        gsl_poly_eval_derivs(x.rightCols<1>().data(), x.rows(), y,
                             &res[0], sizeof(res) / sizeof(res[0]));
    }
    P   = res[0];
    dyP = res[1];

    // Compute x-derivative of baseflow
    dxP = dx.rows()
        ? gsl_poly_eval(dx.rightCols<1>().data(), dx.rows(), y)
        : 0;
}

baseflow_map::baseflow_map()
{
}

void
baseflow_map::conserved(
        const real_t      y,
        real_t*        base,
        real_t*      dybase,
        real_t*      dxbase) const
{
    // FIXME Implement
}

void
baseflow_map::pressure(
        const real_t   y,
        real_t&        P,
        real_t&      dyP,
        real_t&      dxP) const
{
    // FIXME Implement
}

} // namespace suzerain
