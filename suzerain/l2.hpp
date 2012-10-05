//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// l2.hpp: Compute L2 norms and related quantities on parallel fields
// $Id$

#ifndef SUZERAIN_L2_HPP
#define SUZERAIN_L2_HPP

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

/** Holds information on the \f$L^2\f$ norm of a scalar field */
struct L2 {
    double mean2;
    double fluctuating2;
    double total2()      const { return mean2 + fluctuating2;    };
    double total()       const { return std::sqrt(total2());     };
    double mean()        const { return std::sqrt(mean2);        };
    double fluctuating() const { return std::sqrt(fluctuating2); };
};

/**
 * Compute the \f$L^2\f$ norm of all given scalar fields.
 * See writeup/L2.tex for full details.
 */
std::vector<L2>
field_L2(const suzerain::ContiguousState<4,std::complex<double> > &state,
         const suzerain::problem::GridDefinition& grid,
         const suzerain::pencil_grid& dgrid,
         const suzerain::bsplineop& gop);

} // end namespace suzerain

#endif // SUZERAIN_L2_HPP
