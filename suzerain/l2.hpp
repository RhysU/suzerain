//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_L2_HPP
#define SUZERAIN_L2_HPP

/** @file
 * Compute \f$L^2\f$ norms and related quantities on distributed data.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

namespace suzerain {

/**
 * Holds information on the \f$L^2_{xyz}\f$ norm of a scalar field
 * or inner product of two scalar fields.
 */
struct field_L2xyz {
    real_t mean2;
    real_t fluctuating2;
    real_t total2()      const { return mean2 + fluctuating2;    };
    real_t total()       const { return std::sqrt(total2());     };
    real_t mean()        const { return std::sqrt(mean2);        };
    real_t fluctuating() const { return std::sqrt(fluctuating2); };
};

/**
 * Compute the \f$L^2_{xyz}\f$ norm of all given scalar fields.
 * See writeup/L2.tex for full details.
 */
std::vector<field_L2xyz>
compute_field_L2xyz(
        const contiguous_state<4,complex_t> &state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& gop);

} // end namespace suzerain

#endif // SUZERAIN_L2_HPP
