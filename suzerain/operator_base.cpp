//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
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
// operator_base.cpp: useful base class for operator implementation work
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/operator_base.hpp>

namespace suzerain {

OperatorBase::OperatorBase(
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop)
    : one_over_delta_x(grid.N.x() /* !dN.x() */ / grid.L.x()),
      lambda1_x(boost::math::constants::pi<real_t>() * one_over_delta_x),
      lambda2_x(lambda1_x * lambda1_x),
      one_over_delta_z(grid.N.z() /* !dN.z() */ / grid.L.z()),
      lambda1_z(boost::math::constants::pi<real_t>() * one_over_delta_z),
      lambda2_z(lambda1_z * lambda1_z),
      grid(grid),
      dgrid(dgrid),
      bop(bop),
      y_(boost::extents[boost::multi_array_types::extent_range(
              dgrid.local_physical_start.y(),
              dgrid.local_physical_end.y())]),
      one_over_delta_y_(
              boost::extents[boost::multi_array_types::extent_range(
                    dgrid.local_physical_start.y(),
                    dgrid.local_physical_end.y())]),
      lambda1_y_(
              boost::extents[boost::multi_array_types::extent_range(
                    dgrid.local_physical_start.y(),
                    dgrid.local_physical_end.y())]),
      lambda2_y_(
              boost::extents[boost::multi_array_types::extent_range(
                    dgrid.local_physical_start.y(),
                    dgrid.local_physical_end.y())])
{
    const real_t pi = boost::math::constants::pi<real_t>();

    // Compute the B-spline-dependent correction factor to obtain
    // good maximum eigenvalue estimates given wall-normal inhomogeneity.
    // See model document for definitions of C^{(1)} and C^{(2)}.
    real_t C1, Clow1, Chigh1;
    suzerain_bspline_htstretch2_evdeltascale(
            1, b.k(), grid.htdelta, b.n(), &C1, &Clow1, &Chigh1);
    real_t C2, Clow2, Chigh2;
    suzerain_bspline_htstretch2_evdeltascale(
            2, b.k(), grid.htdelta, b.n(), &C2, &Clow2, &Chigh2);

    // In practice, directly using C^{(1)} and C^{(2)} is too aggressive
    // as it requires using inconsistent safety factors when computing
    // convectively- or diffusively-limited test problems.  Softening
    // via taking sqrt(C^{(i)}) seems to permit a single safety factor
    // to address both convective and diffusive stability restrictions.
    using std::sqrt;
    C1 = sqrt(C1); Clow1 = sqrt(Clow1); Chigh1 = sqrt(Chigh1);
    C2 = sqrt(C2); Clow2 = sqrt(Clow2); Chigh2 = sqrt(Chigh2);

    // Compute collocation point-based information local to this rank
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Collocation point locations
        y_[j] = b.collocation_point(j);

        // Inverse spacing to next adjacent collocation point
        one_over_delta_y_[j] = 1 / b.spacing_collocation_point(j);

        // Eigenvalue estimates using collocation point spacing
        const real_t inverse_spacing = one_over_delta_y_[j];

        // Estimating wall-normal first derivative eigenvalue magnitudes
        // Use Clow1 rather than C1 as it is always slightly conservative
        lambda1_y_[j]  = pi * inverse_spacing / Clow1;

        // Estimating wall-normal second derivative eigenvalue magnitudes
        // Use Clow2 rather than C2 as it is always slightly conservative
        lambda2_y_[j]  = pi * inverse_spacing / Clow2;
        lambda2_y_[j] *= lambda2_y_[j];
    }
}

OperatorBase::~OperatorBase()
{
    // NOP
}

} // namespace suzerain
