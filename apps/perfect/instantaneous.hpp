//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_INSTANTANEOUS_HPP
#define SUZERAIN_PERFECT_INSTANTANEOUS_HPP

/** @file
 * Implementation of \ref instantaneous.
 */

#include <suzerain/common.hpp>

namespace suzerain {

namespace perfect {

// Forward declarations
class references;

/**
 * Storage for holding quantities at collocation points computed during
 * nonlinear operator application which are required for statistics sampling
 * purposes.  This is a subset of the information in \ref references.
 *
 * Each reference quantity is a single \e column in a column-major matrix.  This
 * facilitates a stride one operation consuming an single quantity profile.
 */
class instantaneous : private ArrayXXr
{
    typedef ArrayXXr super;

public:

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    instantaneous(): super(20, 0) {} // Magic number from count of quantities

    typedef super::Scalar      Scalar;      ///< Expose underlying scalar type
    typedef super::ColXpr      ColXpr;      ///< Mutable data for quantities
    typedef super::ConstColXpr ConstColXpr; ///< Immutable data for quantities
    using super::data;                      //   Expose raw block of memory
    using super::size;                      //   Expose size of raw memory
    using super::innerStride;               //   Expose inner stride
    using super::outerStride;               //   Expose outer stride

    /** Resize to hold data from \c Ny distinct collocation points. */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        super::setZero(NoChange, Ny);
    }

    ColXpr rho()        { return col( 0); } ///< Quantity \f$\rho         \f$
    ColXpr p()          { return col( 1); } ///< Quantity \f$p            \f$
    ColXpr p2()         { return col( 2); } ///< Quantity \f$p^2          \f$
    ColXpr u()          { return col( 3); } ///< Quantity \f$u_x          \f$
    ColXpr v()          { return col( 4); } ///< Quantity \f$u_y          \f$
    ColXpr w()          { return col( 5); } ///< Quantity \f$u_z          \f$
    ColXpr uu()         { return col( 6); } ///< Quantity \f$u_x u_x      \f$
    ColXpr uv()         { return col( 7); } ///< Quantity \f$u_x u_y      \f$
    ColXpr uw()         { return col( 8); } ///< Quantity \f$u_x u_z      \f$
    ColXpr vv()         { return col( 9); } ///< Quantity \f$u_y u_y      \f$
    ColXpr vw()         { return col(10); } ///< Quantity \f$u_y u_z      \f$
    ColXpr ww()         { return col(11); } ///< Quantity \f$u_z u_z      \f$
    ColXpr rhou()       { return col(12); } ///< Quantity \f$\rho u_x     \f$
    ColXpr rhov()       { return col(13); } ///< Quantity \f$\rho u_y     \f$
    ColXpr rhow()       { return col(14); } ///< Quantity \f$\rho u_z     \f$
    ColXpr rhoE()       { return col(15); } ///< Quantity \f$\rho E       \f$
    ColXpr rhouu()      { return col(16); } ///< Quantity \f$\rho u_x u_x \f$
    ColXpr rhovv()      { return col(17); } ///< Quantity \f$\rho u_y u_y \f$
    ColXpr rhoww()      { return col(18); } ///< Quantity \f$\rho u_z u_z \f$
    ColXpr rhoEE()      { return col(19); } ///< Quantity \f$\rho E   E   \f$

    ConstColXpr rho()        const { return col( 0); } ///< @copydoc rho()
    ConstColXpr p()          const { return col( 1); } ///< @copydoc p()
    ConstColXpr p2()         const { return col( 2); } ///< @copydoc p2()
    ConstColXpr u()          const { return col( 3); } ///< @copydoc u()
    ConstColXpr v()          const { return col( 4); } ///< @copydoc v()
    ConstColXpr w()          const { return col( 5); } ///< @copydoc w()
    ConstColXpr uu()         const { return col( 6); } ///< @copydoc uu()
    ConstColXpr uv()         const { return col( 7); } ///< @copydoc uv()
    ConstColXpr uw()         const { return col( 8); } ///< @copydoc uw()
    ConstColXpr vv()         const { return col( 9); } ///< @copydoc vv()
    ConstColXpr vw()         const { return col(10); } ///< @copydoc vw()
    ConstColXpr ww()         const { return col(11); } ///< @copydoc ww()
    ConstColXpr rhou()       const { return col(12); } ///< @copydoc rhou()
    ConstColXpr rhov()       const { return col(13); } ///< @copydoc rhov()
    ConstColXpr rhow()       const { return col(14); } ///< @copydoc rhow()
    ConstColXpr rhoE()       const { return col(15); } ///< @copydoc rhoE()
    ConstColXpr rhouu()      const { return col(16); } ///< @copydoc rhouu()
    ConstColXpr rhovv()      const { return col(17); } ///< @copydoc rhovv()
    ConstColXpr rhoww()      const { return col(18); } ///< @copydoc rhoww()
    ConstColXpr rhoEE()      const { return col(19); } ///< @copydoc rhoEE()

    /** @} */

    /** Populate all quantities from \c that instance. */
    instantaneous& operator=(const references& that);
};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_INSTANTANEOUS_HPP */
