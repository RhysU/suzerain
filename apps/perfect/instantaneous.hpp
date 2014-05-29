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

#include <suzerain/treatment_constraint.hpp>

#include "slowgrowth.hpp"

namespace suzerain {

// Forward declarations
class bsplineop;
class samples;

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
class instantaneous
    : public virtual constraint::treatment::inputs
    , public virtual slowgrowth::physical_cons
    , public virtual slowgrowth::physical_rqq
    , protected ArrayXXr
{
    typedef ArrayXXr super;

public:

    /** Must call \ref set_zero prior to use. */
    instantaneous();

    /** Virtual destructor to permit use as a base class. */
    virtual ~instantaneous();

    /** Resize to hold data from \c Ny distinct collocation points. */
    virtual void set_zero(int Ny);

    typedef super::ColXpr      ColXpr;      ///< Mutable data for quantities
    typedef super::ConstColXpr ConstColXpr; ///< Immutable data for quantities
    typedef super::Scalar      Scalar;      ///< Expose underlying scalar type
    using super::col;                       //   Expose each quantity by index
    using super::cols;                      //   Expose number of columns
    using super::data;                      //   Expose raw block of memory
    using super::innerStride;               //   Expose inner stride
    using super::operator/=;                //   Permit division-based scaling
    using super::operator*=;                //   Permit multiplicative scaling
    using super::outerStride;               //   Expose outer stride
    using super::rows;                      //   Expose number of rows
    using super::size;                      //   Expose size of raw memory

    /** Logical indices for each scalar quantity that will be stored. */
    struct q { enum {
        rho,        ///< Quantity \f$\rho         \f$
        p,          ///< Quantity \f$p            \f$
        p2,         ///< Quantity \f$p^2          \f$
        u,          ///< Quantity \f$u_x          \f$
        v,          ///< Quantity \f$u_y          \f$
        w,          ///< Quantity \f$u_z          \f$
        uu,         ///< Quantity \f$u_x u_x      \f$
        uv,         ///< Quantity \f$u_x u_y      \f$
        uw,         ///< Quantity \f$u_x u_z      \f$
        vv,         ///< Quantity \f$u_y u_y      \f$
        vw,         ///< Quantity \f$u_y u_z      \f$
        ww,         ///< Quantity \f$u_z u_z      \f$
        rhou,       ///< Quantity \f$\rho u_x     \f$
        rhov,       ///< Quantity \f$\rho u_y     \f$
        rhow,       ///< Quantity \f$\rho u_z     \f$
        rhoE,       ///< Quantity \f$\rho E       \f$
        rhouu,      ///< Quantity \f$\rho u_x u_x \f$
        rhovv,      ///< Quantity \f$\rho u_y u_y \f$
        rhoww,      ///< Quantity \f$\rho u_z u_z \f$
        rhoEE,      ///< Quantity \f$\rho E   E   \f$
        count       ///< Sentry indicating how many quantities are tracked
    }; };

    ColXpr      rho()              { return col(q::rho  ); } ///< @copydoc q::rho  
    ColXpr      p()                { return col(q::p    ); } ///< @copydoc q::p    
    ColXpr      p2()               { return col(q::p2   ); } ///< @copydoc q::p2   
    ColXpr      u()                { return col(q::u    ); } ///< @copydoc q::u    
    ColXpr      v()                { return col(q::v    ); } ///< @copydoc q::v    
    ColXpr      w()                { return col(q::w    ); } ///< @copydoc q::w    
    ColXpr      uu()               { return col(q::uu   ); } ///< @copydoc q::uu   
    ColXpr      uv()               { return col(q::uv   ); } ///< @copydoc q::uv   
    ColXpr      uw()               { return col(q::uw   ); } ///< @copydoc q::uw   
    ColXpr      vv()               { return col(q::vv   ); } ///< @copydoc q::vv   
    ColXpr      vw()               { return col(q::vw   ); } ///< @copydoc q::vw   
    ColXpr      ww()               { return col(q::ww   ); } ///< @copydoc q::ww   
    ColXpr      rhou()             { return col(q::rhou ); } ///< @copydoc q::rhou 
    ColXpr      rhov()             { return col(q::rhov ); } ///< @copydoc q::rhov 
    ColXpr      rhow()             { return col(q::rhow ); } ///< @copydoc q::rhow 
    ColXpr      rhoE()             { return col(q::rhoE ); } ///< @copydoc q::rhoE 
    ColXpr      rhouu()            { return col(q::rhouu); } ///< @copydoc q::rhouu
    ColXpr      rhovv()            { return col(q::rhovv); } ///< @copydoc q::rhovv
    ColXpr      rhoww()            { return col(q::rhoww); } ///< @copydoc q::rhoww
    ColXpr      rhoEE()            { return col(q::rhoEE); } ///< @copydoc q::rhoEE

    ConstColXpr rho()        const { return col(q::rho  ); } ///< @copydoc q::rho   
    ConstColXpr p()          const { return col(q::p    ); } ///< @copydoc q::p     
    ConstColXpr p2()         const { return col(q::p2   ); } ///< @copydoc q::p2    
    ConstColXpr u()          const { return col(q::u    ); } ///< @copydoc q::u     
    ConstColXpr v()          const { return col(q::v    ); } ///< @copydoc q::v     
    ConstColXpr w()          const { return col(q::w    ); } ///< @copydoc q::w     
    ConstColXpr uu()         const { return col(q::uu   ); } ///< @copydoc q::uu    
    ConstColXpr uv()         const { return col(q::uv   ); } ///< @copydoc q::uv    
    ConstColXpr uw()         const { return col(q::uw   ); } ///< @copydoc q::uw    
    ConstColXpr vv()         const { return col(q::vv   ); } ///< @copydoc q::vv    
    ConstColXpr vw()         const { return col(q::vw   ); } ///< @copydoc q::vw    
    ConstColXpr ww()         const { return col(q::ww   ); } ///< @copydoc q::ww    
    ConstColXpr rhou()       const { return col(q::rhou ); } ///< @copydoc q::rhou  
    ConstColXpr rhov()       const { return col(q::rhov ); } ///< @copydoc q::rhov  
    ConstColXpr rhow()       const { return col(q::rhow ); } ///< @copydoc q::rhow  
    ConstColXpr rhoE()       const { return col(q::rhoE ); } ///< @copydoc q::rhoE  
    ConstColXpr rhouu()      const { return col(q::rhouu); } ///< @copydoc q::rhouu 
    ConstColXpr rhovv()      const { return col(q::rhovv); } ///< @copydoc q::rhovv 
    ConstColXpr rhoww()      const { return col(q::rhoww); } ///< @copydoc q::rhoww 
    ConstColXpr rhoEE()      const { return col(q::rhoEE); } ///< @copydoc q::rhoEE 

    /** @} */

    /** Assign all quantities from \c that instance. */
    instantaneous& operator=(const references& that);

    /** Copy all quantities from \c that instance using \c cop. */
    instantaneous& copy_from(const bsplineop& cop,
                             const samples& that);
};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_INSTANTANEOUS_HPP */
