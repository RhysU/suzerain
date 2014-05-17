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

#ifndef SUZERAIN_PERFECT_IMPLICITS_HPP
#define SUZERAIN_PERFECT_IMPLICITS_HPP

/** @file
 * Implementation of \ref implicits.
 */

#include <suzerain/common.hpp>

#include <suzerain/treatment_constraint.hpp>

namespace suzerain {

namespace perfect {

/**
 * Storage for holding implicit quantities at collocation points
 * updated by \ref constraint::treatment usage.
 *
 * Each quantity is a single \e column in a column-major matrix.  This
 * facilitates a stride one operation consuming an single quantity profile.
 */
class implicits
    : public virtual constraint::treatment::outputs
    , private ArrayXXr
{
    typedef ArrayXXr super;

public:

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
        SrhoE,        //   FIXME Redmine #3089
        Srhou,        //   FIXME Redmine #3089
        Srhov,        //   FIXME Redmine #3089
        Srhow,        //   FIXME Redmine #3089
        Srho,         //   FIXME Redmine #3089
        Srhou_dot_u,  //   FIXME Redmine #3089
        fx,           ///< @copydoc constraint::treatment::outputs::fx
        fy,           ///< @copydoc constraint::treatment::outputs::fy
        fz,           ///< @copydoc constraint::treatment::outputs::fz
        f_dot_u,      ///< @copydoc constraint::treatment::outputs::f_dot_u
        qb,           ///< @copydoc constraint::treatment::outputs::qb
        CrhoE,        ///< @copydoc constraint::treatment::outputs::CrhoE
        C2rhoE,       ///< @copydoc constraint::treatment::outputs::C2rhoE
        Crhou,        ///< @copydoc constraint::treatment::outputs::Crhou
        C2rhou,       ///< @copydoc constraint::treatment::outputs::C2rhou
        Crhov,        ///< @copydoc constraint::treatment::outputs::Crhov
        C2rhov,       ///< @copydoc constraint::treatment::outputs::C2rhov
        Crhow,        ///< @copydoc constraint::treatment::outputs::Crhow
        C2rhow,       ///< @copydoc constraint::treatment::outputs::C2rhow
        Crho,         ///< @copydoc constraint::treatment::outputs::Crho
        C2rho,        ///< @copydoc constraint::treatment::outputs::C2rho
        Crhou_dot_u,  ///< @copydoc constraint::treatment::outputs::Crhou_dot_u
        C2rhou_dot_u, ///< @copydoc constraint::treatment::outputs::C2rhou_dot_u
        count         ///< Sentry indicating how many quantities are tracked
    }; };

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    implicits() : super(0, static_cast<super::Index>(q::count))
    {
    }

    /** Resize to hold data from \c Ny distinct collocation points. */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        super::setZero(Ny, NoChange);
    }

    ColXpr      SrhoE()              { return col(q::SrhoE       ); } ///< @copydoc q::SrhoE       
    ColXpr      Srhou()              { return col(q::Srhou       ); } ///< @copydoc q::Srhou       
    ColXpr      Srhov()              { return col(q::Srhov       ); } ///< @copydoc q::Srhov       
    ColXpr      Srhow()              { return col(q::Srhow       ); } ///< @copydoc q::Srhow       
    ColXpr      Srho()               { return col(q::Srho        ); } ///< @copydoc q::Srho        
    ColXpr      Srhou_dot_u()        { return col(q::Srhou_dot_u ); } ///< @copydoc q::Srhou_dot_u 
    ColXpr      fx()                 { return col(q::fx          ); } ///< @copydoc q::fx          
    ColXpr      fy()                 { return col(q::fy          ); } ///< @copydoc q::fy          
    ColXpr      fz()                 { return col(q::fz          ); } ///< @copydoc q::fz          
    ColXpr      f_dot_u()            { return col(q::f_dot_u     ); } ///< @copydoc q::f_dot_u     
    ColXpr      qb()                 { return col(q::qb          ); } ///< @copydoc q::qb          
    ColXpr      CrhoE()              { return col(q::CrhoE       ); } ///< @copydoc q::CrhoE       
    ColXpr      C2rhoE()             { return col(q::C2rhoE      ); } ///< @copydoc q::C2rhoE      
    ColXpr      Crhou()              { return col(q::Crhou       ); } ///< @copydoc q::Crhou       
    ColXpr      C2rhou()             { return col(q::C2rhou      ); } ///< @copydoc q::C2rhou      
    ColXpr      Crhov()              { return col(q::Crhov       ); } ///< @copydoc q::Crhov       
    ColXpr      C2rhov()             { return col(q::C2rhov      ); } ///< @copydoc q::C2rhov      
    ColXpr      Crhow()              { return col(q::Crhow       ); } ///< @copydoc q::Crhow       
    ColXpr      C2rhow()             { return col(q::C2rhow      ); } ///< @copydoc q::C2rhow      
    ColXpr      Crho()               { return col(q::Crho        ); } ///< @copydoc q::Crho        
    ColXpr      C2rho()              { return col(q::C2rho       ); } ///< @copydoc q::C2rho       
    ColXpr      Crhou_dot_u()        { return col(q::Crhou_dot_u ); } ///< @copydoc q::Crhou_dot_u 
    ColXpr      C2rhou_dot_u()       { return col(q::C2rhou_dot_u); } ///< @copydoc q::C2rhou_dot_u

    ConstColXpr SrhoE()        const { return col(q::SrhoE       ); } ///< @copydoc q::SrhoE       
    ConstColXpr Srhou()        const { return col(q::Srhou       ); } ///< @copydoc q::Srhou       
    ConstColXpr Srhov()        const { return col(q::Srhov       ); } ///< @copydoc q::Srhov       
    ConstColXpr Srhow()        const { return col(q::Srhow       ); } ///< @copydoc q::Srhow       
    ConstColXpr Srho()         const { return col(q::Srho        ); } ///< @copydoc q::Srho        
    ConstColXpr Srhou_dot_u()  const { return col(q::Srhou_dot_u ); } ///< @copydoc q::Srhou_dot_u 
    ConstColXpr fx()           const { return col(q::fx          ); } ///< @copydoc q::fx          
    ConstColXpr fy()           const { return col(q::fy          ); } ///< @copydoc q::fy          
    ConstColXpr fz()           const { return col(q::fz          ); } ///< @copydoc q::fz          
    ConstColXpr f_dot_u()      const { return col(q::f_dot_u     ); } ///< @copydoc q::f_dot_u     
    ConstColXpr qb()           const { return col(q::qb          ); } ///< @copydoc q::qb          
    ConstColXpr CrhoE()        const { return col(q::CrhoE       ); } ///< @copydoc q::CrhoE       
    ConstColXpr C2rhoE()       const { return col(q::C2rhoE      ); } ///< @copydoc q::C2rhoE      
    ConstColXpr Crhou()        const { return col(q::Crhou       ); } ///< @copydoc q::Crhou       
    ConstColXpr C2rhou()       const { return col(q::C2rhou      ); } ///< @copydoc q::C2rhou      
    ConstColXpr Crhov()        const { return col(q::Crhov       ); } ///< @copydoc q::Crhov       
    ConstColXpr C2rhov()       const { return col(q::C2rhov      ); } ///< @copydoc q::C2rhov      
    ConstColXpr Crhow()        const { return col(q::Crhow       ); } ///< @copydoc q::Crhow       
    ConstColXpr C2rhow()       const { return col(q::C2rhow      ); } ///< @copydoc q::C2rhow      
    ConstColXpr Crho()         const { return col(q::Crho        ); } ///< @copydoc q::Crho        
    ConstColXpr C2rho()        const { return col(q::C2rho       ); } ///< @copydoc q::C2rho       
    ConstColXpr Crhou_dot_u()  const { return col(q::Crhou_dot_u ); } ///< @copydoc q::Crhou_dot_u 
    ConstColXpr C2rhou_dot_u() const { return col(q::C2rhou_dot_u); } ///< @copydoc q::C2rhou_dot_u
};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_IMPLICITS_HPP */
