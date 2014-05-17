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

#ifndef SUZERAIN_PERFECT_REFERENCES_HPP
#define SUZERAIN_PERFECT_REFERENCES_HPP

/** @file
 * Implementation of \ref references.
 */

#include <suzerain/common.hpp>

// Forward declarations
struct suzerain_rholut_imexop_ref;
struct suzerain_rholut_imexop_refld;

namespace suzerain {

namespace perfect {

/**
 * Storage for holding quantities at collocation points computed during
 * nonlinear operator application which either are required for hybrid linear
 * operator application or for statistics sampling purposes.
 *
 * Each reference quantity is a single \e row in a column-major matrix.  This
 * facilitates a stride one operation loading or writing all reference
 * quantities for a single wall-normal location.
 *
 * @see rholut_imexop.h for details on linearization and the associated
 * reference quantities.
 */
class references : private ArrayXXr
{
    typedef ArrayXXr super;

public:

    typedef super::ConstRowXpr ConstRowXpr; ///< Immutable data for quantities
    typedef super::RowXpr      RowXpr;      ///< Mutable data for quantities
    typedef super::Scalar      Scalar;      ///< Expose underlying scalar type
    using super::cols;                      //   Expose number of columns
    using super::data;                      //   Expose raw block of memory
    using super::innerStride;               //   Expose inner stride
    using super::outerStride;               //   Expose outer stride
    using super::row;                       //   Expose each quantity by index
    using super::rows;                      //   Expose number of rows
    using super::size;                      //   Expose size of raw memory

    /** Logical indices for each scalar quantity that will be stored. */
    struct q { enum {
        rho,        ///< Reference \f$\rho                  \f$
        p,          ///< Reference \f$p                     \f$
        p2,         ///< Quantity  \f$p^2                   \f$
        T,          ///< Reference \f$T                     \f$
        a,          ///< Reference \f$a = \sqrt{T}          \f$
        u,          ///< Reference \f$C^{u_x}               \f$
        v,          ///< Reference \f$C^{u_y}               \f$
        w,          ///< Reference \f$C^{u_z}               \f$
        u2,         ///< Quantity  \f$\vec{u}^2             \f$
        uu,         ///< Reference \f$C^{u_x u_x}           \f$
        uv,         ///< Reference \f$C^{u_x u_y}           \f$
        uw,         ///< Reference \f$C^{u_x u_z}           \f$
        vv,         ///< Reference \f$C^{u_y u_y}           \f$
        vw,         ///< Reference \f$C^{u_y u_z}           \f$
        ww,         ///< Reference \f$C^{u_z u_z}           \f$
        nu,         ///< Reference \f$C^{\nu}               \f$
        nu_u,       ///< Reference \f$C^{\nu u_x}           \f$
        nu_v,       ///< Reference \f$C^{\nu u_y}           \f$
        nu_w,       ///< Reference \f$C^{\nu u_z}           \f$
        nu_u2,      ///< Reference \f$C^{\nu u^2}           \f$
        nu_uu,      ///< Reference \f$C^{\nu u_x u_x}       \f$
        nu_uv,      ///< Reference \f$C^{\nu u_x u_y}       \f$
        nu_uw,      ///< Reference \f$C^{\nu u_x u_z}       \f$
        nu_vv,      ///< Reference \f$C^{\nu u_y u_y}       \f$
        nu_vw,      ///< Reference \f$C^{\nu u_y u_z}       \f$
        nu_ww,      ///< Reference \f$C^{\nu u_z u_z}       \f$
        ex_gradrho, ///< Reference \f$C^{e_x}_{\nabla\rho}  \f$
        ey_gradrho, ///< Reference \f$C^{e_y}_{\nabla\rho}  \f$
        ez_gradrho, ///< Reference \f$C^{e_z}_{\nabla\rho}  \f$
        e_divm,     ///< Reference \f$C^{e}_{\nabla\cdot{}m}\f$
        e_deltarho, ///< Reference \f$C^{e}_{\Delta\rho}    \f$
        rhou,       ///< Quantity  \f$\rho u_x              \f$
        rhov,       ///< Quantity  \f$\rho u_y              \f$
        rhow,       ///< Quantity  \f$\rho u_z              \f$
        rhoE,       ///< Quantity  \f$\rho E                \f$
        rhouu,      ///< Quantity  \f$\rho u_x u_x          \f$
        rhovv,      ///< Quantity  \f$\rho u_y u_y          \f$
        rhoww,      ///< Quantity  \f$\rho u_z u_z          \f$
        rhoEE,      ///< Quantity  \f$\rho E   E            \f$
        count       ///< Sentry indicating how many quantities are tracked
    }; };

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    references() : super(static_cast<super::Index>(q::count), 0)
    {
    }

    /** Resize to hold data from \c Ny distinct collocation points. */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        super::setZero(NoChange, Ny);
    }

    RowXpr      rho()              { return row(q::rho       ); } ///< @copydoc q::rho
    RowXpr      p()                { return row(q::p         ); } ///< @copydoc q::p
    RowXpr      p2()               { return row(q::p2        ); } ///< @copydoc q::p2
    RowXpr      T()                { return row(q::T         ); } ///< @copydoc q::T
    RowXpr      a()                { return row(q::a         ); } ///< @copydoc q::a
    RowXpr      u()                { return row(q::u         ); } ///< @copydoc q::u
    RowXpr      v()                { return row(q::v         ); } ///< @copydoc q::v
    RowXpr      w()                { return row(q::w         ); } ///< @copydoc q::w
    RowXpr      u2()               { return row(q::u2        ); } ///< @copydoc q::u2
    RowXpr      uu()               { return row(q::uu        ); } ///< @copydoc q::uu
    RowXpr      uv()               { return row(q::uv        ); } ///< @copydoc q::uv
    RowXpr      uw()               { return row(q::uw        ); } ///< @copydoc q::uw
    RowXpr      vv()               { return row(q::vv        ); } ///< @copydoc q::vv
    RowXpr      vw()               { return row(q::vw        ); } ///< @copydoc q::vw
    RowXpr      ww()               { return row(q::ww        ); } ///< @copydoc q::ww
    RowXpr      nu()               { return row(q::nu        ); } ///< @copydoc q::nu
    RowXpr      nu_u()             { return row(q::nu_u      ); } ///< @copydoc q::nu_u
    RowXpr      nu_v()             { return row(q::nu_v      ); } ///< @copydoc q::nu_v
    RowXpr      nu_w()             { return row(q::nu_w      ); } ///< @copydoc q::nu_w
    RowXpr      nu_u2()            { return row(q::nu_u2     ); } ///< @copydoc q::nu_u2
    RowXpr      nu_uu()            { return row(q::nu_uu     ); } ///< @copydoc q::nu_uu
    RowXpr      nu_uv()            { return row(q::nu_uv     ); } ///< @copydoc q::nu_uv
    RowXpr      nu_uw()            { return row(q::nu_uw     ); } ///< @copydoc q::nu_uw
    RowXpr      nu_vv()            { return row(q::nu_vv     ); } ///< @copydoc q::nu_vv
    RowXpr      nu_vw()            { return row(q::nu_vw     ); } ///< @copydoc q::nu_vw
    RowXpr      nu_ww()            { return row(q::nu_ww     ); } ///< @copydoc q::nu_ww
    RowXpr      ex_gradrho()       { return row(q::ex_gradrho); } ///< @copydoc q::ex_gradrho
    RowXpr      ey_gradrho()       { return row(q::ey_gradrho); } ///< @copydoc q::ey_gradrho
    RowXpr      ez_gradrho()       { return row(q::ez_gradrho); } ///< @copydoc q::ez_gradrho
    RowXpr      e_divm()           { return row(q::e_divm    ); } ///< @copydoc q::e_divm
    RowXpr      e_deltarho()       { return row(q::e_deltarho); } ///< @copydoc q::e_deltarho
    RowXpr      rhou()             { return row(q::rhou      ); } ///< @copydoc q::rhou
    RowXpr      rhov()             { return row(q::rhov      ); } ///< @copydoc q::rhov
    RowXpr      rhow()             { return row(q::rhow      ); } ///< @copydoc q::rhow
    RowXpr      rhoE()             { return row(q::rhoE      ); } ///< @copydoc q::rhoE
    RowXpr      rhouu()            { return row(q::rhouu     ); } ///< @copydoc q::rhouu
    RowXpr      rhovv()            { return row(q::rhovv     ); } ///< @copydoc q::rhovv
    RowXpr      rhoww()            { return row(q::rhoww     ); } ///< @copydoc q::rhoww
    RowXpr      rhoEE()            { return row(q::rhoEE     ); } ///< @copydoc q::rhoEE

    ConstRowXpr rho()        const { return row(q::rho       ); } ///< @copydoc q::rho
    ConstRowXpr p()          const { return row(q::p         ); } ///< @copydoc q::p
    ConstRowXpr p2()         const { return row(q::p2        ); } ///< @copydoc q::p2
    ConstRowXpr T()          const { return row(q::T         ); } ///< @copydoc q::T
    ConstRowXpr a()          const { return row(q::a         ); } ///< @copydoc q::a
    ConstRowXpr u()          const { return row(q::u         ); } ///< @copydoc q::u
    ConstRowXpr v()          const { return row(q::v         ); } ///< @copydoc q::v
    ConstRowXpr w()          const { return row(q::w         ); } ///< @copydoc q::w
    ConstRowXpr u2()         const { return row(q::u2        ); } ///< @copydoc q::u2
    ConstRowXpr uu()         const { return row(q::uu        ); } ///< @copydoc q::uu
    ConstRowXpr uv()         const { return row(q::uv        ); } ///< @copydoc q::uv
    ConstRowXpr uw()         const { return row(q::uw        ); } ///< @copydoc q::uw
    ConstRowXpr vv()         const { return row(q::vv        ); } ///< @copydoc q::vv
    ConstRowXpr vw()         const { return row(q::vw        ); } ///< @copydoc q::vw
    ConstRowXpr ww()         const { return row(q::ww        ); } ///< @copydoc q::ww
    ConstRowXpr nu()         const { return row(q::nu        ); } ///< @copydoc q::nu
    ConstRowXpr nu_u()       const { return row(q::nu_u      ); } ///< @copydoc q::nu_u
    ConstRowXpr nu_v()       const { return row(q::nu_v      ); } ///< @copydoc q::nu_v
    ConstRowXpr nu_w()       const { return row(q::nu_w      ); } ///< @copydoc q::nu_w
    ConstRowXpr nu_u2()      const { return row(q::nu_u2     ); } ///< @copydoc q::nu_u2
    ConstRowXpr nu_uu()      const { return row(q::nu_uu     ); } ///< @copydoc q::nu_uu
    ConstRowXpr nu_uv()      const { return row(q::nu_uv     ); } ///< @copydoc q::nu_uv
    ConstRowXpr nu_uw()      const { return row(q::nu_uw     ); } ///< @copydoc q::nu_uw
    ConstRowXpr nu_vv()      const { return row(q::nu_vv     ); } ///< @copydoc q::nu_vv
    ConstRowXpr nu_vw()      const { return row(q::nu_vw     ); } ///< @copydoc q::nu_vw
    ConstRowXpr nu_ww()      const { return row(q::nu_ww     ); } ///< @copydoc q::nu_ww
    ConstRowXpr ex_gradrho() const { return row(q::ex_gradrho); } ///< @copydoc q::ex_gradrho
    ConstRowXpr ey_gradrho() const { return row(q::ey_gradrho); } ///< @copydoc q::ey_gradrho
    ConstRowXpr ez_gradrho() const { return row(q::ez_gradrho); } ///< @copydoc q::ez_gradrho
    ConstRowXpr e_divm()     const { return row(q::e_divm    ); } ///< @copydoc q::e_divm
    ConstRowXpr e_deltarho() const { return row(q::e_deltarho); } ///< @copydoc q::e_deltarho
    ConstRowXpr rhou()       const { return row(q::rhou      ); } ///< @copydoc q::rhou
    ConstRowXpr rhov()       const { return row(q::rhov      ); } ///< @copydoc q::rhov
    ConstRowXpr rhow()       const { return row(q::rhow      ); } ///< @copydoc q::rhow
    ConstRowXpr rhoE()       const { return row(q::rhoE      ); } ///< @copydoc q::rhoE
    ConstRowXpr rhouu()      const { return row(q::rhouu     ); } ///< @copydoc q::rhouu
    ConstRowXpr rhovv()      const { return row(q::rhovv     ); } ///< @copydoc q::rhovv
    ConstRowXpr rhoww()      const { return row(q::rhoww     ); } ///< @copydoc q::rhoww
    ConstRowXpr rhoEE()      const { return row(q::rhoEE     ); } ///< @copydoc q::rhoEE

    /** @} */

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_REFERENCES_HPP */
