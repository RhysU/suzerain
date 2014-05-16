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
 * Each reference quantity is a single row in a column-major matrix.  This
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

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    references(): super(39, 0) {} // Magic number from count of quantities

    typedef super::Scalar      Scalar;      ///< Expose underlying scalar type
    typedef super::RowXpr      RowXpr;      ///< Mutable data for quantities
    typedef super::ConstRowXpr ConstRowXpr; ///< Immutable data for quantities
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

    RowXpr rho()        { return row( 0); } ///< Reference \f$\rho                  \f$
    RowXpr p()          { return row( 1); } ///< Reference \f$p                     \f$
    RowXpr p2()         { return row( 2); } ///< Quantity  \f$p^2                   \f$
    RowXpr T()          { return row( 3); } ///< Reference \f$T                     \f$
    RowXpr a()          { return row( 4); } ///< Reference \f$a = \sqrt{T}          \f$
    RowXpr u()          { return row( 5); } ///< Reference \f$C^{u_x}               \f$
    RowXpr v()          { return row( 6); } ///< Reference \f$C^{u_y}               \f$
    RowXpr w()          { return row( 7); } ///< Reference \f$C^{u_z}               \f$
    RowXpr u2()         { return row( 8); } ///< Quantity  \f$\vec{u}^2             \f$
    RowXpr uu()         { return row( 9); } ///< Reference \f$C^{u_x u_x}           \f$
    RowXpr uv()         { return row(10); } ///< Reference \f$C^{u_x u_y}           \f$
    RowXpr uw()         { return row(11); } ///< Reference \f$C^{u_x u_z}           \f$
    RowXpr vv()         { return row(12); } ///< Reference \f$C^{u_y u_y}           \f$
    RowXpr vw()         { return row(13); } ///< Reference \f$C^{u_y u_z}           \f$
    RowXpr ww()         { return row(14); } ///< Reference \f$C^{u_z u_z}           \f$
    RowXpr nu()         { return row(15); } ///< Reference \f$C^{\nu}               \f$
    RowXpr nu_u()       { return row(16); } ///< Reference \f$C^{\nu u_x}           \f$
    RowXpr nu_v()       { return row(17); } ///< Reference \f$C^{\nu u_y}           \f$
    RowXpr nu_w()       { return row(18); } ///< Reference \f$C^{\nu u_z}           \f$
    RowXpr nu_u2()      { return row(19); } ///< Reference \f$C^{\nu u^2}           \f$
    RowXpr nu_uu()      { return row(20); } ///< Reference \f$C^{\nu u_x u_x}       \f$
    RowXpr nu_uv()      { return row(21); } ///< Reference \f$C^{\nu u_x u_y}       \f$
    RowXpr nu_uw()      { return row(22); } ///< Reference \f$C^{\nu u_x u_z}       \f$
    RowXpr nu_vv()      { return row(23); } ///< Reference \f$C^{\nu u_y u_y}       \f$
    RowXpr nu_vw()      { return row(24); } ///< Reference \f$C^{\nu u_y u_z}       \f$
    RowXpr nu_ww()      { return row(25); } ///< Reference \f$C^{\nu u_z u_z}       \f$
    RowXpr ex_gradrho() { return row(26); } ///< Reference \f$C^{e_x}_{\nabla\rho}  \f$
    RowXpr ey_gradrho() { return row(27); } ///< Reference \f$C^{e_y}_{\nabla\rho}  \f$
    RowXpr ez_gradrho() { return row(28); } ///< Reference \f$C^{e_z}_{\nabla\rho}  \f$
    RowXpr e_divm()     { return row(29); } ///< Reference \f$C^{e}_{\nabla\cdot{}m}\f$
    RowXpr e_deltarho() { return row(30); } ///< Reference \f$C^{e}_{\Delta\rho}    \f$
    RowXpr rhou()       { return row(31); } ///< Quantity  \f$\rho u_x              \f$
    RowXpr rhov()       { return row(32); } ///< Quantity  \f$\rho u_y              \f$
    RowXpr rhow()       { return row(33); } ///< Quantity  \f$\rho u_z              \f$
    RowXpr rhoE()       { return row(34); } ///< Quantity  \f$\rho E                \f$
    RowXpr rhouu()      { return row(35); } ///< Quantity  \f$\rho u_x u_x          \f$
    RowXpr rhovv()      { return row(36); } ///< Quantity  \f$\rho u_y u_y          \f$
    RowXpr rhoww()      { return row(37); } ///< Quantity  \f$\rho u_z u_z          \f$
    RowXpr rhoEE()      { return row(38); } ///< Quantity  \f$\rho E   E            \f$

    ConstRowXpr rho()        const { return row( 0); } ///< @copydoc rho()
    ConstRowXpr p()          const { return row( 1); } ///< @copydoc p()
    ConstRowXpr p2()         const { return row( 2); } ///< @copydoc p2()
    ConstRowXpr T()          const { return row( 3); } ///< @copydoc T()
    ConstRowXpr a()          const { return row( 4); } ///< @copydoc a()
    ConstRowXpr u()          const { return row( 5); } ///< @copydoc u()
    ConstRowXpr v()          const { return row( 6); } ///< @copydoc v()
    ConstRowXpr w()          const { return row( 7); } ///< @copydoc w()
    ConstRowXpr u2()         const { return row( 8); } ///< @copydoc u2()
    ConstRowXpr uu()         const { return row( 9); } ///< @copydoc uu()
    ConstRowXpr uv()         const { return row(10); } ///< @copydoc uv()
    ConstRowXpr uw()         const { return row(11); } ///< @copydoc uw()
    ConstRowXpr vv()         const { return row(12); } ///< @copydoc vv()
    ConstRowXpr vw()         const { return row(13); } ///< @copydoc vw()
    ConstRowXpr ww()         const { return row(14); } ///< @copydoc ww()
    ConstRowXpr nu()         const { return row(15); } ///< @copydoc nu()
    ConstRowXpr nu_u()       const { return row(16); } ///< @copydoc nu_u()
    ConstRowXpr nu_v()       const { return row(17); } ///< @copydoc nu_v()
    ConstRowXpr nu_w()       const { return row(18); } ///< @copydoc nu_w()
    ConstRowXpr nu_u2()      const { return row(19); } ///< @copydoc nu_u2()
    ConstRowXpr nu_uu()      const { return row(20); } ///< @copydoc nu_uu()
    ConstRowXpr nu_uv()      const { return row(21); } ///< @copydoc nu_uv()
    ConstRowXpr nu_uw()      const { return row(22); } ///< @copydoc nu_uw()
    ConstRowXpr nu_vv()      const { return row(23); } ///< @copydoc nu_vv()
    ConstRowXpr nu_vw()      const { return row(24); } ///< @copydoc nu_vw()
    ConstRowXpr nu_ww()      const { return row(25); } ///< @copydoc nu_ww()
    ConstRowXpr ex_gradrho() const { return row(26); } ///< @copydoc ex_gradrho()
    ConstRowXpr ey_gradrho() const { return row(27); } ///< @copydoc ey_gradrho()
    ConstRowXpr ez_gradrho() const { return row(28); } ///< @copydoc ez_gradrho()
    ConstRowXpr e_divm()     const { return row(29); } ///< @copydoc e_divm()
    ConstRowXpr e_deltarho() const { return row(30); } ///< @copydoc e_deltarho()
    ConstRowXpr rhou()       const { return row(31); } ///< @copydoc rhou()
    ConstRowXpr rhov()       const { return row(32); } ///< @copydoc rhov()
    ConstRowXpr rhow()       const { return row(33); } ///< @copydoc rhow()
    ConstRowXpr rhoE()       const { return row(34); } ///< @copydoc rhoE()
    ConstRowXpr rhouu()      const { return row(35); } ///< @copydoc rhouu()
    ConstRowXpr rhovv()      const { return row(36); } ///< @copydoc rhovv()
    ConstRowXpr rhoww()      const { return row(37); } ///< @copydoc rhoww()
    ConstRowXpr rhoEE()      const { return row(38); } ///< @copydoc rhoEE()

    /** @} */

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_REFERENCES_HPP */
