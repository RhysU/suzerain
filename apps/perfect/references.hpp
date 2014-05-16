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
    RowXpr ux()         { return row( 5); } ///< Reference \f$C^{u_x}               \f$
    RowXpr uy()         { return row( 6); } ///< Reference \f$C^{u_y}               \f$
    RowXpr uz()         { return row( 7); } ///< Reference \f$C^{u_z}               \f$
    RowXpr u2()         { return row( 8); } ///< Quantity  \f$\vec{u}^2             \f$
    RowXpr uxux()       { return row( 9); } ///< Reference \f$C^{u_x u_x}           \f$
    RowXpr uxuy()       { return row(10); } ///< Reference \f$C^{u_x u_y}           \f$
    RowXpr uxuz()       { return row(11); } ///< Reference \f$C^{u_x u_z}           \f$
    RowXpr uyuy()       { return row(12); } ///< Reference \f$C^{u_y u_y}           \f$
    RowXpr uyuz()       { return row(13); } ///< Reference \f$C^{u_y u_z}           \f$
    RowXpr uzuz()       { return row(14); } ///< Reference \f$C^{u_z u_z}           \f$
    RowXpr nu()         { return row(15); } ///< Reference \f$C^{\nu}               \f$
    RowXpr nu_ux()      { return row(16); } ///< Reference \f$C^{\nu u_x}           \f$
    RowXpr nu_uy()      { return row(17); } ///< Reference \f$C^{\nu u_y}           \f$
    RowXpr nu_uz()      { return row(18); } ///< Reference \f$C^{\nu u_z}           \f$
    RowXpr nu_u2()      { return row(19); } ///< Reference \f$C^{\nu u^2}           \f$
    RowXpr nu_uxux()    { return row(20); } ///< Reference \f$C^{\nu u_x u_x}       \f$
    RowXpr nu_uxuy()    { return row(21); } ///< Reference \f$C^{\nu u_x u_y}       \f$
    RowXpr nu_uxuz()    { return row(22); } ///< Reference \f$C^{\nu u_x u_z}       \f$
    RowXpr nu_uyuy()    { return row(23); } ///< Reference \f$C^{\nu u_y u_y}       \f$
    RowXpr nu_uyuz()    { return row(24); } ///< Reference \f$C^{\nu u_y u_z}       \f$
    RowXpr nu_uzuz()    { return row(25); } ///< Reference \f$C^{\nu u_z u_z}       \f$
    RowXpr ex_gradrho() { return row(26); } ///< Reference \f$C^{e_x}_{\nabla\rho}  \f$
    RowXpr ey_gradrho() { return row(27); } ///< Reference \f$C^{e_y}_{\nabla\rho}  \f$
    RowXpr ez_gradrho() { return row(28); } ///< Reference \f$C^{e_z}_{\nabla\rho}  \f$
    RowXpr e_divm()     { return row(29); } ///< Reference \f$C^{e}_{\nabla\cdot{}m}\f$
    RowXpr e_deltarho() { return row(30); } ///< Reference \f$C^{e}_{\Delta\rho}    \f$
    RowXpr rhoux()      { return row(31); } ///< Quantity  \f$\rho u_x              \f$
    RowXpr rhouy()      { return row(32); } ///< Quantity  \f$\rho u_y              \f$
    RowXpr rhouz()      { return row(33); } ///< Quantity  \f$\rho u_z              \f$
    RowXpr rhoE()       { return row(34); } ///< Quantity  \f$\rho E                \f$
    RowXpr rhouxux()    { return row(35); } ///< Quantity  \f$\rho u_x u_x          \f$
    RowXpr rhouyuy()    { return row(36); } ///< Quantity  \f$\rho u_y u_y          \f$
    RowXpr rhouzuz()    { return row(37); } ///< Quantity  \f$\rho u_z u_z          \f$
    RowXpr rhoEE()      { return row(38); } ///< Quantity  \f$\rho E   E            \f$

    ConstRowXpr rho()        const { return row( 0); } ///< @copydoc rho()
    ConstRowXpr p()          const { return row( 1); } ///< @copydoc p()
    ConstRowXpr p2()         const { return row( 2); } ///< @copydoc p2()
    ConstRowXpr T()          const { return row( 3); } ///< @copydoc T()
    ConstRowXpr a()          const { return row( 4); } ///< @copydoc a()
    ConstRowXpr ux()         const { return row( 5); } ///< @copydoc ux()
    ConstRowXpr uy()         const { return row( 6); } ///< @copydoc uy()
    ConstRowXpr uz()         const { return row( 7); } ///< @copydoc uz()
    ConstRowXpr u2()         const { return row( 8); } ///< @copydoc u2()
    ConstRowXpr uxux()       const { return row( 9); } ///< @copydoc uxux()
    ConstRowXpr uxuy()       const { return row(10); } ///< @copydoc uxuy()
    ConstRowXpr uxuz()       const { return row(11); } ///< @copydoc uxuz()
    ConstRowXpr uyuy()       const { return row(12); } ///< @copydoc uyuy()
    ConstRowXpr uyuz()       const { return row(13); } ///< @copydoc uyuz()
    ConstRowXpr uzuz()       const { return row(14); } ///< @copydoc uzuz()
    ConstRowXpr nu()         const { return row(15); } ///< @copydoc nu()
    ConstRowXpr nu_ux()      const { return row(16); } ///< @copydoc nu_ux()
    ConstRowXpr nu_uy()      const { return row(17); } ///< @copydoc nu_uy()
    ConstRowXpr nu_uz()      const { return row(18); } ///< @copydoc nu_uz()
    ConstRowXpr nu_u2()      const { return row(19); } ///< @copydoc nu_u2()
    ConstRowXpr nu_uxux()    const { return row(20); } ///< @copydoc nu_uxux()
    ConstRowXpr nu_uxuy()    const { return row(21); } ///< @copydoc nu_uxuy()
    ConstRowXpr nu_uxuz()    const { return row(22); } ///< @copydoc nu_uxuz()
    ConstRowXpr nu_uyuy()    const { return row(23); } ///< @copydoc nu_uyuy()
    ConstRowXpr nu_uyuz()    const { return row(24); } ///< @copydoc nu_uyuz()
    ConstRowXpr nu_uzuz()    const { return row(25); } ///< @copydoc nu_uzuz()
    ConstRowXpr ex_gradrho() const { return row(26); } ///< @copydoc ex_gradrho()
    ConstRowXpr ey_gradrho() const { return row(27); } ///< @copydoc ey_gradrho()
    ConstRowXpr ez_gradrho() const { return row(28); } ///< @copydoc ez_gradrho()
    ConstRowXpr e_divm()     const { return row(29); } ///< @copydoc e_divm()
    ConstRowXpr e_deltarho() const { return row(30); } ///< @copydoc e_deltarho()
    ConstRowXpr rhoux()      const { return row(31); } ///< @copydoc rhoux()
    ConstRowXpr rhouy()      const { return row(32); } ///< @copydoc rhouy()
    ConstRowXpr rhouz()      const { return row(33); } ///< @copydoc rhouz()
    ConstRowXpr rhoE()       const { return row(34); } ///< @copydoc rhoE()
    ConstRowXpr rhouxux()    const { return row(35); } ///< @copydoc rhouxux()
    ConstRowXpr rhouyuy()    const { return row(36); } ///< @copydoc rhouyuy()
    ConstRowXpr rhouzuz()    const { return row(37); } ///< @copydoc rhouzuz()
    ConstRowXpr rhoEE()      const { return row(38); } ///< @copydoc rhoEE()

    /** @} */

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_REFERENCES_HPP */
