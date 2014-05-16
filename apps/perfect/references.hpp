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
 * Storage for holding quantities computed during nonlinear operator application
 * which either are required for hybrid linear operator application or for
 * statistics sampling purposes.
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

    /**
     * The reference quantities stored in \c refs are as follows:
     * \li \c ref_rho        Reference \f$\rho                  \f$
     * \li \c ref_p          Reference \f$p                     \f$
     * \li \c ref_T          Reference \f$T                     \f$
     * \li \c ref_a          Reference \f$a = \sqrt{T}          \f$
     * \li \c ref_ux         Reference \f$C^{u_x}               \f$
     * \li \c ref_uy         Reference \f$C^{u_y}               \f$
     * \li \c ref_uz         Reference \f$C^{u_z}               \f$
     * \li \c ref_uxux       Reference \f$C^{u_x u_x}           \f$
     * \li \c ref_uxuy       Reference \f$C^{u_x u_y}           \f$
     * \li \c ref_uxuz       Reference \f$C^{u_x u_z}           \f$
     * \li \c ref_uyuy       Reference \f$C^{u_y u_y}           \f$
     * \li \c ref_uyuz       Reference \f$C^{u_y u_z}           \f$
     * \li \c ref_uzuz       Reference \f$C^{u_z u_z}           \f$
     * \li \c ref_nu         Reference \f$C^{\nu}               \f$
     * \li \c ref_nuux       Reference \f$C^{\nu u_x}           \f$
     * \li \c ref_nuuy       Reference \f$C^{\nu u_y}           \f$
     * \li \c ref_nuuz       Reference \f$C^{\nu u_z}           \f$
     * \li \c ref_nuu2       Reference \f$C^{\nu u^2}           \f$
     * \li \c ref_nuuxux     Reference \f$C^{\nu u_x u_x}       \f$
     * \li \c ref_nuuxuy     Reference \f$C^{\nu u_x u_y}       \f$
     * \li \c ref_nuuxuz     Reference \f$C^{\nu u_x u_z}       \f$
     * \li \c ref_nuuyuy     Reference \f$C^{\nu u_y u_y}       \f$
     * \li \c ref_nuuyuz     Reference \f$C^{\nu u_y u_z}       \f$
     * \li \c ref_nuuzuz     Reference \f$C^{\nu u_z u_z}       \f$
     * \li \c ref_ex_gradrho Reference \f$C^{e_x}_{\nabla\rho}  \f$
     * \li \c ref_ey_gradrho Reference \f$C^{e_y}_{\nabla\rho}  \f$
     * \li \c ref_ez_gradrho Reference \f$C^{e_z}_{\nabla\rho}  \f$
     * \li \c ref_e_divm     Reference \f$C^{e}_{\nabla\cdot{}m}\f$
     * \li \c ref_e_deltarho Reference \f$C^{e}_{\Delta\rho}    \f$
     *
     * The following quantities are also stored, though they are not
     * used for linearization purposes:
     * \li \c ref_p2         Quantity \f$p^2         \f$
     * \li \c ref_rhoux      Quantity \f$\rho u_x    \f$
     * \li \c ref_rhouy      Quantity \f$\rho u_y    \f$
     * \li \c ref_rhouz      Quantity \f$\rho u_z    \f$
     * \li \c ref_rhoE       Quantity \f$\rho E      \f$
     * \li \c ref_rhouxux    Quantity \f$\rho u_x u_x\f$
     * \li \c ref_rhouyuy    Quantity \f$\rho u_y u_y\f$
     * \li \c ref_rhouzuz    Quantity \f$\rho u_z u_z\f$
     * \li \c ref_rhoEE      Quantity \f$\rho E   E  \f$
     * Here, \f$E\f$ denotes \f$e/\rho\f$.  Mean state at collocation points
     * is known on rank zero in wave space, but is additionally tracked
     * here for numerical uniformity with quantities like \c ref_rhoEE.
     *
     * Each reference quantity is a single row within \c refs.  This
     * facilitates a stride one operation loading or writing all reference
     * quantities for a single wall-normal location.
     *
     * @see rholut_imexop.h for details on linearization and the associated
     * reference quantities.
     *
     * @{
     */

    RowXpr      ref_rho()              { return row( 0); }
    RowXpr      ref_p()                { return row( 1); }
    RowXpr      ref_p2()               { return row( 2); }
    RowXpr      ref_T()                { return row( 3); }
    RowXpr      ref_a()                { return row( 4); }
    RowXpr      ref_ux()               { return row( 5); }
    RowXpr      ref_uy()               { return row( 6); }
    RowXpr      ref_uz()               { return row( 7); }
    RowXpr      ref_u2()               { return row( 8); }
    RowXpr      ref_uxux()             { return row( 9); }
    RowXpr      ref_uxuy()             { return row(10); }
    RowXpr      ref_uxuz()             { return row(11); }
    RowXpr      ref_uyuy()             { return row(12); }
    RowXpr      ref_uyuz()             { return row(13); }
    RowXpr      ref_uzuz()             { return row(14); }
    RowXpr      ref_nu()               { return row(15); }
    RowXpr      ref_nuux()             { return row(16); }
    RowXpr      ref_nuuy()             { return row(17); }
    RowXpr      ref_nuuz()             { return row(18); }
    RowXpr      ref_nuu2()             { return row(19); }
    RowXpr      ref_nuuxux()           { return row(20); }
    RowXpr      ref_nuuxuy()           { return row(21); }
    RowXpr      ref_nuuxuz()           { return row(22); }
    RowXpr      ref_nuuyuy()           { return row(23); }
    RowXpr      ref_nuuyuz()           { return row(24); }
    RowXpr      ref_nuuzuz()           { return row(25); }
    RowXpr      ref_ex_gradrho()       { return row(26); }
    RowXpr      ref_ey_gradrho()       { return row(27); }
    RowXpr      ref_ez_gradrho()       { return row(28); }
    RowXpr      ref_e_divm()           { return row(29); }
    RowXpr      ref_e_deltarho()       { return row(30); }
    RowXpr      ref_rhoux()            { return row(31); }
    RowXpr      ref_rhouy()            { return row(32); }
    RowXpr      ref_rhouz()            { return row(33); }
    RowXpr      ref_rhoE()             { return row(34); }
    RowXpr      ref_rhouxux()          { return row(35); }
    RowXpr      ref_rhouyuy()          { return row(36); }
    RowXpr      ref_rhouzuz()          { return row(37); }
    RowXpr      ref_rhoEE()            { return row(38); }

    ConstRowXpr ref_rho()        const { return row( 0); }
    ConstRowXpr ref_p()          const { return row( 1); }
    ConstRowXpr ref_p2()         const { return row( 2); }
    ConstRowXpr ref_T()          const { return row( 3); }
    ConstRowXpr ref_a()          const { return row( 4); }
    ConstRowXpr ref_ux()         const { return row( 5); }
    ConstRowXpr ref_uy()         const { return row( 6); }
    ConstRowXpr ref_uz()         const { return row( 7); }
    ConstRowXpr ref_u2()         const { return row( 8); }
    ConstRowXpr ref_uxux()       const { return row( 9); }
    ConstRowXpr ref_uxuy()       const { return row(10); }
    ConstRowXpr ref_uxuz()       const { return row(11); }
    ConstRowXpr ref_uyuy()       const { return row(12); }
    ConstRowXpr ref_uyuz()       const { return row(13); }
    ConstRowXpr ref_uzuz()       const { return row(14); }
    ConstRowXpr ref_nu()         const { return row(15); }
    ConstRowXpr ref_nuux()       const { return row(16); }
    ConstRowXpr ref_nuuy()       const { return row(17); }
    ConstRowXpr ref_nuuz()       const { return row(18); }
    ConstRowXpr ref_nuu2()       const { return row(19); }
    ConstRowXpr ref_nuuxux()     const { return row(20); }
    ConstRowXpr ref_nuuxuy()     const { return row(21); }
    ConstRowXpr ref_nuuxuz()     const { return row(22); }
    ConstRowXpr ref_nuuyuy()     const { return row(23); }
    ConstRowXpr ref_nuuyuz()     const { return row(24); }
    ConstRowXpr ref_nuuzuz()     const { return row(25); }
    ConstRowXpr ref_ex_gradrho() const { return row(26); }
    ConstRowXpr ref_ey_gradrho() const { return row(27); }
    ConstRowXpr ref_ez_gradrho() const { return row(28); }
    ConstRowXpr ref_e_divm()     const { return row(29); }
    ConstRowXpr ref_e_deltarho() const { return row(30); }
    ConstRowXpr ref_rhoux()      const { return row(31); }
    ConstRowXpr ref_rhouy()      const { return row(32); }
    ConstRowXpr ref_rhouz()      const { return row(33); }
    ConstRowXpr ref_rhoE()       const { return row(34); }
    ConstRowXpr ref_rhouxux()    const { return row(35); }
    ConstRowXpr ref_rhouyuy()    const { return row(36); }
    ConstRowXpr ref_rhouzuz()    const { return row(37); }
    ConstRowXpr ref_rhoEE()      const { return row(38); }

    /** @} */

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_REFERENCES_HPP */
