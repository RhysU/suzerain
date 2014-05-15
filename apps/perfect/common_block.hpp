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

#ifndef SUZERAIN_PERFECT_COMMON_BLOCK_HPP
#define SUZERAIN_PERFECT_COMMON_BLOCK_HPP

/** @file
 * Declarations shared between linearly-implicit and explicit operators.
 */

#include <suzerain/common.hpp>
#include <suzerain/rholut_imexop.h>
#include <suzerain/treatment_constraint.hpp>

#include "linearize_type.hpp"
#include "slowgrowth_type.hpp"

namespace suzerain {

namespace perfect {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application and implicit constraint application which either are required
 * for linear operator application or for statistics sampling purposes.
 */
class operator_common_block
    : public virtual constraint::treatment::inputs
    , public virtual constraint::treatment::outputs
    , public virtual boost::noncopyable
{
public:

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    operator_common_block();

    /**
     * Determines the extent of the implicit treatment
     * by the paired linear and nonlinear operators.
     */
    linearize::type linearization;

    /**
     * Determines the general type of the slow growth treatment by the paired
     * linear and nonlinear operators.  This is insufficient information to
     * nail down the exact treatment in use as it may depend on other settings,
     * e.g. \ref definition_largo in cases when the Largo library is employed.
     */
    slowgrowth::type slow_treatment;

    /**
     * Helper consistently zeroing \c means, \c implicits, and \c refs.
     *
     * Case deliberately differs from the analogous Eigen signature so both
     * consistency with Suzerain naming and so that it is visually distinct.
     */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        means    .setZero(Ny, 20);  // Magic number from \ref means columns
        implicits.setZero(Ny, 23);  // Magic number from \ref implicits columns
        refs     .setZero(39, Ny);  // Magic number from \ref refs rows
    }

    /**
     * The mean quantities, stored as collocation point values in \c means,
     * are updated by the nonlinear operator on each Runge--Kutta substep.
     *
     * The \e nonlinear operator computes the instantaneous spatial (x, z) mean
     * profiles.  The \e linear operator uses this information to compute the
     * implicit \f$f\cdot{}u\f$ and \f$\mathscr{C}_{\rho{}u}\cdot{}u\f$ terms in
     * the total energy equation.  Additional profiles are gathered for slow
     * growth forcing purposes.
     *
     * Each mean quantity is a single column within \c means.  This facilitates
     * operations across the entire wall-normal profile in a stride one
     * fashion.
     *
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    ArrayXXr means;

    ArrayXXr::ColXpr      u()           { return means.col( 0); }
    ArrayXXr::ColXpr      v()           { return means.col( 1); }
    ArrayXXr::ColXpr      w()           { return means.col( 2); }
    ArrayXXr::ColXpr      uu()          { return means.col( 3); }
    ArrayXXr::ColXpr      uv()          { return means.col( 4); }
    ArrayXXr::ColXpr      uw()          { return means.col( 5); }
    ArrayXXr::ColXpr      vv()          { return means.col( 6); }
    ArrayXXr::ColXpr      vw()          { return means.col( 7); }
    ArrayXXr::ColXpr      ww()          { return means.col( 8); }
    ArrayXXr::ColXpr      rho()         { return means.col( 9); }
    ArrayXXr::ColXpr      rhou()        { return means.col(10); }
    ArrayXXr::ColXpr      rhov()        { return means.col(11); }
    ArrayXXr::ColXpr      rhow()        { return means.col(12); }
    ArrayXXr::ColXpr      rhoE()        { return means.col(13); }
    ArrayXXr::ColXpr      rhouu()       { return means.col(14); }
    ArrayXXr::ColXpr      rhovv()       { return means.col(15); }
    ArrayXXr::ColXpr      rhoww()       { return means.col(16); }
    ArrayXXr::ColXpr      rhoEE()       { return means.col(17); }
    ArrayXXr::ColXpr      p()           { return means.col(18); }
    ArrayXXr::ColXpr      p2()          { return means.col(19); }

    ArrayXXr::ConstColXpr u()     const { return means.col( 0); }
    ArrayXXr::ConstColXpr v()     const { return means.col( 1); }
    ArrayXXr::ConstColXpr w()     const { return means.col( 2); }
    ArrayXXr::ConstColXpr uu()    const { return means.col( 3); }
    ArrayXXr::ConstColXpr uv()    const { return means.col( 4); }
    ArrayXXr::ConstColXpr uw()    const { return means.col( 5); }
    ArrayXXr::ConstColXpr vv()    const { return means.col( 6); }
    ArrayXXr::ConstColXpr vw()    const { return means.col( 7); }
    ArrayXXr::ConstColXpr ww()    const { return means.col( 8); }
    ArrayXXr::ConstColXpr rho()   const { return means.col( 9); }
    ArrayXXr::ConstColXpr rhou()  const { return means.col(10); }
    ArrayXXr::ConstColXpr rhov()  const { return means.col(11); }
    ArrayXXr::ConstColXpr rhow()  const { return means.col(12); }
    ArrayXXr::ConstColXpr rhoE()  const { return means.col(13); }
    ArrayXXr::ConstColXpr rhouu() const { return means.col(14); }
    ArrayXXr::ConstColXpr rhovv() const { return means.col(15); }
    ArrayXXr::ConstColXpr rhoww() const { return means.col(16); }
    ArrayXXr::ConstColXpr rhoEE() const { return means.col(17); }
    ArrayXXr::ConstColXpr p()     const { return means.col(18); }
    ArrayXXr::ConstColXpr p2()    const { return means.col(19); }

    /** @} */

    /**
     * The implicit quantities, stored as collocation point values in \ref
     * implicits and updated by \ref constraint::treatment usage.
     *
     * These "Time-step-specific temporal means" are time averages
     * taken across a single time step of quantities which vary on
     * each substep.  Each quantity is a single column within \c means.
     * This facilitates operations across the entire wall-normal profile
     * in a stride one fashion.
     *
     * @see writeups/treatment_channel.tex and writeups/derivation.tex for
     * details on the mean quantities.
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    ArrayXXr implicits;

    ArrayXXr::ColXpr      SrhoE()              { return implicits.col( 0); }
    ArrayXXr::ColXpr      Srhou()              { return implicits.col( 1); }
    ArrayXXr::ColXpr      Srhov()              { return implicits.col( 2); }
    ArrayXXr::ColXpr      Srhow()              { return implicits.col( 3); }
    ArrayXXr::ColXpr      Srho()               { return implicits.col( 4); }
    ArrayXXr::ColXpr      Srhou_dot_u()        { return implicits.col( 5); }
    ArrayXXr::ColXpr      fx()                 { return implicits.col( 6); }
    ArrayXXr::ColXpr      fy()                 { return implicits.col( 7); }
    ArrayXXr::ColXpr      fz()                 { return implicits.col( 8); }
    ArrayXXr::ColXpr      f_dot_u()            { return implicits.col( 9); }
    ArrayXXr::ColXpr      qb()                 { return implicits.col(10); }
    ArrayXXr::ColXpr      CrhoE()              { return implicits.col(11); }
    ArrayXXr::ColXpr      C2rhoE()             { return implicits.col(12); }
    ArrayXXr::ColXpr      Crhou()              { return implicits.col(13); }
    ArrayXXr::ColXpr      C2rhou()             { return implicits.col(14); }
    ArrayXXr::ColXpr      Crhov()              { return implicits.col(15); }
    ArrayXXr::ColXpr      C2rhov()             { return implicits.col(16); }
    ArrayXXr::ColXpr      Crhow()              { return implicits.col(17); }
    ArrayXXr::ColXpr      C2rhow()             { return implicits.col(18); }
    ArrayXXr::ColXpr      Crho()               { return implicits.col(19); }
    ArrayXXr::ColXpr      C2rho()              { return implicits.col(20); }
    ArrayXXr::ColXpr      Crhou_dot_u()        { return implicits.col(21); }
    ArrayXXr::ColXpr      C2rhou_dot_u()       { return implicits.col(22); }

    ArrayXXr::ConstColXpr SrhoE()        const { return implicits.col( 0); }
    ArrayXXr::ConstColXpr Srhou()        const { return implicits.col( 1); }
    ArrayXXr::ConstColXpr Srhov()        const { return implicits.col( 2); }
    ArrayXXr::ConstColXpr Srhow()        const { return implicits.col( 3); }
    ArrayXXr::ConstColXpr Srho()         const { return implicits.col( 4); }
    ArrayXXr::ConstColXpr Srhou_dot_u()  const { return implicits.col( 5); }
    ArrayXXr::ConstColXpr fx()           const { return implicits.col( 6); }
    ArrayXXr::ConstColXpr fy()           const { return implicits.col( 7); }
    ArrayXXr::ConstColXpr fz()           const { return implicits.col( 8); }
    ArrayXXr::ConstColXpr f_dot_u()      const { return implicits.col( 9); }
    ArrayXXr::ConstColXpr qb()           const { return implicits.col(10); }
    ArrayXXr::ConstColXpr CrhoE()        const { return implicits.col(11); }
    ArrayXXr::ConstColXpr C2rhoE()       const { return implicits.col(12); }
    ArrayXXr::ConstColXpr Crhou()        const { return implicits.col(13); }
    ArrayXXr::ConstColXpr C2rhou()       const { return implicits.col(14); }
    ArrayXXr::ConstColXpr Crhov()        const { return implicits.col(15); }
    ArrayXXr::ConstColXpr C2rhov()       const { return implicits.col(16); }
    ArrayXXr::ConstColXpr Crhow()        const { return implicits.col(17); }
    ArrayXXr::ConstColXpr C2rhow()       const { return implicits.col(18); }
    ArrayXXr::ConstColXpr Crho()         const { return implicits.col(19); }
    ArrayXXr::ConstColXpr C2rho()        const { return implicits.col(20); }
    ArrayXXr::ConstColXpr Crhou_dot_u()  const { return implicits.col(21); }
    ArrayXXr::ConstColXpr C2rhou_dot_u() const { return implicits.col(22); }

    /** @} */

    /**
     *
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

    /** Column-major storage housing all mean quantities (one per row). */
    ArrayXXr refs;

    ArrayXXr::RowXpr      ref_rho()              { return refs.row( 0); }
    ArrayXXr::RowXpr      ref_p()                { return refs.row( 1); }
    ArrayXXr::RowXpr      ref_p2()               { return refs.row( 2); }
    ArrayXXr::RowXpr      ref_T()                { return refs.row( 3); }
    ArrayXXr::RowXpr      ref_a()                { return refs.row( 4); }
    ArrayXXr::RowXpr      ref_ux()               { return refs.row( 5); }
    ArrayXXr::RowXpr      ref_uy()               { return refs.row( 6); }
    ArrayXXr::RowXpr      ref_uz()               { return refs.row( 7); }
    ArrayXXr::RowXpr      ref_u2()               { return refs.row( 8); }
    ArrayXXr::RowXpr      ref_uxux()             { return refs.row( 9); }
    ArrayXXr::RowXpr      ref_uxuy()             { return refs.row(10); }
    ArrayXXr::RowXpr      ref_uxuz()             { return refs.row(11); }
    ArrayXXr::RowXpr      ref_uyuy()             { return refs.row(12); }
    ArrayXXr::RowXpr      ref_uyuz()             { return refs.row(13); }
    ArrayXXr::RowXpr      ref_uzuz()             { return refs.row(14); }
    ArrayXXr::RowXpr      ref_nu()               { return refs.row(15); }
    ArrayXXr::RowXpr      ref_nuux()             { return refs.row(16); }
    ArrayXXr::RowXpr      ref_nuuy()             { return refs.row(17); }
    ArrayXXr::RowXpr      ref_nuuz()             { return refs.row(18); }
    ArrayXXr::RowXpr      ref_nuu2()             { return refs.row(19); }
    ArrayXXr::RowXpr      ref_nuuxux()           { return refs.row(20); }
    ArrayXXr::RowXpr      ref_nuuxuy()           { return refs.row(21); }
    ArrayXXr::RowXpr      ref_nuuxuz()           { return refs.row(22); }
    ArrayXXr::RowXpr      ref_nuuyuy()           { return refs.row(23); }
    ArrayXXr::RowXpr      ref_nuuyuz()           { return refs.row(24); }
    ArrayXXr::RowXpr      ref_nuuzuz()           { return refs.row(25); }
    ArrayXXr::RowXpr      ref_ex_gradrho()       { return refs.row(26); }
    ArrayXXr::RowXpr      ref_ey_gradrho()       { return refs.row(27); }
    ArrayXXr::RowXpr      ref_ez_gradrho()       { return refs.row(28); }
    ArrayXXr::RowXpr      ref_e_divm()           { return refs.row(29); }
    ArrayXXr::RowXpr      ref_e_deltarho()       { return refs.row(30); }
    ArrayXXr::RowXpr      ref_rhoux()            { return refs.row(31); }
    ArrayXXr::RowXpr      ref_rhouy()            { return refs.row(32); }
    ArrayXXr::RowXpr      ref_rhouz()            { return refs.row(33); }
    ArrayXXr::RowXpr      ref_rhoE()             { return refs.row(34); }
    ArrayXXr::RowXpr      ref_rhouxux()          { return refs.row(35); }
    ArrayXXr::RowXpr      ref_rhouyuy()          { return refs.row(36); }
    ArrayXXr::RowXpr      ref_rhouzuz()          { return refs.row(37); }
    ArrayXXr::RowXpr      ref_rhoEE()            { return refs.row(38); }

    ArrayXXr::ConstRowXpr ref_rho()        const { return refs.row( 0); }
    ArrayXXr::ConstRowXpr ref_p()          const { return refs.row( 1); }
    ArrayXXr::ConstRowXpr ref_p2()         const { return refs.row( 2); }
    ArrayXXr::ConstRowXpr ref_T()          const { return refs.row( 3); }
    ArrayXXr::ConstRowXpr ref_a()          const { return refs.row( 4); }
    ArrayXXr::ConstRowXpr ref_ux()         const { return refs.row( 5); }
    ArrayXXr::ConstRowXpr ref_uy()         const { return refs.row( 6); }
    ArrayXXr::ConstRowXpr ref_uz()         const { return refs.row( 7); }
    ArrayXXr::ConstRowXpr ref_u2()         const { return refs.row( 8); }
    ArrayXXr::ConstRowXpr ref_uxux()       const { return refs.row( 9); }
    ArrayXXr::ConstRowXpr ref_uxuy()       const { return refs.row(10); }
    ArrayXXr::ConstRowXpr ref_uxuz()       const { return refs.row(11); }
    ArrayXXr::ConstRowXpr ref_uyuy()       const { return refs.row(12); }
    ArrayXXr::ConstRowXpr ref_uyuz()       const { return refs.row(13); }
    ArrayXXr::ConstRowXpr ref_uzuz()       const { return refs.row(14); }
    ArrayXXr::ConstRowXpr ref_nu()         const { return refs.row(15); }
    ArrayXXr::ConstRowXpr ref_nuux()       const { return refs.row(16); }
    ArrayXXr::ConstRowXpr ref_nuuy()       const { return refs.row(17); }
    ArrayXXr::ConstRowXpr ref_nuuz()       const { return refs.row(18); }
    ArrayXXr::ConstRowXpr ref_nuu2()       const { return refs.row(19); }
    ArrayXXr::ConstRowXpr ref_nuuxux()     const { return refs.row(20); }
    ArrayXXr::ConstRowXpr ref_nuuxuy()     const { return refs.row(21); }
    ArrayXXr::ConstRowXpr ref_nuuxuz()     const { return refs.row(22); }
    ArrayXXr::ConstRowXpr ref_nuuyuy()     const { return refs.row(23); }
    ArrayXXr::ConstRowXpr ref_nuuyuz()     const { return refs.row(24); }
    ArrayXXr::ConstRowXpr ref_nuuzuz()     const { return refs.row(25); }
    ArrayXXr::ConstRowXpr ref_ex_gradrho() const { return refs.row(26); }
    ArrayXXr::ConstRowXpr ref_ey_gradrho() const { return refs.row(27); }
    ArrayXXr::ConstRowXpr ref_ez_gradrho() const { return refs.row(28); }
    ArrayXXr::ConstRowXpr ref_e_divm()     const { return refs.row(29); }
    ArrayXXr::ConstRowXpr ref_e_deltarho() const { return refs.row(30); }
    ArrayXXr::ConstRowXpr ref_rhoux()      const { return refs.row(31); }
    ArrayXXr::ConstRowXpr ref_rhouy()      const { return refs.row(32); }
    ArrayXXr::ConstRowXpr ref_rhouz()      const { return refs.row(33); }
    ArrayXXr::ConstRowXpr ref_rhoE()       const { return refs.row(34); }
    ArrayXXr::ConstRowXpr ref_rhouxux()    const { return refs.row(35); }
    ArrayXXr::ConstRowXpr ref_rhouyuy()    const { return refs.row(36); }
    ArrayXXr::ConstRowXpr ref_rhouzuz()    const { return refs.row(37); }
    ArrayXXr::ConstRowXpr ref_rhoEE()      const { return refs.row(38); }

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

    /** @} */

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_COMMON_BLOCK_HPP */
