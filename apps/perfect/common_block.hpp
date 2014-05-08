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
{
    /** Type of the contiguous storage housing all mean quantities */
    typedef Array<real_t, Dynamic, 16, ColMajor> means_type;

    /** Type of the contiguous storage housing all implicit quantities */
    typedef Array<real_t, Dynamic, 23, ColMajor> implicits_type;

    /** Type of the contiguous storage housing all reference quantities */
    typedef Array<real_t, 35, Dynamic, ColMajor> refs_type;

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
     * The mean quantities, stored as collocation point values in \c means,
     * are as follows:
     *
     * \li \c u  The \e nonlinear operator computes the instantaneous spatial
     *     (x, z) mean streamwise velocity profile.  The \e linear operator
     *     then uses the information to compute the implicit \f$f\cdot{}u\f$
     *     and \f$\mathscr{C}_{\rho{}u}\cdot{}u\f$ terms in the total energy
     *     equation.
     * \li \c v  Treated identically to \c u.
     * \li \c w  Treated identically to \c u.
     * \li \c uu Gathered by the \e nonlinear operator for miscellaneous use,
     *           for example \ref treatment.
     * \li \c uv Treated identically to \c uu.
     * \li \c uw Treated identically to \c uu.
     * \li \c vv Treated identically to \c uu.
     * \li \c vw Treated identically to \c uu.
     * \li \c ww Treated identically to \c uu.
     *
     * Each mean quantity is a single column within \c means.  This facilitates
     * operations across the entire wall-normal profile in a stride one
     * fashion.
     *
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    means_type means;

    /** Type returned by the non-const mean quantity accessors. */
    typedef means_type::ColXpr mean_type;

    mean_type       u()           { return means.col( 0); }
    mean_type       v()           { return means.col( 1); }
    mean_type       w()           { return means.col( 2); }
    mean_type       uu()          { return means.col( 3); }
    mean_type       uv()          { return means.col( 4); }
    mean_type       uw()          { return means.col( 5); }
    mean_type       vv()          { return means.col( 6); }
    mean_type       vw()          { return means.col( 7); }
    mean_type       ww()          { return means.col( 8); }
    mean_type       rho()         { return means.col( 9); }
    mean_type       rhouu()       { return means.col(10); }
    mean_type       rhovv()       { return means.col(11); }
    mean_type       rhoww()       { return means.col(12); }
    mean_type       rhoEE()       { return means.col(13); }
    mean_type       p()           { return means.col(14); }
    mean_type       p2()          { return means.col(15); }

    /** Type returned by the const mean quantity accessors. */
    typedef means_type::ConstColXpr const_mean_type;

    const_mean_type u()     const { return means.col( 0); }
    const_mean_type v()     const { return means.col( 1); }
    const_mean_type w()     const { return means.col( 2); }
    const_mean_type uu()    const { return means.col( 3); }
    const_mean_type uv()    const { return means.col( 4); }
    const_mean_type uw()    const { return means.col( 5); }
    const_mean_type vv()    const { return means.col( 6); }
    const_mean_type vw()    const { return means.col( 7); }
    const_mean_type ww()    const { return means.col( 8); }
    const_mean_type rho()   const { return means.col( 9); }
    const_mean_type rhouu() const { return means.col(10); }
    const_mean_type rhovv() const { return means.col(11); }
    const_mean_type rhoww() const { return means.col(12); }
    const_mean_type rhoEE() const { return means.col(13); }
    const_mean_type p()     const { return means.col(14); }
    const_mean_type p2()    const { return means.col(15); }

    /** @} */

    /**
     * The implicit quantities, stored as collocation point values in \c
     * implicits, are as follows:
     *
     * \li \c SrhoE The \e nonlinear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{S}_{\rho{}E}\f$ term in
     *     the energy equation.
     * \li \c Srhou The \e nonlinear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{S}_{\rho{}u}\f$ term
     *     in the streamwise (x) momentum equation.
     * \li \c Srhov The \e nonlinear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{S}_{\rho{}v}\f$ term
     *     in the wall-normal (y) momentum equation.
     * \li \c Srhow The \e nonlinear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{S}_{\rho{}w}\f$ term
     *     in the spanwise momentum equation.
     * \li \c Srho The \e nonlinear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{S}_{\rho{}}\f$ term
     *     in the density equation.
     * \li \c Srhou_dot_u The \e nonlinear operator accumulates the
     *     time-step-specific temporal mean of the implicit
     *     \f$\mathscr{S}_{\rho{}u}\cdot{}u\f$ term in the energy equation.
     * \li \c fx The \e linear operator accumulates the time-step-specific
     *     temporal mean streamwise (x) component of the implicit \f$f\f$
     *     term in the momentum equation.
     * \li \c fy The \e linear operator accumulates the time-step-specific
     *     temporal mean wall-normal (y) component of the implicit \f$f\f$
     *     term in the momentum equation.
     * \li \c fz The \e linear operator accumulates the time-step-specific
     *     temporal mean spanwise (z) component of the implicit \f$f\f$
     *     term in the momentum equation.
     * \li \c f_dot_u The \e linear operator accumulates the
     *     time-step-specific temporal mean of the implicit \f$f\cdot{}u\f$
     *     term in the energy equation.
     * \li \c qb The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$q_b\f$ term in the energy equation.
     * \li \c CrhoE The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}E}\f$ term in
     *     the energy equation.
     * \li \c C2rhoE The \e linear operator accumulates the time-step-specific
     *     temporal mean of the quantity \f$\mathscr{C}^2_{\rho{}E}\f$.
     * \li \c Crhou The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}u}\f$ term in
     *     the streamwise (x) momentum equation.
     * \li \c C2rhou The \e linear operator accumulates the time-step-specific
     *     temporal mean of the quantity \f$\mathscr{C}^2_{\rho{}u}\f$.
     * \li \c Crhov The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}v}\f$ term in
     *     the wall-normal (y) momentum equation.
     * \li \c C2rhov The \e linear operator accumulates the time-step-specific
     *     temporal mean of the quantity \f$\mathscr{C}^2_{\rho{}v}\f$.
     * \li \c Crhow The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}w}\f$ term in
     *     the spanwise momentum equation.
     * \li \c C2rhow The \e linear operator accumulates the time-step-specific
     *     temporal mean of the quantity \f$\mathscr{C}^2_{\rho{}w}\f$.
     * \li \c Crho The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}}\f$ term in the
     *     density equation.
     * \li \c C2rho The \e linear operator accumulates the time-step-specific
     *     temporal mean of the quantity \f$\mathscr{C}^2_{\rho{}}\f$.
     * \li \c Crhou_dot_u The \e linear operator accumulates the
     *     time-step-specific temporal mean of the implicit
     *     \f$\mathscr{C}_{\rho{}u}\cdot{}u\f$ term in the energy equation.
     * \li \c C2rhou_dot_u The \e linear operator accumulates the
     *     time-step-specific temporal mean of the quantity
     *     \f$\mathscr{C}^2_{\rho{}u}\cdot{}u\f$.
     *
     * "Time-step-specific temporal means" are time averages taken across a
     * single time step of quantities which vary on each substep.  Each
     * quantity is a single column within \c means.  This facilitates
     * operations across the entire wall-normal profile in a stride one
     * fashion.
     *
     * @see writeups/treatment_channel.tex and writeups/derivation.tex for
     * details on the mean quantities.
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    implicits_type implicits;

    /** Type returned by the non-const implicit quantity accessors. */
    typedef implicits_type::ColXpr implicit_type;

    implicit_type       SrhoE()              { return implicits.col( 0); }
    implicit_type       Srhou()              { return implicits.col( 1); }
    implicit_type       Srhov()              { return implicits.col( 2); }
    implicit_type       Srhow()              { return implicits.col( 3); }
    implicit_type       Srho()               { return implicits.col( 4); }
    implicit_type       Srhou_dot_u()        { return implicits.col( 5); }
    implicit_type       fx()                 { return implicits.col( 6); }
    implicit_type       fy()                 { return implicits.col( 7); }
    implicit_type       fz()                 { return implicits.col( 8); }
    implicit_type       f_dot_u()            { return implicits.col( 9); }
    implicit_type       qb()                 { return implicits.col(10); }
    implicit_type       CrhoE()              { return implicits.col(11); }
    implicit_type       C2rhoE()             { return implicits.col(12); }
    implicit_type       Crhou()              { return implicits.col(13); }
    implicit_type       C2rhou()             { return implicits.col(14); }
    implicit_type       Crhov()              { return implicits.col(15); }
    implicit_type       C2rhov()             { return implicits.col(16); }
    implicit_type       Crhow()              { return implicits.col(17); }
    implicit_type       C2rhow()             { return implicits.col(18); }
    implicit_type       Crho()               { return implicits.col(19); }
    implicit_type       C2rho()              { return implicits.col(20); }
    implicit_type       Crhou_dot_u()        { return implicits.col(21); }
    implicit_type       C2rhou_dot_u()       { return implicits.col(22); }

    /** Type returned by the const implicit quantity accessors. */
    typedef implicits_type::ConstColXpr const_implicit_type;

    const_implicit_type SrhoE()        const { return implicits.col( 0); }
    const_implicit_type Srhou()        const { return implicits.col( 1); }
    const_implicit_type Srhov()        const { return implicits.col( 2); }
    const_implicit_type Srhow()        const { return implicits.col( 3); }
    const_implicit_type Srho()         const { return implicits.col( 4); }
    const_implicit_type Srhou_dot_u()  const { return implicits.col( 5); }
    const_implicit_type fx()           const { return implicits.col( 6); }
    const_implicit_type fy()           const { return implicits.col( 7); }
    const_implicit_type fz()           const { return implicits.col( 8); }
    const_implicit_type f_dot_u()      const { return implicits.col( 9); }
    const_implicit_type qb()           const { return implicits.col(10); }
    const_implicit_type CrhoE()        const { return implicits.col(11); }
    const_implicit_type C2rhoE()       const { return implicits.col(12); }
    const_implicit_type Crhou()        const { return implicits.col(13); }
    const_implicit_type C2rhou()       const { return implicits.col(14); }
    const_implicit_type Crhov()        const { return implicits.col(15); }
    const_implicit_type C2rhov()       const { return implicits.col(16); }
    const_implicit_type Crhow()        const { return implicits.col(17); }
    const_implicit_type C2rhow()       const { return implicits.col(18); }
    const_implicit_type Crho()         const { return implicits.col(19); }
    const_implicit_type C2rho()        const { return implicits.col(20); }
    const_implicit_type Crhou_dot_u()  const { return implicits.col(21); }
    const_implicit_type C2rhou_dot_u() const { return implicits.col(22); }

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
     * \li \c ref_rhouxux    Quantity \f$\rho u_x u_x\f$
     * \li \c ref_rhouyuy    Quantity \f$\rho u_y u_y\f$
     * \li \c ref_rhouzuz    Quantity \f$\rho u_z u_z\f$
     * \li \c ref_rhoEE      Quantity \f$\rho E   E  \f$
     * Here, \f$E\f$ denotes \f$e/\rho\f$.
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
    refs_type refs;

    /** Type returned by the non-const reference quantity accessors. */
    typedef refs_type::RowXpr ref_type;

    ref_type       ref_rho()              { return refs.row( 0); }
    ref_type       ref_p()                { return refs.row( 1); }
    ref_type       ref_p2()               { return refs.row( 2); }
    ref_type       ref_T()                { return refs.row( 3); }
    ref_type       ref_a()                { return refs.row( 4); }
    ref_type       ref_ux()               { return refs.row( 5); }
    ref_type       ref_uy()               { return refs.row( 6); }
    ref_type       ref_uz()               { return refs.row( 7); }
    ref_type       ref_u2()               { return refs.row( 8); }
    ref_type       ref_uxux()             { return refs.row( 9); }
    ref_type       ref_uxuy()             { return refs.row(10); }
    ref_type       ref_uxuz()             { return refs.row(11); }
    ref_type       ref_uyuy()             { return refs.row(12); }
    ref_type       ref_uyuz()             { return refs.row(13); }
    ref_type       ref_uzuz()             { return refs.row(14); }
    ref_type       ref_nu()               { return refs.row(15); }
    ref_type       ref_nuux()             { return refs.row(16); }
    ref_type       ref_nuuy()             { return refs.row(17); }
    ref_type       ref_nuuz()             { return refs.row(18); }
    ref_type       ref_nuu2()             { return refs.row(19); }
    ref_type       ref_nuuxux()           { return refs.row(20); }
    ref_type       ref_nuuxuy()           { return refs.row(21); }
    ref_type       ref_nuuxuz()           { return refs.row(22); }
    ref_type       ref_nuuyuy()           { return refs.row(23); }
    ref_type       ref_nuuyuz()           { return refs.row(24); }
    ref_type       ref_nuuzuz()           { return refs.row(25); }
    ref_type       ref_ex_gradrho()       { return refs.row(26); }
    ref_type       ref_ey_gradrho()       { return refs.row(27); }
    ref_type       ref_ez_gradrho()       { return refs.row(28); }
    ref_type       ref_e_divm()           { return refs.row(29); }
    ref_type       ref_e_deltarho()       { return refs.row(30); }
    ref_type       ref_rhouxux()          { return refs.row(31); }
    ref_type       ref_rhouyuy()          { return refs.row(32); }
    ref_type       ref_rhouzuz()          { return refs.row(33); }
    ref_type       ref_rhoEE()            { return refs.row(34); }

    /** Type returned by the const reference quantity accessors. */
    typedef refs_type::ConstRowXpr const_ref_type;

    const_ref_type ref_rho()        const { return refs.row( 0); }
    const_ref_type ref_p()          const { return refs.row( 1); }
    const_ref_type ref_p2()         const { return refs.row( 2); }
    const_ref_type ref_T()          const { return refs.row( 3); }
    const_ref_type ref_a()          const { return refs.row( 4); }
    const_ref_type ref_ux()         const { return refs.row( 5); }
    const_ref_type ref_uy()         const { return refs.row( 6); }
    const_ref_type ref_uz()         const { return refs.row( 7); }
    const_ref_type ref_u2()         const { return refs.row( 8); }
    const_ref_type ref_uxux()       const { return refs.row( 9); }
    const_ref_type ref_uxuy()       const { return refs.row(10); }
    const_ref_type ref_uxuz()       const { return refs.row(11); }
    const_ref_type ref_uyuy()       const { return refs.row(12); }
    const_ref_type ref_uyuz()       const { return refs.row(13); }
    const_ref_type ref_uzuz()       const { return refs.row(14); }
    const_ref_type ref_nu()         const { return refs.row(15); }
    const_ref_type ref_nuux()       const { return refs.row(16); }
    const_ref_type ref_nuuy()       const { return refs.row(17); }
    const_ref_type ref_nuuz()       const { return refs.row(18); }
    const_ref_type ref_nuu2()       const { return refs.row(19); }
    const_ref_type ref_nuuxux()     const { return refs.row(20); }
    const_ref_type ref_nuuxuy()     const { return refs.row(21); }
    const_ref_type ref_nuuxuz()     const { return refs.row(22); }
    const_ref_type ref_nuuyuy()     const { return refs.row(23); }
    const_ref_type ref_nuuyuz()     const { return refs.row(24); }
    const_ref_type ref_nuuzuz()     const { return refs.row(25); }
    const_ref_type ref_ex_gradrho() const { return refs.row(26); }
    const_ref_type ref_ey_gradrho() const { return refs.row(27); }
    const_ref_type ref_ez_gradrho() const { return refs.row(28); }
    const_ref_type ref_e_divm()     const { return refs.row(29); }
    const_ref_type ref_e_deltarho() const { return refs.row(30); }
    const_ref_type ref_rhouxux()    const { return refs.row(31); }
    const_ref_type ref_rhouyuy()    const { return refs.row(32); }
    const_ref_type ref_rhouzuz()    const { return refs.row(33); }
    const_ref_type ref_rhoEE()      const { return refs.row(34); }

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld);

    /** @} */

    /**
     * Helper consistently zeroing both \c means and \c refs.
     *
     * Case deliberately differs from the analogous Eigen signature so both
     * consistency with Suzerain naming and so that it is visually distinct.
     */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        means.setZero(Ny, means_type::ColsAtCompileTime);
        implicits.setZero(Ny, implicits_type::ColsAtCompileTime);
        refs.setZero(refs_type::RowsAtCompileTime, Ny);
    }

private:

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    operator_common_block(const operator_common_block&);
    operator_common_block& operator=(const operator_common_block&);
};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_COMMON_BLOCK_HPP */
