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

#ifndef SUZERAIN_REACTING_OPERATOR_NONLINEAR_FWD_HPP
#define SUZERAIN_REACTING_OPERATOR_NONLINEAR_FWD_HPP

/** @file
 * Declarations of Nonlinear Navier--Stokes spatial operators
 * implemented within operator_nonlinear.hpp.
 */

#include <suzerain/lowstorage.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/reacting_imexop.h>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timers.h>
#include <suzerain/treatment_constraint.hpp>

#include "reacting.hpp"
#include "definition_filter.hpp"

namespace suzerain {

// Forward declarations
class specification_largo;

namespace reacting {

/** Provides scoping semantics for linearize::type */
namespace linearize {

/** What type of hybrid implicit/explicit linearization is employed? */
enum type {
    none,  ///< No linearization implying a fully implicit treatment
    rhome_y  ///< Linearization of density, momentum, and total energy
};

} // namespace linearize

/** Provides scoping semantics for slowgrowth::type */
namespace slowgrowth {

/** What slow growth sources are employed? */
enum type {
    none   ///< No slow growth sources
};

} // namespace slowgrowth

/** Provides scoping semantics for filter::type */
namespace filter {

/** What filter sources are employed? */
enum type {
    none,
    cook,
    viscous
};

} // namespace filter


/**
 * Storage for holding quantities computed during nonlinear operator
 * application which either are required for linear operator application or for
 * statistics sampling purposes.
 */
class operator_common_block
    : public virtual constraint::treatment::inputs
    , public virtual constraint::treatment::outputs
    , public virtual boost::noncopyable
{

public:

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    operator_common_block() {}

    /**
     * Determines the extent of the implicit treatment
     * by the paired linear and nonlinear operators.
     */
    linearize::type linearization;

    /**
     * Determines the extent of the slow growth treatment
     * by the paired linear and nonlinear operators.
     */
    slowgrowth::type slow_treatment;

    /**
     * Determines the filter operator to be used.
     */
    filter::type filter_treatment;

    /**
     * Number of species
     */
    std::size_t Ns;

    /**
     * Helper consistently zeroing both \c means and \c refs.
     *
     * Case deliberately differs from the analogous Eigen signature so both
     * consistency with Suzerain naming and so that it is visually distinct.
     */
    template<typename Index>
    void set_zero(const Index& Ny)
    {
        means.setZero(Ny, 34);     // Magic number from \ref means columns
        refs .setZero(Nref(), Ny); // Magic number computed by \ref Nref()
    }

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
     * \li \c uu Treated identically to \c u.
     * \li \c uv Treated identically to \c u.
     * \li \c uw Treated identically to \c u.
     * \li \c vv Treated identically to \c u.
     * \li \c vw Treated identically to \c u.
     * \li \c ww Treated identically to \c u.
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
     * "Time-step-specific temporal means" are time averages taken across
     * a single time step of quantities which vary on each substep.
     * As the substeps are all of equal length, a simple running mean
     * is reset on substep zero and then accumulated.
     *
     * Each mean quantity is a single column within \c means.  This facilitates
     * operations across the entire wall-normal profile in a stride one
     * fashion.
     *
     * @see writeups/treatment_channel.tex and writeups/derivation.tex for
     * details on the mean quantities.
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    ArrayXXr means;

    ArrayXXr::ColXpr      u()                  { return means.col( 0); }
    ArrayXXr::ColXpr      v()                  { return means.col( 1); }
    ArrayXXr::ColXpr      w()                  { return means.col( 2); }
    ArrayXXr::ColXpr      uu()                 { return means.col( 3); }
    ArrayXXr::ColXpr      uv()                 { return means.col( 4); }
    ArrayXXr::ColXpr      uw()                 { return means.col( 5); }
    ArrayXXr::ColXpr      vv()                 { return means.col( 6); }
    ArrayXXr::ColXpr      vw()                 { return means.col( 7); }
    ArrayXXr::ColXpr      ww()                 { return means.col( 8); }
    ArrayXXr::ColXpr      SrhoE()              { return means.col( 9); }
    ArrayXXr::ColXpr      Srhou()              { return means.col(10); }
    ArrayXXr::ColXpr      Srhov()              { return means.col(11); }
    ArrayXXr::ColXpr      Srhow()              { return means.col(12); }
    ArrayXXr::ColXpr      Srho()               { return means.col(13); }
    ArrayXXr::ColXpr      Srhou_dot_u()        { return means.col(14); }
    ArrayXXr::ColXpr      fx()                 { return means.col(15); }
    ArrayXXr::ColXpr      fy()                 { return means.col(16); }
    ArrayXXr::ColXpr      fz()                 { return means.col(17); }
    ArrayXXr::ColXpr      f_dot_u()            { return means.col(18); }
    ArrayXXr::ColXpr      qb()                 { return means.col(19); }
    ArrayXXr::ColXpr      CrhoE()              { return means.col(20); }
    ArrayXXr::ColXpr      C2rhoE()             { return means.col(21); }
    ArrayXXr::ColXpr      Crhou()              { return means.col(22); }
    ArrayXXr::ColXpr      C2rhou()             { return means.col(23); }
    ArrayXXr::ColXpr      Crhov()              { return means.col(24); }
    ArrayXXr::ColXpr      C2rhov()             { return means.col(25); }
    ArrayXXr::ColXpr      Crhow()              { return means.col(26); }
    ArrayXXr::ColXpr      C2rhow()             { return means.col(27); }
    ArrayXXr::ColXpr      Crho()               { return means.col(28); }
    ArrayXXr::ColXpr      C2rho()              { return means.col(29); }
    ArrayXXr::ColXpr      Crhou_dot_u()        { return means.col(30); }
    ArrayXXr::ColXpr      C2rhou_dot_u()       { return means.col(31); }
    ArrayXXr::ColXpr      p()                  { return means.col(32); }
    ArrayXXr::ColXpr      p2()                 { return means.col(33); }

    ArrayXXr::ConstColXpr u()            const { return means.col( 0); }
    ArrayXXr::ConstColXpr v()            const { return means.col( 1); }
    ArrayXXr::ConstColXpr w()            const { return means.col( 2); }
    ArrayXXr::ConstColXpr uu()           const { return means.col( 3); }
    ArrayXXr::ConstColXpr uv()           const { return means.col( 4); }
    ArrayXXr::ConstColXpr uw()           const { return means.col( 5); }
    ArrayXXr::ConstColXpr vv()           const { return means.col( 6); }
    ArrayXXr::ConstColXpr vw()           const { return means.col( 7); }
    ArrayXXr::ConstColXpr ww()           const { return means.col( 8); }
    ArrayXXr::ConstColXpr SrhoE()        const { return means.col( 9); }
    ArrayXXr::ConstColXpr Srhou()        const { return means.col(10); }
    ArrayXXr::ConstColXpr Srhov()        const { return means.col(11); }
    ArrayXXr::ConstColXpr Srhow()        const { return means.col(12); }
    ArrayXXr::ConstColXpr Srho()         const { return means.col(13); }
    ArrayXXr::ConstColXpr Srhou_dot_u()  const { return means.col(14); }
    ArrayXXr::ConstColXpr fx()           const { return means.col(15); }
    ArrayXXr::ConstColXpr fy()           const { return means.col(16); }
    ArrayXXr::ConstColXpr fz()           const { return means.col(17); }
    ArrayXXr::ConstColXpr f_dot_u()      const { return means.col(18); }
    ArrayXXr::ConstColXpr qb()           const { return means.col(19); }
    ArrayXXr::ConstColXpr CrhoE()        const { return means.col(20); }
    ArrayXXr::ConstColXpr C2rhoE()       const { return means.col(21); }
    ArrayXXr::ConstColXpr Crhou()        const { return means.col(22); }
    ArrayXXr::ConstColXpr C2rhou()       const { return means.col(23); }
    ArrayXXr::ConstColXpr Crhov()        const { return means.col(24); }
    ArrayXXr::ConstColXpr C2rhov()       const { return means.col(25); }
    ArrayXXr::ConstColXpr Crhow()        const { return means.col(26); }
    ArrayXXr::ConstColXpr C2rhow()       const { return means.col(27); }
    ArrayXXr::ConstColXpr Crho()         const { return means.col(28); }
    ArrayXXr::ConstColXpr C2rho()        const { return means.col(29); }
    ArrayXXr::ConstColXpr Crhou_dot_u()  const { return means.col(30); }
    ArrayXXr::ConstColXpr C2rhou_dot_u() const { return means.col(31); }
    ArrayXXr::ConstColXpr p()            const { return means.col(32); }
    ArrayXXr::ConstColXpr p2()           const { return means.col(33); }

    /** @} */

    // FIXME: Check doxygen to make sure this looks as it is supposed to
    /**
     *
     * The reference quantities stored in \c refs are as follows:
     * \li \c ref_ux         Reference \f$C^{u_x}               \f$
     * \li \c ref_uy         Reference \f$C^{u_y}               \f$
     * \li \c ref_uz         Reference \f$C^{u_z}               \f$
     * \li \c ref_Cmy_rho    Reference \f$C^{my_rho}            \f$
     * \li \c ref_Ce_rho     Reference \f$C^{e_rho}             \f$
     * \li \c ref_Ce_rv      Reference \f$C^{e_rv}              \f$
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

    ArrayXXr::RowXpr      ref_ux()                  {return refs.row( 0       );}
    ArrayXXr::RowXpr      ref_uy()                  {return refs.row( 1       );}
    ArrayXXr::RowXpr      ref_uz()                  {return refs.row( 2       );}
    ArrayXXr::RowXpr      ref_uxux()                {return refs.row( 3       );}
    ArrayXXr::RowXpr      ref_uxuy()                {return refs.row( 4       );}
    ArrayXXr::RowXpr      ref_uxuz()                {return refs.row( 5       );}
    ArrayXXr::RowXpr      ref_uyuy()                {return refs.row( 6       );}
    ArrayXXr::RowXpr      ref_uyuz()                {return refs.row( 7       );}
    ArrayXXr::RowXpr      ref_uzuz()                {return refs.row( 8       );}
    ArrayXXr::RowXpr      ref_p_ru()                {return refs.row( 9       );}
    ArrayXXr::RowXpr      ref_p_rw()                {return refs.row(10       );}
    ArrayXXr::RowXpr      ref_p_rE()                {return refs.row(11       );}
    ArrayXXr::RowXpr      ref_vp_ru()               {return refs.row(12       );}
    ArrayXXr::RowXpr      ref_vp_rw()               {return refs.row(13       );}
    ArrayXXr::RowXpr      ref_vp_rE()               {return refs.row(14       );}
    ArrayXXr::RowXpr      ref_Cmy_rho()             {return refs.row(15       );}
    ArrayXXr::RowXpr      ref_Ce_rho()              {return refs.row(16       );}
    ArrayXXr::RowXpr      ref_Ce_rv()               {return refs.row(17       );}
    ArrayXXr::RowXpr      ref_nu()                  {return refs.row(18       );}
    ArrayXXr::RowXpr      ref_korCv()               {return refs.row(19       );}
    ArrayXXr::RowXpr      ref_Ds()                  {return refs.row(20       );}
    ArrayXXr::RowXpr      ref_T()                   {return refs.row(21       );}
    ArrayXXr::RowXpr      ref_gamma()               {return refs.row(22       );}
    ArrayXXr::RowXpr      ref_a()                   {return refs.row(23       );}
    ArrayXXr::RowXpr      ref_p()                   {return refs.row(24       );}
    ArrayXXr::RowXpr      ref_p2()                  {return refs.row(25       );}
    ArrayXXr::RowXpr      ref_cs(const int s)       {return refs.row(26     +s);}
    ArrayXXr::RowXpr      ref_es(const int s)       {return refs.row(26+  Ns+s);}
    ArrayXXr::RowXpr      ref_hs(const int s)       {return refs.row(26+2*Ns+s);}

    ArrayXXr::ConstRowXpr ref_ux()            const {return refs.row( 0       );}
    ArrayXXr::ConstRowXpr ref_uy()            const {return refs.row( 1       );}
    ArrayXXr::ConstRowXpr ref_uz()            const {return refs.row( 2       );}
    ArrayXXr::ConstRowXpr ref_uxux()          const {return refs.row( 3       );}
    ArrayXXr::ConstRowXpr ref_uxuy()          const {return refs.row( 4       );}
    ArrayXXr::ConstRowXpr ref_uxuz()          const {return refs.row( 5       );}
    ArrayXXr::ConstRowXpr ref_uyuy()          const {return refs.row( 6       );}
    ArrayXXr::ConstRowXpr ref_uyuz()          const {return refs.row( 7       );}
    ArrayXXr::ConstRowXpr ref_uzuz()          const {return refs.row( 8       );}
    ArrayXXr::ConstRowXpr ref_p_ru()          const {return refs.row( 9       );}
    ArrayXXr::ConstRowXpr ref_p_rw()          const {return refs.row(10       );}
    ArrayXXr::ConstRowXpr ref_p_rE()          const {return refs.row(11       );}
    ArrayXXr::ConstRowXpr ref_vp_ru()         const {return refs.row(12       );}
    ArrayXXr::ConstRowXpr ref_vp_rw()         const {return refs.row(13       );}
    ArrayXXr::ConstRowXpr ref_vp_rE()         const {return refs.row(14       );}
    ArrayXXr::ConstRowXpr ref_Cmy_rho()       const {return refs.row(15       );}
    ArrayXXr::ConstRowXpr ref_Ce_rho()        const {return refs.row(16       );}
    ArrayXXr::ConstRowXpr ref_Ce_rv()         const {return refs.row(17       );}
    ArrayXXr::ConstRowXpr ref_nu()            const {return refs.row(18       );}
    ArrayXXr::ConstRowXpr ref_korCv()         const {return refs.row(19       );}
    ArrayXXr::ConstRowXpr ref_Ds()            const {return refs.row(20       );}
    ArrayXXr::ConstRowXpr ref_T()             const {return refs.row(21       );}
    ArrayXXr::ConstRowXpr ref_gamma()         const {return refs.row(22       );}
    ArrayXXr::ConstRowXpr ref_a()             const {return refs.row(23       );}
    ArrayXXr::ConstRowXpr ref_p()             const {return refs.row(24       );}
    ArrayXXr::ConstRowXpr ref_p2()            const {return refs.row(25       );}
    ArrayXXr::ConstRowXpr ref_cs(const int s) const {return refs.row(26     +s);}
    ArrayXXr::ConstRowXpr ref_es(const int s) const {return refs.row(26+  Ns+s);}
    ArrayXXr::ConstRowXpr ref_hs(const int s) const {return refs.row(26+2*Ns+s);}

    std::size_t Nref() const { return 26+3*this->Ns; }

    /** Prepare data for use by implicit operator API in reacting_imexop.h. */
    void imexop_ref(suzerain_reacting_imexop_ref   &ref,
                    suzerain_reacting_imexop_refld &ld)
    {
        // FIXME: update reacting_imexop.h to be consistent with above

        ref.ux         = ref_ux().data();
        ref.uy         = ref_uy().data();
        ref.uz         = ref_uz().data();
        ref.uxuy       = ref_uxuy().data();
        ref.uyuz       = ref_uyuz().data();
        ref.p_ru       = ref_p_ru().data();
        ref.p_rw       = ref_p_rw().data();
        ref.p_rE       = ref_p_rE().data();
        ref.vp_ru      = ref_vp_ru().data();
        ref.vp_rw      = ref_vp_rw().data();
        ref.vp_rE      = ref_vp_rE().data();
        ref.Cmy_rho    = ref_Cmy_rho().data();
        ref.Ce_rho     = ref_Ce_rho().data();
        ref.Ce_rv      = ref_Ce_rv().data();
        ref.nu         = ref_nu().data();
        ref.korCv      = ref_korCv().data();
        ref.Ds         = ref_Ds().data();
        ref.T          = ref_T().data();
        ref.gamma      = ref_gamma().data();
        ref.a          = ref_a().data();
        ref.p          = ref_p().data();
        ref.p2         = ref_p2().data();

        const int inc = refs.colStride();
        ld.ux         = inc;
        ld.uy         = inc;
        ld.uz         = inc;
        ld.uxuy       = inc;
        ld.uyuz       = inc;
        ld.p_ru       = inc;
        ld.p_rw       = inc;
        ld.p_rE       = inc;
        ld.vp_ru      = inc;
        ld.vp_rw      = inc;
        ld.vp_rE      = inc;
        ld.Cmy_rho    = inc;
        ld.Ce_rho     = inc;
        ld.Ce_rv      = inc;
        ld.nu         = inc;
        ld.korCv      = inc;
        ld.Ds         = inc;
        ld.T          = inc;
        ld.gamma      = inc;
        ld.a          = inc;
        ld.p          = inc;
        ld.p2         = inc;
    }

    /** @} */

    /**
     * Vector with species specific total energy.
     */
    // FIXME: Temporarily placing here variables for reference values
    // To be used by Giles BCs.
    // Is this the right place for it?
    real_t   P_ref;
    real_t   T_ref;
    real_t   u_ref;
    real_t   v_ref;
    real_t   w_ref;
    real_t   rho_ref;
    real_t   a_ref;
    real_t   gamma_ref;
    real_t   R_ref;
    real_t   Cv_ref;
    VectorXr cs_ref;
    VectorXr etots_ref;
    VectorXr etots_upper;

    /**
     * Store chemistry sources at the top plane.
     */
     MatrixXXc chemsrcw;

};


/**
 * A complete Navier&ndash;Stokes \c apply_operator implementation.  The
 * implementation is provided as a common building block for
 * <tt>lowstorage::operator_nonlinear< contiguous_state<4,complex_t> ></tt>
 * subclasses allowing varying numbers of passive scalars or varying hybrid
 * implicit/explicit treatment.  Such subclasses feature an overwhelming amount
 * of redundancy and are error prone to create.  This implementation allows
 * writing the "futsy" bits once and then sharing the logic repeatedly.
 * Templating allows for compile-time branching amongst different
 * implementation choices rather than paying runtime cost for such flexibility.
 *
 * Computation follows the "Numerical Considerations" section of
 * writeups/derivation.tex.  No boundary conditions are applied.  The
 * instantaneous wall-normal velocity is averaged across the streamwise and
 * spanwise directions and stored into <tt>common.u()</tt>.
 *
 * \param alpha Parameter controlling bulk viscosity
 *        per \f$\alpha\f$ in namespace \ref rholut.
 * \param beta Temperature power law exponent
 *        per \f$\beta\f$ in namespace \ref rholut.
 * \param gamma Constant ratio of specific heats
 *        per \f$\gamma\f$ in namespace \ref rholut.
 * \param Ma Mach number
 *        per \f$\textrm{Ma}\f$ in namespace \ref rholut.
 * \param Pr Prandtl number
 *        per \f$\textrm{Pr}\f$ in namespace \ref rholut.
 * \param Re Reynolds number
 *        per \f$\textrm{Re}\f$ in namespace \ref rholut.
 * \param o Provides access to discretization and parallel decomposition
 *        operational details.
 * \param common Shared storage for interaction with an linear_operator
 *        implementation providing forcing and boundary conditions.
 * \param fsdef Definitions for filter source.
 * \param sgdef Definitions for slow growth (largo) source.
 * \param msoln If \c msoln evaluates to \c true in a boolean context,
 *        then it will be used to provide manufactured forcing terms.
 * \param time Simulation time at which the operator should be applied.
 *        This allows time-dependent forcing (e.g. from \c msoln).
 * \param swave State to which the operator should be applied.  On
 *        entry, it must be coefficients in the X, Y, and Z directions.
 *        on exit, it must be coefficients in the X and Z directions but
 *        collocation point values in the Y direction.
 * \param method Low-storage timestepping scheme used to compute a stable
 *        time step when <tt>ZerothSubstep == true</tt>.
 *
 * \tparam ZerothSubstep Should one-time activities taking place at the
 *         beginning of a Runge-Kutta step be performed?  Examples include
 *         computing a stable time step size and also computing reference
 *         quantities for linearization.
 * \tparam Linearize What type of hybrid implicit/explicit linearization
 *         is employed?
 * \tparam ManufacturedSolution What manufactured solution should be used to
 *         provide additional forcing (when enabled)?
 *
 * @return A vector of stable timestep sizes according to different criteria
 *         per lowstorage::operator_nonlinear::apply_operator.
 *
 * @see lowstorage::operator_nonlinear for the (slighly different)
 *      interface that an actual operator would provide.
 */
template <bool ZerothSubstep,
          linearize::type Linearize,
          filter::type Filter,
          class ManufacturedSolution,
          class ConstitutiveModels>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const operator_base &o,
            operator_common_block &common,
            const definition_filter &fsdef,
            const specification_largo &sgdef,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const ConstitutiveModels& cmods,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method);

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_REACTING_OPERATOR_NONLINEAR_FWD_HPP */
