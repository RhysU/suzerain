//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_REACTING_NONLINEAR_OPERATOR_FWD_HPP
#define SUZERAIN_REACTING_NONLINEAR_OPERATOR_FWD_HPP

/** @file
 * Declarations of Nonlinear Navier--Stokes spatial operators
 * implemented within nonlinear_operator.hpp.
 */

#include <suzerain/operator_base.hpp>
#include <suzerain/reacting_imexop.h>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timers.h>

#include "reacting.hpp"
#include "filter_definition.hpp"

namespace suzerain {

namespace reacting {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application which either are required for linear operator application or for
 * statistics sampling purposes.
 */
class operator_common_block
{
    /** Type of the contiguous storage housing all mean quantities */
    typedef Array<real_t, Dynamic, 18, ColMajor> means_type;


    // FIXME: Set size correctly
    /** Type of the contiguous storage housing all reference quantities */
    typedef Array<real_t, 14, Dynamic, ColMajor> refs_type;

public:

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    operator_common_block() {}

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
     * \li \c f The \e linear operator accumulates the time-step-specific
     *     temporal mean streamwise (x) component of the implicit \f$f\f$
     *     term in the momentum equation.
     * \li \c f_dot_u The \e linear operator accumulates the
     *     time-step-specific temporal mean of the implicit \f$f\cdot{}u\f$
     *     term in the energy equation.
     * \li \c qb The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$q_b\f$ term in the energy equation.
     * \li \c CrhoE The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}E}\f$ term in
     *     the energy equation.
     * \li \c Crhou The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}u}\f$ term in
     *     the streamwise (x) momentum equation.
     * \li \c Crhov The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}v}\f$ term in
     *     the wall-normal (y) momentum equation.
     * \li \c Crhow The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}w}\f$ term in
     *     the spanwise momentum equation.
     * \li \c Crho The \e linear operator accumulates the time-step-specific
     *     temporal mean of the implicit \f$\mathscr{C}_{\rho{}}\f$ term in the
     *     density equation.
     * \li \c Crhou_dot_u The \e linear operator accumulates the
     *     time-step-specific temporal mean of the implicit
     *     \f$\mathscr{C}_{\rho{}u}\cdot{}u\f$ term in the energy equation.
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
     * @see writeups/channel_treatment.tex and writeups/derivation.tex for
     * details on the mean quantities.
     * @{
     */

    /** Column-major storage housing all mean quantities (one per column). */
    means_type means;

    means_type::ColXpr      u()                 { return means.col( 0); }
    means_type::ColXpr      v()                 { return means.col( 1); }
    means_type::ColXpr      w()                 { return means.col( 2); }
    means_type::ColXpr      SrhoE()             { return means.col( 3); }
    means_type::ColXpr      Srhou()             { return means.col( 4); }
    means_type::ColXpr      Srhov()             { return means.col( 5); }
    means_type::ColXpr      Srhow()             { return means.col( 6); }
    means_type::ColXpr      Srho()              { return means.col( 7); }
    means_type::ColXpr      Srhou_dot_u()       { return means.col( 8); }
    means_type::ColXpr      f()                 { return means.col( 9); }
    means_type::ColXpr      f_dot_u()           { return means.col(10); }
    means_type::ColXpr      qb()                { return means.col(11); }
    means_type::ColXpr      CrhoE()             { return means.col(12); }
    means_type::ColXpr      Crhou()             { return means.col(13); }
    means_type::ColXpr      Crhov()             { return means.col(14); }
    means_type::ColXpr      Crhow()             { return means.col(15); }
    means_type::ColXpr      Crho()              { return means.col(16); }
    means_type::ColXpr      Crhou_dot_u()       { return means.col(17); }

    means_type::ConstColXpr u()           const { return means.col( 0); }
    means_type::ConstColXpr v()           const { return means.col( 1); }
    means_type::ConstColXpr w()           const { return means.col( 2); }
    means_type::ConstColXpr SrhoE()       const { return means.col( 3); }
    means_type::ConstColXpr Srhou()       const { return means.col( 4); }
    means_type::ConstColXpr Srhov()       const { return means.col( 5); }
    means_type::ConstColXpr Srhow()       const { return means.col( 6); }
    means_type::ConstColXpr Srho()        const { return means.col( 7); }
    means_type::ConstColXpr Srhou_dot_u() const { return means.col( 8); }
    means_type::ConstColXpr f()           const { return means.col( 9); }
    means_type::ConstColXpr f_dot_u()     const { return means.col(10); }
    means_type::ConstColXpr qb()          const { return means.col(11); }
    means_type::ConstColXpr CrhoE()       const { return means.col(12); }
    means_type::ConstColXpr Crhou()       const { return means.col(13); }
    means_type::ConstColXpr Crhov()       const { return means.col(14); }
    means_type::ConstColXpr Crhow()       const { return means.col(15); }
    means_type::ConstColXpr Crho()        const { return means.col(16); }
    means_type::ConstColXpr Crhou_dot_u() const { return means.col(17); }

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
    refs_type refs;

    refs_type::RowXpr      ref_ux()               { return refs.row( 0); }
    refs_type::RowXpr      ref_uy()               { return refs.row( 1); }
    refs_type::RowXpr      ref_uz()               { return refs.row( 2); }
    refs_type::RowXpr      ref_uxuy()             { return refs.row( 3); }
    refs_type::RowXpr      ref_uzuy()             { return refs.row( 4); }
    refs_type::RowXpr      ref_p_ru()             { return refs.row( 5); }
    refs_type::RowXpr      ref_p_rw()             { return refs.row( 6); }
    refs_type::RowXpr      ref_p_rE()             { return refs.row( 7); }
    refs_type::RowXpr      ref_vp_ru()            { return refs.row( 8); }
    refs_type::RowXpr      ref_vp_rw()            { return refs.row( 9); }
    refs_type::RowXpr      ref_vp_rE()            { return refs.row(10); }
    refs_type::RowXpr      ref_Cmy_rho()          { return refs.row(11); }
    refs_type::RowXpr      ref_Ce_rho()           { return refs.row(12); }
    refs_type::RowXpr      ref_Ce_rv()            { return refs.row(13); }

    refs_type::ConstRowXpr ref_ux()         const { return refs.row( 0); }
    refs_type::ConstRowXpr ref_uy()         const { return refs.row( 1); }
    refs_type::ConstRowXpr ref_uz()         const { return refs.row( 2); }
    refs_type::ConstRowXpr ref_uxuy()       const { return refs.row( 3); }
    refs_type::ConstRowXpr ref_uzuy()       const { return refs.row( 4); }
    refs_type::ConstRowXpr ref_p_ru()       const { return refs.row( 5); }
    refs_type::ConstRowXpr ref_p_rw()       const { return refs.row( 6); }
    refs_type::ConstRowXpr ref_p_rE()       const { return refs.row( 7); }
    refs_type::ConstRowXpr ref_vp_ru()      const { return refs.row( 8); }
    refs_type::ConstRowXpr ref_vp_rw()      const { return refs.row( 9); }
    refs_type::ConstRowXpr ref_vp_rE()      const { return refs.row(10); }
    refs_type::ConstRowXpr ref_Cmy_rho()    const { return refs.row(11); }
    refs_type::ConstRowXpr ref_Ce_rho()     const { return refs.row(12); }
    refs_type::ConstRowXpr ref_Ce_rv()      const { return refs.row(13); }
    
    /** Prepare data for use by implicit operator API in reacting_imexop.h. */
    void imexop_ref(suzerain_reacting_imexop_ref   &ref,
                    suzerain_reacting_imexop_refld &ld)
    {
        // FIXME: update reacting_imexop.h to be consistent with above

        ref.ux         = ref_ux().data();
        ref.uy         = ref_uy().data();
        ref.uz         = ref_uz().data();
        ref.uxuy       = ref_uxuy().data();
        ref.uzuy       = ref_uzuy().data();
        ref.p_ru       = ref_p_ru().data();
        ref.p_rw       = ref_p_rw().data();
        ref.p_rE       = ref_p_rE().data();
        ref.vp_ru      = ref_vp_ru().data();
        ref.vp_rw      = ref_vp_rw().data();
        ref.vp_rE      = ref_vp_rE().data();
        ref.Cmy_rho    = ref_Cmy_rho().data();
        ref.Ce_rho     = ref_Ce_rho().data();
        ref.Ce_rv      = ref_Ce_rv().data();

        const int inc = refs.colStride();
        ld.ux         = inc;
        ld.uy         = inc;
        ld.uz         = inc;
        ld.uxuy       = inc;
        ld.uzuy       = inc;
        ld.p_ru       = inc;
        ld.p_rw       = inc;
        ld.p_rE       = inc;
        ld.vp_ru      = inc;
        ld.vp_rw      = inc;
        ld.vp_rE      = inc;
        ld.Cmy_rho    = inc;
        ld.Ce_rho     = inc;
        ld.Ce_rv      = inc;
    }

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
        refs.setZero(refs_type::RowsAtCompileTime, Ny);
    }

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    operator_common_block(const operator_common_block&);
    operator_common_block& operator=(const operator_common_block&);
};

/** Provides scoping semantics for linearize::type */
namespace linearize {

/** What type of hybrid implicit/explicit linearization is employed? */
enum type {
    none,  ///< No linearization implying a fully implicit treatment
    rhome  ///< Linearization of density, momentum, and total energy
};

} // namespace linearize

/** Provides scoping semantics for slowgrowth::type */
namespace slowgrowth {

/** What slow growth sources are employed? */
enum type {
    none   ///< No slow growth sources
};

} // namespace slowgrowth

/**
 * A complete Navier&ndash;Stokes \c apply_operator implementation.  The
 * implementation is provided as a common building block for
 * <tt>timestepper::nonlinear_operator< contiguous_state<4,complex_t> ></tt>
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
 * \param msoln If \c msoln evaluates to \c true in a boolean context,
 *        then it will be used to provide manufactured forcing terms.
 * \param time Simulation time at which the operator should be applied.
 *        This allows time-dependent forcing (e.g. from \c msoln).
 * \param swave State to which the operator should be applied.  On
 *        entry, it must be coefficients in the X, Y, and Z directions.
 *        on exit, it must be coefficients in the X and Z directions but
 *        collocation point values in the Y direction.
 * \param evmaxmag_real Maximum real eigenvalue magnitude used for
 *        stable time step computation when <tt>ZerothSubstep == true</tt>.
 * \param evmaxmag_imag Maximum imaginary eigenvalue magnitude used for
 *        stable time step computation when <tt>ZerothSubstep == true</tt>.
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
 *         per timestepper::nonlinear_operator::apply_operator.
 *
 * @see timestepper::nonlinear_operator for the (slighly different)
 *      interface that an actual operator would provide.
 */
template <bool ZerothSubstep,
          linearize::type Linearize,
          class ManufacturedSolution,
	  class ConstitutiveModels>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const operator_base &o,
            operator_common_block &common,
            const filter_definition &fsdef,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const ConstitutiveModels& cmods,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag);

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_REACTING_NONLINEAR_OPERATOR_FWD_HPP */
