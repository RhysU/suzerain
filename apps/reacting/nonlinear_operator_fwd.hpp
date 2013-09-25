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

#include <suzerain/lowstorage.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/reacting_imexop.h>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/largo_definition.hpp>
#include <suzerain/timers.h>

#include "reacting.hpp"
#include "filter_definition.hpp"

namespace suzerain {

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
{
    /** Type of the contiguous storage housing all mean quantities */
    typedef Array<real_t, Dynamic, 20, ColMajor> means_type;


    /** Type of the contiguous storage housing all reference quantities */
    typedef Array<real_t, Dynamic, Dynamic, ColMajor> refs_type;

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

    /** Type returned by the non-const mean quantity accessors. */
    typedef means_type::ColXpr mean_type;

    mean_type       u()                 { return means.col( 0); }
    mean_type       v()                 { return means.col( 1); }
    mean_type       w()                 { return means.col( 2); }
    mean_type       SrhoE()             { return means.col( 3); }
    mean_type       Srhou()             { return means.col( 4); }
    mean_type       Srhov()             { return means.col( 5); }
    mean_type       Srhow()             { return means.col( 6); }
    mean_type       Srho()              { return means.col( 7); }
    mean_type       Srhou_dot_u()       { return means.col( 8); }
    mean_type       fx()                { return means.col( 9); }
    mean_type       fy()                { return means.col(10); }
    mean_type       fz()                { return means.col(11); }
    mean_type       f_dot_u()           { return means.col(12); }
    mean_type       qb()                { return means.col(13); }
    mean_type       CrhoE()             { return means.col(14); }
    mean_type       Crhou()             { return means.col(15); }
    mean_type       Crhov()             { return means.col(16); }
    mean_type       Crhow()             { return means.col(17); }
    mean_type       Crho()              { return means.col(18); }
    mean_type       Crhou_dot_u()       { return means.col(19); }

    /** Type returned by the const mean quantity accessors. */
    typedef means_type::ConstColXpr const_mean_type;

    const_mean_type u()           const { return means.col( 0); }
    const_mean_type v()           const { return means.col( 1); }
    const_mean_type w()           const { return means.col( 2); }
    const_mean_type SrhoE()       const { return means.col( 3); }
    const_mean_type Srhou()       const { return means.col( 4); }
    const_mean_type Srhov()       const { return means.col( 5); }
    const_mean_type Srhow()       const { return means.col( 6); }
    const_mean_type Srho()        const { return means.col( 7); }
    const_mean_type Srhou_dot_u() const { return means.col( 8); }
    const_mean_type fx()          const { return means.col( 9); }
    const_mean_type fy()          const { return means.col(10); }
    const_mean_type fz()          const { return means.col(11); }
    const_mean_type f_dot_u()     const { return means.col(12); }
    const_mean_type qb()          const { return means.col(13); }
    const_mean_type CrhoE()       const { return means.col(14); }
    const_mean_type Crhou()       const { return means.col(15); }
    const_mean_type Crhov()       const { return means.col(16); }
    const_mean_type Crhow()       const { return means.col(17); }
    const_mean_type Crho()        const { return means.col(18); }
    const_mean_type Crhou_dot_u() const { return means.col(19); }

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

    /** Type returned by the non-const reference quantity accessors. */
    typedef refs_type::RowXpr ref_type;

    ref_type       ref_ux()                  { return refs.row( 0       );}
    ref_type       ref_uy()                  { return refs.row( 1       );}
    ref_type       ref_uz()                  { return refs.row( 2       );}
    ref_type       ref_uxuy()                { return refs.row( 3       );}
    ref_type       ref_uzuy()                { return refs.row( 4       );}
    ref_type       ref_p_ru()                { return refs.row( 5       );}
    ref_type       ref_p_rw()                { return refs.row( 6       );}
    ref_type       ref_p_rE()                { return refs.row( 7       );}
    ref_type       ref_vp_ru()               { return refs.row( 8       );}
    ref_type       ref_vp_rw()               { return refs.row( 9       );}
    ref_type       ref_vp_rE()               { return refs.row(10       );}
    ref_type       ref_Cmy_rho()             { return refs.row(11       );}
    ref_type       ref_Ce_rho()              { return refs.row(12       );}
    ref_type       ref_Ce_rv()               { return refs.row(13       );}
    ref_type       ref_nu()                  { return refs.row(14       );}
    ref_type       ref_korCv()               { return refs.row(15       );}
    ref_type       ref_Ds()                  { return refs.row(16       );}
    ref_type       ref_T()                   { return refs.row(17       );}
    ref_type       ref_gamma()               { return refs.row(18       );}
    ref_type       ref_a()                   { return refs.row(19       );}
    ref_type       ref_cs(const int s)       { return refs.row(20     +s);}
    ref_type       ref_es(const int s)       { return refs.row(20+  Ns+s);}
    ref_type       ref_hs(const int s)       { return refs.row(20+2*Ns+s);}

    /** Type returned by the const reference quantity accessors. */
    typedef refs_type::ConstRowXpr const_ref_type;

    const_ref_type ref_ux()            const {return refs.row( 0       );}
    const_ref_type ref_uy()            const {return refs.row( 1       );}
    const_ref_type ref_uz()            const {return refs.row( 2       );}
    const_ref_type ref_uxuy()          const {return refs.row( 3       );}
    const_ref_type ref_uzuy()          const {return refs.row( 4       );}
    const_ref_type ref_p_ru()          const {return refs.row( 5       );}
    const_ref_type ref_p_rw()          const {return refs.row( 6       );}
    const_ref_type ref_p_rE()          const {return refs.row( 7       );}
    const_ref_type ref_vp_ru()         const {return refs.row( 8       );}
    const_ref_type ref_vp_rw()         const {return refs.row( 9       );}
    const_ref_type ref_vp_rE()         const {return refs.row(10       );}
    const_ref_type ref_Cmy_rho()       const {return refs.row(11       );}
    const_ref_type ref_Ce_rho()        const {return refs.row(12       );}
    const_ref_type ref_Ce_rv()         const {return refs.row(13       );}
    const_ref_type ref_nu()            const {return refs.row(14       );}
    const_ref_type ref_korCv()         const {return refs.row(15       );}
    const_ref_type ref_Ds()            const {return refs.row(16       );}
    const_ref_type ref_T()             const {return refs.row(17       );}
    const_ref_type ref_gamma()         const {return refs.row(18       );}
    const_ref_type ref_a()             const {return refs.row(19       );}
    const_ref_type ref_cs(const int s) const {return refs.row(20     +s);}
    const_ref_type ref_es(const int s) const {return refs.row(20+  Ns+s);}
    const_ref_type ref_hs(const int s) const {return refs.row(20+2*Ns+s);}

    std::size_t Nref() const { return 20+3*this->Ns; }

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
        ref.nu         = ref_nu().data();
        ref.korCv      = ref_korCv().data();
        ref.Ds         = ref_Ds().data();
        ref.T          = ref_T().data();
        ref.gamma      = ref_gamma().data();
        ref.a          = ref_a().data();

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
        ld.nu         = inc;
        ld.korCv      = inc;
        ld.Ds         = inc;
        ld.T          = inc;
        ld.gamma      = inc;
        ld.a          = inc;
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
        //refs.setZero(refs_type::RowsAtCompileTime, Ny);
        refs.setZero(this->Nref(), Ny);
    }


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

private:

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    operator_common_block(const operator_common_block&);
    operator_common_block& operator=(const operator_common_block&);
};


/**
 * A complete Navier&ndash;Stokes \c apply_operator implementation.  The
 * implementation is provided as a common building block for
 * <tt>lowstorage::nonlinear_operator< contiguous_state<4,complex_t> ></tt>
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
 *         per lowstorage::nonlinear_operator::apply_operator.
 *
 * @see lowstorage::nonlinear_operator for the (slighly different)
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
            const filter_definition &fsdef,
            support::largo_definition &sgdef,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const ConstitutiveModels& cmods,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method);

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_REACTING_NONLINEAR_OPERATOR_FWD_HPP */
