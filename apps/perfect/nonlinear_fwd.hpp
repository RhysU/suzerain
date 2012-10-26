//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// nonlinear_fwd.hpp: Building blocks for nonlinear Navier--Stokes operators
// $Id$

#ifndef NONLINEAR_FWD_HPP
#define NONLINEAR_FWD_HPP

#include <suzerain/operator_base.hpp>
#include <suzerain/rholut_imexop.h>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timers.h>

#include "perfect.hpp"

namespace suzerain { namespace perfect {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application which either are required for linear operator application or for
 * statistics sampling purposes.
 */
class OperatorCommonBlock
{
public:

    OperatorCommonBlock() {}

private:

    /** Type of the contiguous storage housing all mean quantities */
    typedef Eigen::Array<real_t, Eigen::Dynamic,  4, Eigen::ColMajor> means_t;

    /** Type of the contiguous storage housing all reference quantities */
    typedef Eigen::Array<real_t, 26, Eigen::Dynamic, Eigen::ColMajor> refs_t;

public:

    /**
     * The mean quantities stored in \c means are as follows:
     * \li \c u  The nonlinear operator computes the instantaneous spatial
     *     (x, z) mean streamwise velocity profile as collocation point values.
     *     The linear operator then uses the information to compute the
     *     implicit \f$f\cdot{}u\f$ term in the total energy equation.
     * \li \c f The linear operator accumulates the time-step-specific
     *     temporal mean streamwise (x) component of the implicit \f$f\f$
     *     term in the momentum equation stored as coefficients.
     * \li \c f_dot_u The linear operator accumulates the
     *     time-step-specific temporal mean the implicit \f$f\cdot{}u\f$
     *     term in the energy equation stored as coefficients.
     * \li \c qb The linear operator accumulates the time-step-specific
     *     temporal mean the implicit \f$q_b\f$ term in the energy equation
     *     stored as coefficients.
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
    means_t means;

    means_t::ColXpr      u()             { return means.col(0); }
    means_t::ColXpr      f()             { return means.col(1); }
    means_t::ColXpr      f_dot_u()       { return means.col(2); }
    means_t::ColXpr      qb()            { return means.col(3); }

    means_t::ConstColXpr u()       const { return means.col(0); }
    means_t::ConstColXpr f()       const { return means.col(1); }
    means_t::ConstColXpr f_dot_u() const { return means.col(2); }
    means_t::ConstColXpr qb()      const { return means.col(3); }

    /** @} */

    /**
     *
     * The reference quantities stored in \c refs are as follows:
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
    refs_t refs;

    refs_t::RowXpr      ref_ux()               { return refs.row( 0); }
    refs_t::RowXpr      ref_uy()               { return refs.row( 1); }
    refs_t::RowXpr      ref_uz()               { return refs.row( 2); }
    refs_t::RowXpr      ref_u2()               { return refs.row( 3); }
    refs_t::RowXpr      ref_uxux()             { return refs.row( 4); }
    refs_t::RowXpr      ref_uxuy()             { return refs.row( 5); }
    refs_t::RowXpr      ref_uxuz()             { return refs.row( 6); }
    refs_t::RowXpr      ref_uyuy()             { return refs.row( 7); }
    refs_t::RowXpr      ref_uyuz()             { return refs.row( 8); }
    refs_t::RowXpr      ref_uzuz()             { return refs.row( 9); }
    refs_t::RowXpr      ref_nu()               { return refs.row(10); }
    refs_t::RowXpr      ref_nuux()             { return refs.row(11); }
    refs_t::RowXpr      ref_nuuy()             { return refs.row(12); }
    refs_t::RowXpr      ref_nuuz()             { return refs.row(13); }
    refs_t::RowXpr      ref_nuu2()             { return refs.row(14); }
    refs_t::RowXpr      ref_nuuxux()           { return refs.row(15); }
    refs_t::RowXpr      ref_nuuxuy()           { return refs.row(16); }
    refs_t::RowXpr      ref_nuuxuz()           { return refs.row(17); }
    refs_t::RowXpr      ref_nuuyuy()           { return refs.row(18); }
    refs_t::RowXpr      ref_nuuyuz()           { return refs.row(19); }
    refs_t::RowXpr      ref_nuuzuz()           { return refs.row(20); }
    refs_t::RowXpr      ref_ex_gradrho()       { return refs.row(21); }
    refs_t::RowXpr      ref_ey_gradrho()       { return refs.row(22); }
    refs_t::RowXpr      ref_ez_gradrho()       { return refs.row(23); }
    refs_t::RowXpr      ref_e_divm()           { return refs.row(24); }
    refs_t::RowXpr      ref_e_deltarho()       { return refs.row(25); }

    refs_t::ConstRowXpr ref_ux()         const { return refs.row( 0); }
    refs_t::ConstRowXpr ref_uy()         const { return refs.row( 1); }
    refs_t::ConstRowXpr ref_uz()         const { return refs.row( 2); }
    refs_t::ConstRowXpr ref_u2()         const { return refs.row( 3); }
    refs_t::ConstRowXpr ref_uxux()       const { return refs.row( 4); }
    refs_t::ConstRowXpr ref_uxuy()       const { return refs.row( 5); }
    refs_t::ConstRowXpr ref_uxuz()       const { return refs.row( 6); }
    refs_t::ConstRowXpr ref_uyuy()       const { return refs.row( 7); }
    refs_t::ConstRowXpr ref_uyuz()       const { return refs.row( 8); }
    refs_t::ConstRowXpr ref_uzuz()       const { return refs.row( 9); }
    refs_t::ConstRowXpr ref_nu()         const { return refs.row(10); }
    refs_t::ConstRowXpr ref_nuux()       const { return refs.row(11); }
    refs_t::ConstRowXpr ref_nuuy()       const { return refs.row(12); }
    refs_t::ConstRowXpr ref_nuuz()       const { return refs.row(13); }
    refs_t::ConstRowXpr ref_nuu2()       const { return refs.row(14); }
    refs_t::ConstRowXpr ref_nuuxux()     const { return refs.row(15); }
    refs_t::ConstRowXpr ref_nuuxuy()     const { return refs.row(16); }
    refs_t::ConstRowXpr ref_nuuxuz()     const { return refs.row(17); }
    refs_t::ConstRowXpr ref_nuuyuy()     const { return refs.row(18); }
    refs_t::ConstRowXpr ref_nuuyuz()     const { return refs.row(19); }
    refs_t::ConstRowXpr ref_nuuzuz()     const { return refs.row(20); }
    refs_t::ConstRowXpr ref_ex_gradrho() const { return refs.row(21); }
    refs_t::ConstRowXpr ref_ey_gradrho() const { return refs.row(22); }
    refs_t::ConstRowXpr ref_ez_gradrho() const { return refs.row(23); }
    refs_t::ConstRowXpr ref_e_divm()     const { return refs.row(24); }
    refs_t::ConstRowXpr ref_e_deltarho() const { return refs.row(25); }

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld)
    {
        ref.ux         = ref_ux().data();
        ref.uy         = ref_uy().data();
        ref.uz         = ref_uz().data();
        ref.u2         = ref_u2().data();
        ref.uxux       = ref_uxux().data();
        ref.uxuy       = ref_uxuy().data();
        ref.uxuz       = ref_uxuz().data();
        ref.uyuy       = ref_uyuy().data();
        ref.uyuz       = ref_uyuz().data();
        ref.uzuz       = ref_uzuz().data();
        ref.nu         = ref_nu().data();
        ref.nuux       = ref_nuux().data();
        ref.nuuy       = ref_nuuy().data();
        ref.nuuz       = ref_nuuz().data();
        ref.nuu2       = ref_nuu2().data();
        ref.nuuxux     = ref_nuuxux().data();
        ref.nuuxuy     = ref_nuuxuy().data();
        ref.nuuxuz     = ref_nuuxuz().data();
        ref.nuuyuy     = ref_nuuyuy().data();
        ref.nuuyuz     = ref_nuuyuz().data();
        ref.nuuzuz     = ref_nuuzuz().data();
        ref.ex_gradrho = ref_ex_gradrho().data();
        ref.ey_gradrho = ref_ey_gradrho().data();
        ref.ez_gradrho = ref_ez_gradrho().data();
        ref.e_divm     = ref_e_divm().data();
        ref.e_deltarho = ref_e_deltarho().data();

        const int inc = refs.colStride();
        ld.ux         = inc;
        ld.uy         = inc;
        ld.uz         = inc;
        ld.u2         = inc;
        ld.uxux       = inc;
        ld.uxuy       = inc;
        ld.uxuz       = inc;
        ld.uyuy       = inc;
        ld.uyuz       = inc;
        ld.uzuz       = inc;
        ld.nu         = inc;
        ld.nuux       = inc;
        ld.nuuy       = inc;
        ld.nuuz       = inc;
        ld.nuu2       = inc;
        ld.nuuxux     = inc;
        ld.nuuxuy     = inc;
        ld.nuuxuz     = inc;
        ld.nuuyuy     = inc;
        ld.nuuyuz     = inc;
        ld.nuuzuz     = inc;
        ld.ex_gradrho = inc;
        ld.ey_gradrho = inc;
        ld.ez_gradrho = inc;
        ld.e_divm     = inc;
        ld.e_deltarho = inc;
    }

    /** @} */

    /** Helper consistently zeroing both \c means and \c refs */
    template<typename Index>
    void setZero(const Index& Ny)
    {
        means.setZero(Ny, means_t::ColsAtCompileTime);
        refs.setZero(refs_t::RowsAtCompileTime, Ny);
    }

private:

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    OperatorCommonBlock(const OperatorCommonBlock&);
    OperatorCommonBlock& operator=(const OperatorCommonBlock&);
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
 * A complete Navier&ndash;Stokes \c applyOperator implementation.  The
 * implementation is provided as a common building block for
 * <tt>timestepper::INonlinearOperator< ContiguousState<4,complex_t> ></tt>
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
 * \param common Shared storage for interaction with an ILinearOperator
 *        implementation providing forcing and boundary conditions.
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
 *         per timestepper::INonlinearOperator::applyOperator.
 *
 * @see timestepper::INonlinearOperator for the (slighly different)
 *      interface that an actual operator would provide.
 */
template <bool ZerothSubstep,
          linearize::type Linearize,
          class ManufacturedSolution>
std::vector<real_t> applyNonlinearOperator(
            const real_t alpha,
            const real_t beta,
            const real_t gamma,
            const real_t Ma,
            const real_t Pr,
            const real_t Re,
            const OperatorBase &o,
            OperatorCommonBlock &common,
            const boost::shared_ptr<const ManufacturedSolution>& msoln,
            const real_t time,
            ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag);

} /* namespace perfect */ } /* namespace suzerain */

#endif  /* NONLINEAR_FWD_HPP */
