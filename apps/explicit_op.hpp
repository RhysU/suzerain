//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// explicit_op.hpp: Operators for channel_explicit
// $Id$

#ifndef EXPLICIT_OP_HPP
#define EXPLICIT_OP_HPP

#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/rholut_imexop.h>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state_fwd.hpp>

#include "precision.hpp"
#include "channel.hpp"

#pragma warning(disable:383 1572)

namespace channel {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application which either are required for linear operator application or for
 * statistics sampling purposes.
 */
class OperatorCommonBlock
{
public:
    // See http://eigen.tuxfamily.org/dox-devel/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OperatorCommonBlock() {}

private:

    /** Type of the contiguous storage housing all mean quantities */
    typedef Eigen::Array<real_t, Eigen::Dynamic,  4, Eigen::ColMajor> means_t;

    /** Type of the contiguous storage housing all reference quantities */
    typedef Eigen::Array<real_t, 13, Eigen::Dynamic, Eigen::ColMajor> refs_t;

public:

    /**
     * The mean quantities stored in \c means are as follows:
     * \li \c u  The nonlinear operator computes the instantaneous spatial
     *     (x, z) mean streamwise velocity profile as collocation point
     *     values.  The linear operator then uses the information to
     *     compute the
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
     * \li \c ref_nu         Reference \f$C^{\nu}               \f$
     * \li \c ref_ux         Reference \f$C^{u_x}               \f$
     * \li \c ref_uy         Reference \f$C^{u_y}               \f$
     * \li \c ref_uz         Reference \f$C^{u_z}               \f$
     * \li \c ref_nuux       Reference \f$C^{\nu u_x}           \f$
     * \li \c ref_nuuy       Reference \f$C^{\nu u_y}           \f$
     * \li \c ref_nuuz       Reference \f$C^{\nu u_z}           \f$
     * \li \c ref_m_gradrho  Reference \f$C^{m}_{\nabla\rho}    \f$
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

    refs_t::RowXpr      ref_nu()               { return refs.row( 0); }
    refs_t::RowXpr      ref_ux()               { return refs.row( 1); }
    refs_t::RowXpr      ref_uy()               { return refs.row( 2); }
    refs_t::RowXpr      ref_uz()               { return refs.row( 3); }
    refs_t::RowXpr      ref_nuux()             { return refs.row( 4); }
    refs_t::RowXpr      ref_nuuy()             { return refs.row( 5); }
    refs_t::RowXpr      ref_nuuz()             { return refs.row( 6); }
    refs_t::RowXpr      ref_m_gradrho()        { return refs.row( 7); }
    refs_t::RowXpr      ref_ex_gradrho()       { return refs.row( 8); }
    refs_t::RowXpr      ref_ey_gradrho()       { return refs.row( 9); }
    refs_t::RowXpr      ref_ez_gradrho()       { return refs.row(10); }
    refs_t::RowXpr      ref_e_divm()           { return refs.row(11); }
    refs_t::RowXpr      ref_e_deltarho()       { return refs.row(12); }

    refs_t::ConstRowXpr ref_nu()         const { return refs.row( 0); }
    refs_t::ConstRowXpr ref_ux()         const { return refs.row( 1); }
    refs_t::ConstRowXpr ref_uy()         const { return refs.row( 2); }
    refs_t::ConstRowXpr ref_uz()         const { return refs.row( 3); }
    refs_t::ConstRowXpr ref_nuux()       const { return refs.row( 4); }
    refs_t::ConstRowXpr ref_nuuy()       const { return refs.row( 5); }
    refs_t::ConstRowXpr ref_nuuz()       const { return refs.row( 6); }
    refs_t::ConstRowXpr ref_m_gradrho()  const { return refs.row( 7); }
    refs_t::ConstRowXpr ref_ex_gradrho() const { return refs.row( 8); }
    refs_t::ConstRowXpr ref_ey_gradrho() const { return refs.row( 9); }
    refs_t::ConstRowXpr ref_ez_gradrho() const { return refs.row(10); }
    refs_t::ConstRowXpr ref_e_divm()     const { return refs.row(11); }
    refs_t::ConstRowXpr ref_e_deltarho() const { return refs.row(12); }

    /** Prepare data for use by implicit operator API in rholut_imexop.h. */
    void imexop_ref(suzerain_rholut_imexop_ref   &ref,
                    suzerain_rholut_imexop_refld &ld)
    {
        ld.nu         = 1; ref.nu         = ref_nu().data();
        ld.ux         = 1; ref.ux         = ref_ux().data();
        ld.uy         = 1; ref.uy         = ref_uy().data();
        ld.uz         = 1; ref.uz         = ref_uz().data();
        ld.nuux       = 1; ref.nuux       = ref_nuux().data();
        ld.nuuy       = 1; ref.nuuy       = ref_nuuy().data();
        ld.nuuz       = 1; ref.nuuz       = ref_nuuz().data();
        ld.m_gradrho  = 1; ref.m_gradrho  = ref_m_gradrho().data();
        ld.ex_gradrho = 1; ref.ex_gradrho = ref_ex_gradrho().data();
        ld.ey_gradrho = 1; ref.ey_gradrho = ref_ey_gradrho().data();
        ld.ez_gradrho = 1; ref.ez_gradrho = ref_ez_gradrho().data();
        ld.e_divm     = 1; ref.e_divm     = ref_e_divm().data();
        ld.e_deltarho = 1; ref.e_deltarho = ref_e_deltarho().data();
    }

    /** @} */

    /** Helper consistently resizing both \c means and \c refs */
    template<typename Index>
    void resize(const Index& Ny)
    {
        means.resize(Ny, means_t::ColsAtCompileTime);
        refs.resize(refs_t::RowsAtCompileTime,  Ny);
    }

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

/**
 * A boundary-condition agnostic, fully explicit Navier&ndash;Stokes operator.
 *
 * During \ref applyOperator the instantaneous wall-normal velocity is averaged
 * across the streamwise and spanwise directions and stored into
 * OperatorCommonBlock::u() using an instance provided at construction time.
 */
class NonlinearOperator
    : public suzerain::OperatorBase<real_t>,
      public suzerain::timestepper::INonlinearOperator<
            suzerain::ContiguousState<4,complex_t>
      >
{
public:

    typedef suzerain::ContiguousState<4,complex_t> state_type;

    NonlinearOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common,
            const boost::shared_ptr<
                  const channel::manufactured_solution>& msoln);

    virtual std::vector<real_t> applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
    {
        // Dispatch to implementation paying nothing for substep-related ifs
        if (substep_index == 0) {
            return applyOperator<true >(time, swave, evmaxmag_real, evmaxmag_imag);
        } else {
            return applyOperator<false>(time, swave, evmaxmag_real, evmaxmag_imag);
        }
    }

protected:

    /** Houses data additionally required for some linear operators */
    OperatorCommonBlock &common;

    /** Holds optional manufactured solution forcing details */
    const boost::shared_ptr<const channel::manufactured_solution> msoln;

private:

    // Internal implementation of INonlinearOperator::applyOperator logic
    template< bool zeroth_substep >
    std::vector<real_t> applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag) const;

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    NonlinearOperator(const NonlinearOperator&);
    NonlinearOperator& operator=(const NonlinearOperator&);

};

/** An operator which applies or inverts a B-spline mass matrix */
class BsplineMassOperator
  : public suzerain::OperatorBase<real_t>,
    public suzerain::timestepper::lowstorage::ILinearOperator<
        suzerain::multi_array::ref<complex_t,4>,
        suzerain::ContiguousState<4,complex_t>
    >
{
public:

    BsplineMassOperator(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop);

    virtual void applyMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void accumulateMassPlusScaledOperator(
             const complex_t &phi,
             const suzerain::multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             suzerain::ContiguousState<4,complex_t> &output,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void invertMassPlusScaledOperator(
             const complex_t &phi,
             suzerain::multi_array::ref<complex_t,4> &state,
             const component delta_t,
             const std::size_t substep_index,
             const real_t iota) const;

private:

     /** Precomputed mass matrix factorization */
    suzerain::bsplineop_luz massluz;

};

/**
 * A mass operator that forces bulk momentum and provides no slip, isothermal
 * walls.  It requires interoperation with NonlinearOperator via
 * OperatorCommonBlock.
 *
 * During \ref invertMassPlusScaledOperator implicit momentum forcing is
 * applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from OperatorCommonBlock::u().
 *
 * Also during \ref invertMassPlusScaledOperator, OperatorCommonBlock::f(),
 * OperatorCommonBlock::f_dot_u(), and OperatorCommonBlock::qb() are
 * accumulated using ILowStorageMethod::iota().
 */
class BsplineMassOperatorIsothermal
  : public BsplineMassOperator
{
public:

    typedef BsplineMassOperator base;

    BsplineMassOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common);

    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            suzerain::multi_array::ref<complex_t,4> &state,
            const component delta_t,
            const std::size_t substep_index,
            const real_t iota) const;

protected:

    /** Precomputed integration coefficients */
    Eigen::VectorXr bulkcoeff;

    /** Houses data required for \ref invertMassPlusScaledOperator */
    OperatorCommonBlock &common;

};

} // namespace channel

#endif  /* EXPLICIT_OP_HPP */
