/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * rholut_imexop.h: hybrid implicit/explicit operator apply and solve
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_RHOLUT_IMEXOP_H__
#define __SUZERAIN_RHOLUT_IMEXOP_H__

#include <suzerain/bsmbsm.h>
#include <suzerain/bsplineop.h>

/** @file
 * Provides implicit operator apply and solve routines for a single wall-normal
 * pencil of state information.  Meant to be used in conjunction with compute
 * kernels found in rholut.hpp.  Coded in C99 to facilitate profiling.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Scenario-like constants used for operator formation. */
typedef struct {
    const double Re;    /**< \f$\mbox{Re} = \frac{\rho_0 u_0 l_0 }{\mu_0}\f$ */
    const double Pr;    /**< \f$\mbox{Pr} = \frac{\mu_0 C_{p}}{\kappa_0}\f$ */
    const double Ma;    /**< \f$\mbox{Ma} = \frac{u_0}{a_0}\f$ */
    const double alpha; /**< \f$\alpha\f$ such that \f$\mu_B = \alpha\mu\f$ */
    const double gamma; /**< \f$\gamma\f$ such that \f$\rho T = p \gamma \f$ */
} suzerain_rholut_imexop_scenario;

/** Wall-normal diagonal reference matrices used for operator formation. */
typedef struct {
    double *nu;         /**< \f$ C^{\nu}                \f$ */
    double *ux;         /**< \f$ C^{u_x}                \f$ */
    double *uy;         /**< \f$ C^{u_y}                \f$ */
    double *uz;         /**< \f$ C^{u_z}                \f$ */
    double *nuux;       /**< \f$ C^{\nu u_x}            \f$ */
    double *nuuy;       /**< \f$ C^{\nu u_y}            \f$ */
    double *nuuz;       /**< \f$ C^{\nu u_z}            \f$ */
    double *m_gradrho;  /**< \f$ C^{m}_{\nabla\rho}     \f$ */
    double *ex_gradrho; /**< \f$ C^{e_x}_{\nabla\rho}   \f$ */
    double *ey_gradrho; /**< \f$ C^{e_y}_{\nabla\rho}   \f$ */
    double *ez_gradrho; /**< \f$ C^{e_z}_{\nabla\rho}   \f$ */
    double *e_divm;     /**< \f$ C^{e}_{\nabla\cdot{}m} \f$ */
    double *e_deltarho; /**< \f$ C^{e}_{\Delta\rho}     \f$ */
} suzerain_rholut_imexop_ref;

/** Strides between elements in a \ref suzerain_rholut_imexop_ref. */
typedef struct {
    int nu;             /**< \copydoc suzerain_rholut_imexop_ref::nu         */
    int ux;             /**< \copydoc suzerain_rholut_imexop_ref::ux         */
    int uy;             /**< \copydoc suzerain_rholut_imexop_ref::uy         */
    int uz;             /**< \copydoc suzerain_rholut_imexop_ref::uz         */
    int nuux;           /**< \copydoc suzerain_rholut_imexop_ref::nuux       */
    int nuuy;           /**< \copydoc suzerain_rholut_imexop_ref::nuuy       */
    int nuuz;           /**< \copydoc suzerain_rholut_imexop_ref::nuuz       */
    int m_gradrho;      /**< \copydoc suzerain_rholut_imexop_ref::m_gradrho  */
    int ex_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ex_gradrho */
    int ey_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ey_gradrho */
    int ez_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ez_gradrho */
    int e_divm;         /**< \copydoc suzerain_rholut_imexop_ref::e_divm     */
    int e_deltarho;     /**< \copydoc suzerain_rholut_imexop_ref::e_deltarho */
} suzerain_rholut_imexop_refld;

/**
 * Apply linear implicit operator \f$y \leftarrow{} \left(M + \varphi{}L\right)
 * x + \beta{} y\f$.  The matrix \f$L\f$ is a function of the provided
 * wavenumbers, scenario parameters, and wall-normal reference quantities.  The
 * problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  Input and output state variables must be stride one.
 *
 * @param[in] phi Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in] km X direction wavenumber \f$k_m = 2\pi{}m/L_x\f$
 *               (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] kn Z direction wavenumber \f$k_n = 2\pi{}n/L_z\f$
 *               (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] s  Scenario parameters used to form the operator.
 * @param[in] r  Reference quantities used to form the operator.
 * @param[in] ld Strides between reference quantity values.
 * @param[in] w  B-spline workspace providing discrete operators.
 * @param[in] in_rho    Wall-normal input data for \f$\rho{}\f$.
 * @param[in] in_rhou   Wall-normal input data for \f$\rho{}u\f$.
 * @param[in] in_rhov   Wall-normal input data for \f$\rho{}v\f$.
 * @param[in] in_rhow   Wall-normal input data for \f$\rho{}w\f$.
 * @param[in] in_rhoe   Wall-normal input data for \f$\rho{}e\f$.
 * @param[in] beta      Accumulation scaling factor \f$\beta\f$.
 * @param[out] out_rho  Wall-normal output data for \f$\rho{}\f$.
 * @param[out] out_rhou Wall-normal output data for \f$\rho{}u\f$.
 * @param[out] out_rhov Wall-normal output data for \f$\rho{}v\f$.
 * @param[out] out_rhow Wall-normal output data for \f$\rho{}w\f$.
 * @param[out] out_rhoe Wall-normal output data for \f$\rho{}e\f$.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for full details.
 */
void
suzerain_rholut_imexop_apply(
        const double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const double (*in_rho )[2],
        const double (*in_rhou)[2],
        const double (*in_rhov)[2],
        const double (*in_rhow)[2],
        const double (*in_rhoe)[2],
        const double beta,
        double (*out_rho )[2],
        double (*out_rhou)[2],
        double (*out_rhov)[2],
        double (*out_rhow)[2],
        double (*out_rhoe)[2]);

/**
 * Pack \f$\left(M + \varphi{}L\right)\f$ into the corresponding locations
 * within contiguous storage of \f$P\left(M + \varphi{}L\right)P^{\mbox{T}}\f$
 * for use by the BLAS' <tt>gbmv</tt> or LAPACK's <tt>gbsvx</tt>.  The matrix
 * \f$L\f$ is a function of the provided wavenumbers, scenario parameters, and
 * wall-normal reference quantities.
 *
 * The problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  On entry, \c papt must be of at least size
 * <tt>w->n*(5*(w->max_kl + w->max_ku + 2) - 1)</tt>.  Boundary conditions,
 * which are \em not applied, will require using information about the
 * permutation returned in \c A.
 *
 * @param[in] phi   Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in] km    X direction wavenumber \f$k_m = 2\pi{}m/L_x\f$
 *                  (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] kn    Z direction wavenumber \f$k_n = 2\pi{}n/L_z\f$
 *                  (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] s     Scenario parameters used to form the operator.
 * @param[in] r     Reference quantities used to form the operator.
 * @param[in] ld    Strides between reference quantity values.
 * @param[in] w     B-spline workspace providing discrete operators.
 * @param[in] ndx_rho  Order of contiguous density data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhou Order of contiguous streamwise momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhov Order of contiguous wall-normal momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhow Order of contiguous spanwise momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhoe Order of contiguous total energy data within a
 *                     globally contiguous state vector.
 * @param[out] A    Storage details for the BSMBSM matrix \c papt.
 * @param[out] papt Band storage of renumbered matrix \f$PAP^{\mbox{T}}\f$
 *                  which will have <tt>A->KL</tt> and <tt>A->KU</tt>
 *                  diagonals and leading dimension <tt>A->LD</tt>.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for full details.
 * @see suzerain_bsmbsm_zaPxpby() for how to permute state to match \c papt.
 */
void
suzerain_rholut_imexop_packc(
        const double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int ndx_rho,
        const int ndx_rhou,
        const int ndx_rhov,
        const int ndx_rhow,
        const int ndx_rhoe,
        suzerain_bsmbsm *A,
        double (*papt)[2]);

/**
 * Pack \f$\left(M + \varphi{}L\right)\f$ into the corresponding locations
 * within contiguous, LU factorization-ready storage of \f$P\left(M +
 * \varphi{}L\right)P^{\mbox{T}}\f$ for use by the LAPACK's <tt>gbtrf</tt> or
 * <tt>gbsv</tt>.  The matrix \f$L\f$ is a function of the provided
 * wavenumbers, scenario parameters, and wall-normal reference quantities.
 *
 * The problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  On entry, \c papt must be of at least size
 * <tt>w->n*(10*(w->max_kl+1) + 5*(w->max_ku+1) - 2)</tt>.  Boundary
 * conditions, which are \em not applied, will require using information about
 * the permutation returned in \c A taking care that in accordance with
 * <tt>gbtrf</tt> the operator starts at row <tt>A->KL</tt>.
 *
 * @param[in] phi   Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in] km    X direction wavenumber \f$k_m = 2\pi{}m/L_x\f$
 *                  (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] kn    Z direction wavenumber \f$k_n = 2\pi{}n/L_z\f$
 *                  (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in] s     Scenario parameters used to form the operator.
 * @param[in] r     Reference quantities used to form the operator.
 * @param[in] ld    Strides between reference quantity values.
 * @param[in] w     B-spline workspace providing discrete operators.
 * @param[in] ndx_rho  Order of contiguous density data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhou Order of contiguous streamwise momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhov Order of contiguous wall-normal momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhow Order of contiguous spanwise momentum data within a
 *                     globally contiguous state vector.
 * @param[in] ndx_rhoe Order of contiguous total energy data within a
 *                     globally contiguous state vector.
 * @param[out] A    Storage details for the BSMBSM matrix \c papt.
 * @param[out] papt Band storage of renumbered matrix \f$PAP^{\mbox{T}}\f$
 *                  which will have <tt>A->KL</tt> and <tt>A->KU</tt>
 *                  diagonals and leading dimension <tt>A->LD + A->KL</tt>.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for full details.
 * @see suzerain_bsmbsm_zaPxpby() for how to permute state to match \c papt.
 */
void
suzerain_rholut_imexop_packf(
        const double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int ndx_rho,
        const int ndx_rhou,
        const int ndx_rhov,
        const int ndx_rhow,
        const int ndx_rhoe,
        suzerain_bsmbsm *A,
        double (*papt)[2]);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_RHOLUT_IMEXOP_H__ */
