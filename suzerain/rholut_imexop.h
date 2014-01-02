/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#ifndef SUZERAIN_RHOLUT_IMEXOP_H
#define SUZERAIN_RHOLUT_IMEXOP_H

/** @file
 * Provides implicit operator apply and solve routines for a single wall-normal
 * pencil of state information.  The operator is described in "Hybrid
 * implicit/explicit" treatment within the model document.  These routines are
 * meant to be used in conjunction with compute kernels found in rholut.hpp.
 * The logic is coded in C99 to facilitate profiling.
 */

#include <suzerain/bsmbsm.h>
#include <suzerain/bsplineop.h>
#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Scenario-like constants used for operator formation. */
typedef struct {
    double Re;    /**< \f$\mbox{Re} = \frac{\rho_0 u_0 l_0 }{\mu_0}\f$ */
    double Pr;    /**< \f$\mbox{Pr} = \frac{\mu_0 C_{p}}{\kappa_0}\f$ */
    double Ma;    /**< \f$\mbox{Ma} = \frac{u_0}{a_0}\f$ */
    double alpha; /**< \f$\alpha\f$ such that \f$\mu_B = \alpha\mu\f$ */
    double gamma; /**< \f$\gamma\f$ such that \f$\rho T = p \gamma \f$ */
} suzerain_rholut_imexop_scenario;

/** Wall-normal diagonal reference matrices used for operator formation. */
typedef struct {
    double *ux;         /**< \f$ C^{u_x}                \f$ */
    double *uy;         /**< \f$ C^{u_y}                \f$ */
    double *uz;         /**< \f$ C^{u_z}                \f$ */
    double *u2;         /**< \f$ C^{u^2}                \f$ */
    double *uxux;       /**< \f$ C^{u_x u_x}            \f$ */
    double *uxuy;       /**< \f$ C^{u_x u_y}            \f$ */
    double *uxuz;       /**< \f$ C^{u_x u_z}            \f$ */
    double *uyuy;       /**< \f$ C^{u_y u_y}            \f$ */
    double *uyuz;       /**< \f$ C^{u_y u_z}            \f$ */
    double *uzuz;       /**< \f$ C^{u_z u_z}            \f$ */
    double *nu;         /**< \f$ C^{\nu}                \f$ */
    double *nuux;       /**< \f$ C^{\nu u_x}            \f$ */
    double *nuuy;       /**< \f$ C^{\nu u_y}            \f$ */
    double *nuuz;       /**< \f$ C^{\nu u_z}            \f$ */
    double *nuu2;       /**< \f$ C^{\nu u^2}            \f$ */
    double *nuuxux;     /**< \f$ C^{\nu u_x u_x}        \f$ */
    double *nuuxuy;     /**< \f$ C^{\nu u_x u_y}        \f$ */
    double *nuuxuz;     /**< \f$ C^{\nu u_x u_z}        \f$ */
    double *nuuyuy;     /**< \f$ C^{\nu u_y u_y}        \f$ */
    double *nuuyuz;     /**< \f$ C^{\nu u_y u_z}        \f$ */
    double *nuuzuz;     /**< \f$ C^{\nu u_z u_z}        \f$ */
    double *ex_gradrho; /**< \f$ C^{e_x}_{\nabla\rho}   \f$ */
    double *ey_gradrho; /**< \f$ C^{e_y}_{\nabla\rho}   \f$ */
    double *ez_gradrho; /**< \f$ C^{e_z}_{\nabla\rho}   \f$ */
    double *e_divm;     /**< \f$ C^{e}_{\nabla\cdot{}m} \f$ */
    double *e_deltarho; /**< \f$ C^{e}_{\Delta\rho}     \f$ */
} suzerain_rholut_imexop_ref;

/** Strides between elements in a \ref suzerain_rholut_imexop_ref. */
typedef struct {
    int ux;             /**< \copydoc suzerain_rholut_imexop_ref::ux         */
    int uy;             /**< \copydoc suzerain_rholut_imexop_ref::uy         */
    int uz;             /**< \copydoc suzerain_rholut_imexop_ref::uz         */
    int u2;             /**< \copydoc suzerain_rholut_imexop_ref::u2         */
    int uxux;           /**< \copydoc suzerain_rholut_imexop_ref::uxux       */
    int uxuy;           /**< \copydoc suzerain_rholut_imexop_ref::uxuy       */
    int uxuz;           /**< \copydoc suzerain_rholut_imexop_ref::uxuz       */
    int uyuy;           /**< \copydoc suzerain_rholut_imexop_ref::uyuy       */
    int uyuz;           /**< \copydoc suzerain_rholut_imexop_ref::uyuz       */
    int uzuz;           /**< \copydoc suzerain_rholut_imexop_ref::uzuz       */
    int nu;             /**< \copydoc suzerain_rholut_imexop_ref::nu         */
    int nuux;           /**< \copydoc suzerain_rholut_imexop_ref::nuux       */
    int nuuy;           /**< \copydoc suzerain_rholut_imexop_ref::nuuy       */
    int nuuz;           /**< \copydoc suzerain_rholut_imexop_ref::nuuz       */
    int nuu2;           /**< \copydoc suzerain_rholut_imexop_ref::nuu2       */
    int nuuxux;         /**< \copydoc suzerain_rholut_imexop_ref::nuuxux     */
    int nuuxuy;         /**< \copydoc suzerain_rholut_imexop_ref::nuuxuy     */
    int nuuxuz;         /**< \copydoc suzerain_rholut_imexop_ref::nuuxuz     */
    int nuuyuy;         /**< \copydoc suzerain_rholut_imexop_ref::nuuyuy     */
    int nuuyuz;         /**< \copydoc suzerain_rholut_imexop_ref::nuuyuz     */
    int nuuzuz;         /**< \copydoc suzerain_rholut_imexop_ref::nuuzuz     */
    int ex_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ex_gradrho */
    int ey_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ey_gradrho */
    int ez_gradrho;     /**< \copydoc suzerain_rholut_imexop_ref::ez_gradrho */
    int e_divm;         /**< \copydoc suzerain_rholut_imexop_ref::e_divm     */
    int e_deltarho;     /**< \copydoc suzerain_rholut_imexop_ref::e_deltarho */
} suzerain_rholut_imexop_refld;

/**
 * Accumulate the linear implicit operator application \f$y \leftarrow{}
 * \left(M + \varphi{}L\right) x + \beta{} y\f$.  The matrix \f$L\f$ is a
 * function of the provided wavenumbers, scenario parameters, and wall-normal
 * reference quantities.  The problem size and discrete operators are taken
 * from the provided B-spline workspace \c w.  Input and output state variables
 * must be stride one.
 *
 * Providing a \c NULL value for any input and its corresponding output (e.g.
 * \c in_rho and \c out_rho) omits the corresponding part of the operator.
 * This may be useful for omitting part of the operator (e.g. the density
 * equation).
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
 * @param[in]  in_rho_E  Wall-normal input data for \f$\rho{}E\f$.
 * @param[in]  in_rho_w  Wall-normal input data for \f$\rho{}w\f$.
 * @param[in]  in_rho_v  Wall-normal input data for \f$\rho{}v\f$.
 * @param[in]  in_rho_u  Wall-normal input data for \f$\rho{}u\f$.
 * @param[in]  in_rho    Wall-normal input data for \f$\rho{}\f$.
 * @param[in]  beta      Accumulation scaling factor \f$\beta\f$.
 * @param[out] out_rho_E Wall-normal output data for \f$\rho{}E\f$.
 * @param[out] out_rho_w Wall-normal output data for \f$\rho{}w\f$.
 * @param[out] out_rho_v Wall-normal output data for \f$\rho{}v\f$.
 * @param[out] out_rho_u Wall-normal output data for \f$\rho{}u\f$.
 * @param[out] out_rho   Wall-normal output data for \f$\rho{}\f$.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for details.
 */
void
suzerain_rholut_imexop_accumulate(
        const complex_double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const complex_double *in_rho_E,
        const complex_double *in_rho_w,
        const complex_double *in_rho_v,
        const complex_double *in_rho_u,
        const complex_double *in_rho,
        const complex_double beta,
        complex_double *out_rho_E,
        complex_double *out_rho_u,
        complex_double *out_rho_v,
        complex_double *out_rho_w,
        complex_double *out_rho );

/**
 * Accumulate the ``zero-zero'' mode 'linear implicit operator application \f$y
 * \leftarrow{} \left(M + \varphi{}L\right) x + \beta{} y\f$.  This method is
 * equivalent to calling suzerain_rholut_imexop_accumulate() using <tt>km ==
 * 0</tt> and <tt>kn == 0</tt> but should be appreciably faster.
 *
 * @see The documentation of suzerain_rholut_imexop_accumulate() for details.
 */
void
suzerain_rholut_imexop_accumulate00(
        const complex_double phi,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const complex_double *in_rho_E,
        const complex_double *in_rho_w,
        const complex_double *in_rho_v,
        const complex_double *in_rho_u,
        const complex_double *in_rho,
        const complex_double beta,
        complex_double *out_rho_E,
        complex_double *out_rho_u,
        complex_double *out_rho_v,
        complex_double *out_rho_w,
        complex_double *out_rho );

/**
 * Pack \f$\left(M + \varphi{}L\right)^{\mbox{T}}\f$ into the corresponding
 * locations within contiguous storage of \f$P\left(M +
 * \varphi{}L\right)^{\mbox{T}}P^{\mbox{T}}\f$ for use by the BLAS'
 * <tt>gbmv</tt> or LAPACK's <tt>gbsvx</tt>.  The matrix \f$L\f$ is a function
 * of the provided wavenumbers, scenario parameters, and wall-normal reference
 * quantities.  Notice that the \e transpose of the operator is assembled.
 *
 * Supplying a negative value for state order parameter (i.e. <tt>rho</tt>,
 * <tt>rho_u</tt>, <tt>rho_v</tt>, <tt>rho_w</tt>, and <tt>rho_E</tt>) will
 * omit the corresponding rows and columns from the formed matrix.  This may be
 * useful for omitting part of the operator (e.g. the density equation).  All
 * nonnegative order parameters must form a unique, contiguous set starting
 * from zero.
 *
 * The problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  On entry, \c patpt must be of at least size
 * <tt>(nneg*w->n)*(nneg*(w->max_kl + w->max_ku + 2) - 1)</tt> where \c nneg
 * is the number of nonnegative state order parameters.  Boundary conditions,
 * which are \em not applied, will require using information about the
 * permutation returned in \c A_T.
 *
 * @param[in]  phi   Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in]  km    X direction wavenumber \f$k_m = 2\pi{}m/L_x\f$
 *                   (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in]  kn    Z direction wavenumber \f$k_n = 2\pi{}n/L_z\f$
 *                   (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in]  s     Scenario parameters used to form the operator.
 * @param[in]  r     Reference quantities used to form the operator.
 * @param[in]  ld    Strides between reference quantity values.
 * @param[in]  w     B-spline workspace providing discrete operators.
 * @param[in]  rho_E Order of contiguous total energy data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_u Order of contiguous streamwise momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_v Order of contiguous wall-normal momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_w Order of contiguous spanwise momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho   Order of contiguous density data within a
 *                   globally contiguous state vector.
 * @param[in]  buf   Working storage of size at least
 *                   <tt>w->n*(w->max_kl + 1 + w->max_ku)</tt>.
 * @param[out] A_T   Storage details for the BSMBSM matrix \c patpt.
 * @param[out] patpt Band storage of renumbered matrix
 *                   \f$PA^{\mbox{T}}P^{\mbox{T}}\f$ which will have
 *                   <tt>A_T->KL</tt> and <tt>A_T->KU</tt> diagonals and
 *                   leading dimension <tt>A_T->LD</tt>.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for details.
 * @see suzerain_bsmbsm_zaPxpby() for how to permute state to match \c patpt.
 */
void
suzerain_rholut_imexop_packc(
        const complex_double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int rho_E,
        const int rho_w,
        const int rho_v,
        const int rho_u,
        const int rho,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double * const patpt);

/**
 * Pack the ``zero-zero'' mode portion of '\f$\left(M +
 * \varphi{}L\right)^{\mbox{T}}\f$ into the corresponding locations within
 * contiguous storage of \f$P\left(M +
 * \varphi{}L\right)^{\mbox{T}}P^{\mbox{T}}\f$ for use by the BLAS'
 * <tt>gbmv</tt> or LAPACK's <tt>gbsvx</tt>.  This method is equivalent to
 * calling suzerain_rholut_imexop_packc() using <tt>km == 0</tt> and <tt>kn ==
 * 0</tt> but should be appreciably faster.
 *
 * @see The documentation of suzerain_rholut_imexop_packc() for details.
 */
void
suzerain_rholut_imexop_packc00(
        const complex_double phi,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int rho_E,
        const int rho_w,
        const int rho_v,
        const int rho_u,
        const int rho,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double * const patpt);

/**
 * Pack \f$\left(M + \varphi{}L\right)^{\mbox{T}}\f$ into the corresponding
 * locations within contiguous, LU factorization-ready storage of \f$P\left(M +
 * \varphi{}L\right)^{\mbox{T}}P^{\mbox{T}}\f$ for use by the LAPACK's
 * <tt>gbtrf</tt> or <tt>gbsv</tt>.  The matrix \f$L\f$ is a function of the
 * provided wavenumbers, scenario parameters, and wall-normal reference
 * quantities.  Notice that the \e transpose of the operator is assembled.
 *
 * Supplying a negative value for state order parameter (i.e. <tt>rho</tt>,
 * <tt>rho_u</tt>, <tt>rho_v</tt>, <tt>rho_w</tt>, and <tt>rho_E</tt>) will
 * omit the corresponding rows and columns from the formed matrix.  This may be
 * useful for omitting part of the operator (e.g. the density equation).  All
 * nonnegative order parameters must form a unique, contiguous set starting
 * from zero.
 *
 * The problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  On entry, \c patpt must be of at least size
 * <tt>(nneg*w->n)*(nneg*(2*w->max_kl + w->max_ku + 3) - 2)</tt> where
 * <tt>nneg</tt> is the number of nonnegative state order parameters.  Boundary
 * conditions, which are \em not applied, will require using information about
 * the permutation returned in \c A_T taking care that in accordance with
 * <tt>gbtrf</tt> the operator starts at row <tt>A_T->KL</tt>.
 *
 * @param[in]  phi   Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in]  km    X direction wavenumber \f$k_m = 2\pi{}m/L_x\f$
 *                   (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in]  kn    Z direction wavenumber \f$k_n = 2\pi{}n/L_z\f$
 *                   (from, for example, \ref suzerain_inorder_wavenumber).
 * @param[in]  s     Scenario parameters used to form the operator.
 * @param[in]  r     Reference quantities used to form the operator.
 * @param[in]  ld    Strides between reference quantity values.
 * @param[in]  w     B-spline workspace providing discrete operators.
 * @param[in]  rho_E Order of contiguous total energy data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_u Order of contiguous streamwise momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_v Order of contiguous wall-normal momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho_w Order of contiguous spanwise momentum data within a
 *                   globally contiguous state vector.
 * @param[in]  rho   Order of contiguous density data within a
 *                   globally contiguous state vector.
 * @param[in]  buf   Working storage of size at least
 *                   <tt>w->n*(w->max_kl + 1 + w->max_ku)</tt>.
 * @param[out] A_T   Storage details for the BSMBSM matrix \c patpt.
 * @param[out] patpt Band storage of renumbered matrix
 *                   \f$PA^{\mbox{T}}P^{\mbox{T}}\f$ which will have
 *                   <tt>A_T->KL</tt> and <tt>A_T->KU</tt> diagonals and
 *                   leading dimension <tt>A_T->LD + A_T->KL</tt>.
 *
 * @see Model documentation in <tt>writeups/derivation.tex</tt> for details.
 * @see suzerain_bsmbsm_zaPxpby() for how to permute state to match \c patpt.
 */
void
suzerain_rholut_imexop_packf(
        const complex_double phi,
        const double km,
        const double kn,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int rho_E,
        const int rho_w,
        const int rho_v,
        const int rho_u,
        const int rho,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double *const patpt);

/**
 * Pack the ``zero-zero'' mode portion of '\f$\left(M +
 * \varphi{}L\right)^{\mbox{T}}\f$ into the corresponding locations within
 * contiguous, LU factorization-ready storage of \f$P\left(M +
 * \varphi{}L\right)^{\mbox{T}}P^{\mbox{T}}\f$ for use by the LAPACK's
 * <tt>gbtrf</tt> or <tt>gbsv</tt>.  This method is equivalent to calling
 * suzerain_rholut_imexop_packf() using <tt>km == 0</tt> and <tt>kn == 0</tt>
 * but should be appreciably faster.
 *
 * @see The documentation of suzerain_rholut_imexop_packf() for details.
 */
void
suzerain_rholut_imexop_packf00(
        const complex_double phi,
        const suzerain_rholut_imexop_scenario * const s,
        const suzerain_rholut_imexop_ref      * const r,
        const suzerain_rholut_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        const int rho_E,
        const int rho_w,
        const int rho_v,
        const int rho_u,
        const int rho,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double *const patpt);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RHOLUT_IMEXOP_H */
