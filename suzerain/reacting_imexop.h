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
 * implicit/explicit" treatment within the model document.  The logic is coded
 * in C99 to facilitate profiling.
 */

#include <suzerain/bsmbsm.h>
#include <suzerain/bsplineop.h>
#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

// FIXME: Non of this scenario info belongs anymore except alpha
/** Scenario-like constants used for operator formation. */
typedef struct {
    double alpha; /**< \f$\alpha\f$ such that \f$\mu_B = \alpha\mu\f$ */
} suzerain_reacting_imexop_scenario;

/** Wall-normal diagonal reference matrices used for operator formation. */
typedef struct {
    double *ux;         /**< \f$ C^{u_x}                \f$ */
    double *uy;         /**< \f$ C^{u_y}                \f$ */
    double *uz;         /**< \f$ C^{u_z}                \f$ */
    double *uxuy;       // TODO: Doxygen
    double *uzuy;
    double *p_ru;
    double *p_rw;
    double *p_rE;
    double *vp_ru;
    double *vp_rw;
    double *vp_rE;
    double *Cmy_rho;
    double *Ce_rho;
    double *Ce_rv;
    double *nu;
    double *korCv;
    double *Ds;
    double *T;
    double *gamma;
    double *a;
    double *p;
    double *p2;
} suzerain_reacting_imexop_ref;

/** Strides between elements in a \ref suzerain_reacting_imexop_ref. */
typedef struct {
    int ux;             /**< \copydoc suzerain_reacting_imexop_ref::ux         */
    int uy;             /**< \copydoc suzerain_reacting_imexop_ref::uy         */
    int uz;             /**< \copydoc suzerain_reacting_imexop_ref::uz         */
    int uxuy;           // TODO: Doxygen
    int uzuy;
    int p_ru;
    int p_rw;
    int p_rE;
    int vp_ru;
    int vp_rw;
    int vp_rE;
    int Cmy_rho;
    int Ce_rho;
    int Ce_rv;
    int nu;
    int korCv;
    int Ds;
    int T;
    int gamma;
    int a;
    int p;
    int p2;
} suzerain_reacting_imexop_refld;

/**
 * Accumulate the linear implicit operator application \f$y
 * \leftarrow{} \left(M + \varphi{}L\right) x + \beta{} y\f$ for the
 * flow equations.  The matrix \f$L\f$ is a function of the provided
 * wavenumbers, scenario parameters, and wall-normal reference
 * quantities.  The problem size and discrete operators are taken from
 * the provided B-spline workspace \c w.  Input and output state
 * variables must be stride one.
 *
 * Providing a \c NULL value for any input and its corresponding output (e.g.
 * \c in_rho and \c out_rho) omits the corresponding part of the operator.
 * This may be useful for omitting part of the operator (e.g. the density
 * equation).
 *
 * @param[in] phi Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in] s  Scenario parameters used to form the operator.
 * @param[in] r  Reference quantities used to form the operator.
 * @param[in] ld Strides between reference quantity values.
 * @param[in] w  B-spline workspace providing discrete operators.
 * @param[in] imagzero  Treat imaginary part of all input data as if
 *                      it were uniformly zero?  That is, assume
                        <tt>Im(in_rho{,u,v,w,E}) == 0</tt>.
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
 */
void
suzerain_reacting_flow_imexop_accumulate(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace      * const w,
        const int imagzero,
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
 * Accumulate the linear implicit operator application \f$y
 * \leftarrow{} \left(M + \varphi{}L\right) x + \beta{} y\f$ for the
 * sth species equation.  The matrix \f$L\f$ is a function of the
 * provided scenario parameters and wall-normal reference quantities.
 * The problem size and discrete operators are taken from the provided
 * B-spline workspace \c w.  Input and output state variables must be
 * stride one.
 *
 * @param[in] phi Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in] s  Scenario parameters used to form the operator.
 * @param[in] r  Reference quantities used to form the operator.
 * @param[in] ld Strides between reference quantity values.
 * @param[in] w  B-spline workspace providing discrete operators.
 * @param[in] imagzero  Treat imaginary part of all input data as if
 *                      it were uniformly zero?  That is, assume
                        <tt>Im(in_rho{,u,v,w,E}) == 0</tt>.
 * @param[in]  in_rho_s  Wall-normal input data for \f$\rho{}u\f$.
 * @param[in]  beta      Accumulation scaling factor \f$\beta\f$.
 * @param[out] out_rho_s Wall-normal output data for \f$\rho{}u\f$.
 *
 */
void
suzerain_reacting_species_imexop_accumulate(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace      * const w,
        const int imagzero,
        const complex_double *in_rho_s,
        const complex_double beta,
        complex_double *out_rho_s );

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
suzerain_reacting_flow_imexop_packc(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace      * const w,
        const int rho_E,
        const int rho_w,
        const int rho_v,
        const int rho_u,
        const int rho,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double * const patpt);

/**
 * Pack \f$\left(M + \varphi{}L\right)^{\mbox{T}}\f$ for the species
 * equations into the corresponding locations within contiguous
 * storage of \f$P\left(M +
 * \varphi{}L\right)^{\mbox{T}}P^{\mbox{T}}\f$ for use by the BLAS'
 * <tt>gbmv</tt> or LAPACK's <tt>gbsvx</tt>.  The matrix \f$L\f$ is a
 * function of the provided wavenumbers, scenario parameters, and
 * wall-normal reference quantities.  Notice that the \e transpose of
 * the operator is assembled.
 *
 * The problem size and discrete operators are taken from the provided B-spline
 * workspace \c w.  On entry, \c patpt must be of at least size
 * <tt>(nneg*w->n)*(nneg*(w->max_kl + w->max_ku + 2) - 1)</tt> where \c nneg
 * is the number of nonnegative state order parameters.  Boundary conditions,
 * which are \em not applied, will require using information about the
 * permutation returned in \c A_T.
 *
 * @param[in]  phi   Factor \f$\varphi\f$ used in forming \f$M+\varphi{}L\f$.
 * @param[in]  s     Scenario parameters used to form the operator.
 * @param[in]  r     Reference quantities used to form the operator.
 * @param[in]  ld    Strides between reference quantity values.
 * @param[in]  w     B-spline workspace providing discrete operators.
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
suzerain_reacting_species_imexop_packc(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace      * const w,
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
suzerain_reacting_flow_imexop_packf(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
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
suzerain_reacting_species_imexop_packf(
        const complex_double phi,
        const suzerain_reacting_imexop_scenario * const s,
        const suzerain_reacting_imexop_ref      * const r,
        const suzerain_reacting_imexop_refld    * const ld,
        const suzerain_bsplineop_workspace    * const w,
        complex_double * const buf,
        suzerain_bsmbsm * const A_T,
        complex_double *const patpt);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_RHOLUT_IMEXOP_H */
