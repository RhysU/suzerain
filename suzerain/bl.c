//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

/** @file
 * @copydoc bl.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bl.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sys.h>

#include <suzerain/common.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/pre_gsl.h>

static inline double square(double x) { return x*x; }

// Defensively NaN assuming p points to a type containing only doubles
#define FILL_WITH_NANS(p)                                    \
    for (size_t i = 0; i < sizeof(*p)/sizeof(double); ++i) { \
        ((double *) p)[i] = GSL_NAN;                         \
    }

int
suzerain_bl_compute_viscous(
        const double code_Re,
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous)
{
    FILL_WITH_NANS(viscous);

    // Compute quantities of interest in "code units" divided by  [units    ]
    viscous->tau_w    = wall->mu * wall->u__y;                 // [\mu u / l]
    viscous->u_tau    = sqrt(viscous->tau_w / wall->rho);      // [u        ]
    viscous->delta_nu = wall->mu / wall->rho / viscous->u_tau; // [l        ]

    // Adjust for code Reynolds number per writeups/viscous_nondim.pdf which
    // details why and how how these scaling factors appear.  Provided that
    // they are addressed here, they don't "leak" into other computations.
    const double sqrt_code_Re = sqrt(code_Re);
    viscous->tau_w    /= code_Re;
    viscous->u_tau    /= sqrt_code_Re;
    viscous->delta_nu /= sqrt_code_Re;

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_find_edge(
    const double * coeffs_H0,
    const double lowerbnd,
    const double upperbnd,
    double * location,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Everything hinges on the second derivative of H0 crossing 0.0 */
    enum { nderiv = 2, threshold = 0 };

    /* Ensure the basis has at least non-constant derivatives of interest */
    if (SUZERAIN_UNLIKELY(w->k < nderiv + 2)) {
        SUZERAIN_ERROR("Basis must have non-constant second derivatives",
                       SUZERAIN_EINVAL);
    }

    /* Assume failure until proven otherwise. */
    *location  = GSL_NAN;
    int status = SUZERAIN_EFAILED;

    /* Searching for bounding intervals within breakpoints.                  */
    /* gsl_interp_bsearch has exactly the semantics that we want, so use it. */
    assert(w->knots->stride == 1);
    const double * const breakpts = w->knots->data+w->k-1; // gsl_bspline_nbreak
    const size_t ilo = gsl_interp_bsearch(breakpts, lowerbnd, 0,   w->nbreak);
    const size_t ihi = gsl_interp_bsearch(breakpts, upperbnd, ilo, w->nbreak);

    /* Start by evaluating function on the trimmed lower bound... */
    double flower = GSL_NAN, lower = GSL_MAX(lowerbnd, breakpts[ilo]);
    suzerain_bspline_linear_combination(
            nderiv, coeffs_H0, 1, &lower, &flower, 0, dB, w, dw);

    /* ...look breakpoint-by-breakpoint for upwards crossing of threshold... */
    double fupper = GSL_NAN, upper = GSL_NAN;
    for (size_t i = ilo+1; i <= ihi; ++i, lower = upper, flower = fupper) {

        /* ...(honoring any requested upper bound)... */
        upper = breakpts[i];
        if (upper > upperbnd) break;
        suzerain_bspline_linear_combination(
                nderiv, coeffs_H0, 1, &upper, &fupper, 0, dB, w, dw);

        /* ...when found, polish the crossing bracket into a location... */
        if (flower < threshold && fupper >= threshold) {
            const double epsabs = 10*GSL_DBL_EPSILON;
            status = suzerain_bspline_crossing(
                    nderiv, coeffs_H0, threshold, &lower, &upper, 100,
                    epsabs, /*not relative*/ 0, location, dB, w, dw);
            break;
        }

    }
    /* ...if never found, location remains NAN and status reflects failure. */

    return status;
}

int
suzerain_bl_displacement_thickness(
    const double * coeffs_rhou,
    double * delta1,
    gsl_vector *Bk,
    gsl_bspline_workspace *w,
    const gsl_integration_glfixed_table * tbl)
{
    /* Integrand has same piecewise polynomial order as basis. */
    if (SUZERAIN_UNLIKELY(tbl->n < (w->k + 1)/2)) {
        *delta1 = GSL_NAN;
        SUZERAIN_ERROR("Gaussian quadrature table order too low",
                       SUZERAIN_ESANITY);
    }

    /* Momentum at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double edge_rhou = coeffs_rhou[w->n - 1];

    /* Accumulate the breakpoint-by-breakpoint contributions into *delta1 */
    *delta1 = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_eval_nonzero(xj, Bk, &kstart, &kend, w);

            /* Accumulate basis linear combinations to evaluate rhou */
            double rhou = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                rhou += coeffs_rhou[k] * gsl_vector_get(Bk, k - kstart);
            }

            /* Then use integrand scaled by weight to accumulate result */
            *delta1 += wj * (1 - rhou / edge_rhou);
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_momentum_thickness(
    const double * coeffs_rhou,
    const double * coeffs_u,
    double * delta2,
    gsl_vector *Bk,
    gsl_bspline_workspace *w,
    const gsl_integration_glfixed_table * tbl)
{
    /* Integrand has twice the piecewise polynomial order of the basis. */
    /* That is, solve 2*(k - 1) = 2*n - 1 for number of Gauss points n. */
    if (SUZERAIN_UNLIKELY(tbl->n < w->k)) {
        *delta2 = GSL_NAN;
        SUZERAIN_ERROR("Gaussian quadrature table order too low",
                       SUZERAIN_ESANITY);
    }

    /* State at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double edge_rhou = coeffs_rhou[w->n - 1];
    const double edge_u    = coeffs_u   [w->n - 1];

    /* Accumulate the breakpoint-by-breakpoint contributions into *delta2 */
    *delta2 = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_eval_nonzero(xj, Bk, &kstart, &kend, w);

            /* Accumulate basis linear combinations to evaluate rhou, u */
            double rhou = 0;
            double u    = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                const double B = gsl_vector_get(Bk, k - kstart);
                rhou += coeffs_rhou[k] * B;
                u    += coeffs_u[k]    * B;
            }

            /* Then use integrand scaled by weight to accumulate result */
            *delta2 += wj * (rhou / edge_rhou) * (1 - u / edge_u);
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_energy_thickness(
    const double * coeffs_ke,
    const double * coeffs_rhou,
    double * delta3,
    gsl_vector *Bk,
    gsl_bspline_workspace *w,
    const gsl_integration_glfixed_table * tbl)
{
    // The form of the energy thickness definition is nothing but the momentum
    // thickness with ke = u^2/2 replacing u.  Notice argument order changed.
    return suzerain_bl_momentum_thickness(
            coeffs_rhou, coeffs_ke, delta3, Bk, w, tbl);
}

int
suzerain_bl_enthalpy_thickness(
    const double code_Ma,
    const double * coeffs_H0,
    const double * coeffs_ke,
    const double * coeffs_rhou,
    double * deltah,
    gsl_vector * Bk,
    gsl_bspline_workspace * w,
    const gsl_integration_glfixed_table * tbl)
{
    /* Integrand has twice the piecewise polynomial order of the basis. */
    /* That is, solve 2*(k - 1) = 2*n - 1 for number of Gauss points n. */
    if (SUZERAIN_UNLIKELY(tbl->n < w->k)) {
        *deltah = GSL_NAN;
        SUZERAIN_ERROR("Gaussian quadrature table order too low",
                       SUZERAIN_ESANITY);
    }

    /* State at wall taken from first B_spline collocation point */
    /* which happens to be the value of the first coefficient */
    const double wall_H0 = coeffs_H0[0];
    const double wall_ke = coeffs_ke[0];
    const double wall_h  = wall_H0 - code_Ma*code_Ma*wall_ke;

    /* State at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double edge_H0   = coeffs_H0  [w->n - 1];
    const double edge_ke   = coeffs_ke  [w->n - 1];
    const double edge_rhou = coeffs_rhou[w->n - 1];
    const double edge_h    = edge_H0 - code_Ma*code_Ma*edge_ke;

    /* Accumulate the breakpoint-by-breakpoint contributions into *deltah */
    *deltah = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_eval_nonzero(xj, Bk, &kstart, &kend, w);

            /* Accumulate basis linear combinations to evaluate rhou, u */
            double H0   = 0;
            double ke   = 0;
            double rhou = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                const double B = gsl_vector_get(Bk, k - kstart);
                H0   += coeffs_H0  [k] * B;
                ke   += coeffs_ke  [k] * B;
                rhou += coeffs_rhou[k] * B;
            }

            /* Then use integrand scaled by weight to accumulate result */
            const double h = H0 - code_Ma*code_Ma*ke;
            *deltah += wj * (rhou/edge_rhou) * (edge_h - h)/(edge_h - wall_h);
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_compute_thicknesses(
    const double code_Ma,
    const double * coeffs_H0,
    const double * coeffs_ke,
    const double * coeffs_rhou,
    const double * coeffs_u,
    suzerain_bl_thicknesses * thick,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw)
{
    FILL_WITH_NANS(thick);

    // Prepare repeatedly used resources
    int status                         = SUZERAIN_SUCCESS;
    gsl_matrix *dB                     = NULL;
    gsl_integration_glfixed_table *tbl = NULL;
    if (NULL == (dB = gsl_matrix_alloc(w->k, 3))) {
        SUZERAIN_ERROR_REPORT("failed to allocate dB",
                              (status = SUZERAIN_ENOMEM));
    }
    gsl_vector_view Bk = gsl_matrix_column(dB, 0);
    if (NULL == (tbl = gsl_integration_glfixed_table_alloc(w->k))) {
        SUZERAIN_ERROR_REPORT("failed to allocate tbl",
                              (status = SUZERAIN_ENOMEM));
    }
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_displacement_thickness(
            coeffs_rhou, &thick->delta1, &Bk.vector, w, tbl);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    // As \delta_{99} should always be outside \delta_1 and the latter
    // is more robust than the former, use \delta_1 as a lower bound on
    // where we might find \delta_{99}.
    status = suzerain_bl_find_edge(
            coeffs_H0, thick->delta1, gsl_bspline_breakpoint(w->nbreak-1, w),
            &thick->delta, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_momentum_thickness(
            coeffs_rhou, coeffs_u, &thick->delta2, &Bk.vector, w, tbl);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_energy_thickness(
            coeffs_ke, coeffs_rhou, &thick->delta3, &Bk.vector, w, tbl);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_enthalpy_thickness(
            code_Ma, coeffs_H0, coeffs_ke, coeffs_rhou,
            &thick->deltah, &Bk.vector, w, tbl);
    /* Done regardless of status */

done:

    gsl_integration_glfixed_table_free(tbl);
    gsl_matrix_free(dB);
    return status;
}

int
suzerain_bl_compute_reynolds(
    const double code_Re,
    const suzerain_bl_local       * const edge,
    const suzerain_bl_thicknesses * const thick,
          suzerain_bl_reynolds    * const reynolds)
{
    FILL_WITH_NANS(reynolds);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the final line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_bl_viscous struct declaration to check their scaling.
    reynolds->delta  = edge->rho * edge->u * thick->delta  / edge->mu
                     * code_Re;
    reynolds->delta1 = edge->rho * edge->u * thick->delta1 / edge->mu
                     * code_Re;
    reynolds->delta2 = edge->rho * edge->u * thick->delta2 / edge->mu
                     * code_Re;
    reynolds->delta3 = edge->rho * edge->u * thick->delta3 / edge->mu
                     * code_Re;
    reynolds->deltah = edge->rho * edge->u * thick->deltah / edge->mu
                     * code_Re;

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_compute_qoi(
        const double code_Ma,
        const double code_Re,
        const suzerain_bl_local       * const wall,
        const suzerain_bl_viscous     * const viscous,
        const suzerain_bl_local       * const edge,
        const suzerain_bl_thicknesses * const thick,
              suzerain_bl_qoi         * const qoi)
{
    FILL_WITH_NANS(qoi);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the final line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_bl_viscous struct declaration to check their scaling.
    qoi->cf          = 2 * viscous->tau_w / edge->rho / square(edge->u)
                     * 1;
    qoi->gamma_e     = edge->gamma
                     * 1;
    qoi->Ma_e        = edge->u / edge->a
                     * code_Ma;
    qoi->Ma_tau      = viscous->u_tau / wall->a
                     * code_Ma;
    qoi->Pr_w        = wall->Pr;
    qoi->Bq          = - (wall->mu * wall->T__y)
                     / (qoi->Pr_w * wall->rho * viscous->u_tau * wall->T)
                     / code_Re;
    qoi->ratio_nu    = (edge->mu / edge->rho) / (wall->mu / wall->rho)
                     * 1;
    qoi->ratio_rho   = edge->rho / wall->rho
                     * 1;
    qoi->ratio_T     = edge->T / wall->T
                     * 1;
    qoi->shapefactor = thick->delta1 / thick->delta2
                     * 1;
    qoi->v_wallplus  = wall->v / viscous->u_tau
                     * 1;

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_compute_pg(
        const double code_Ma,
        const double code_Re,
        const suzerain_bl_local       * const wall,
        const suzerain_bl_viscous     * const viscous,
        const suzerain_bl_local       * const edge,
        const double                          edge_p__x,
        const double                          edge_u__x,
        const suzerain_bl_thicknesses * const thick,
              suzerain_bl_pg          * const pg)
{
    FILL_WITH_NANS(pg);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the second line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_bl_viscous struct declaration to check their scaling.
    pg->Clauser      = thick->delta1 / viscous->tau_w * edge_p__x
                     / square(code_Ma);
    pg->Launder_e    = edge->mu * edge_u__x / edge->rho / square(edge->u)
                     / code_Re;
    pg->Launder_w    = wall->mu * edge_u__x / edge->rho / square(edge->u)
                     / code_Re;
    pg->Lambda_n     = - thick->delta / viscous->tau_w * edge_p__x
                     / square(code_Ma);
    pg->p_ex         = thick->delta / edge->rho / square(edge->u) * edge_p__x
                     / square(code_Ma);
    pg->Pohlhausen   = square(thick->delta) * edge->rho / edge->mu * edge_u__x
                     * code_Re;

    return SUZERAIN_SUCCESS;
}

// Parameters for baseflow-ready integrand_thickness_displacement
typedef struct {
    double                  inner_cutoff;  // Must be first member
    const double          * vis_rhou;      // Always stride one
    int                     inv_stride;    // Stride for inv_rhou
    const double          * inv_rhou;      // May have stride zero
    gsl_vector            * Bk;
    gsl_bspline_workspace * w;
} params_thickness_displacement;

// Per writeups/thicknesses.pdf, when p->inner_cutoff == \delta^\ast this
// integral evaluated over [0, \infty) should be zero.
static
double integrand_thickness_displacement(
        double y, void * params)
{
    const params_thickness_displacement * const p
            = (params_thickness_displacement *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {
        integrand = 0.0;
        const int inv_stride = p->inv_stride;
        if (y >= p->inner_cutoff) {
            for (size_t i = istart; i <= iend; ++i) {
                integrand += (p->inv_rhou[i*inv_stride] - p->vis_rhou[i])
                           * gsl_vector_get(p->Bk, i-istart);
            }
        } else {
            for (size_t i = istart; i <= iend; ++i) {
                integrand -= p->vis_rhou[i]
                           * gsl_vector_get(p->Bk, i-istart);
            }
        }
    }
    return integrand;
}

// Parameters for baseflow-ready integrand_reynolds_displacement
typedef struct {
    const double          *vis_rhou;     // Always stride one
    int                    inv_stride;   // Stride for inv_rhou
    const double          *inv_rhou;     // May have stride zero
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_reynolds_displacement;

// Per writeups/thicknesses.pdf, this integral evaluated over
// [0, \infty) should produce \mu * Re_{\delta^\ast}
static
double integrand_reynolds_displacement(
        double y, void * params)
{
    const params_reynolds_displacement * const p
            = (params_reynolds_displacement *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {
        integrand = 0.0;
        const int inv_stride = p->inv_stride;
        for (size_t i = istart; i <= iend; ++i) {
            integrand += (p->inv_rhou[i*inv_stride] - p->vis_rhou[i])
                       * gsl_vector_get(p->Bk, i-istart);
        }
    }
    return integrand;
}

// Parameters for baseflow-ready integrand_thickness_momentum
typedef struct {
    double                 inner_cutoff;  // Must be first member
    const double          *vis_rhou;      // Always stride one
    const double          *vis_u;         // Always stride one
    int                    inv_stride;    // Stride for inv_rhou, inv_u
    const double          *inv_rhou;      // Strided per inv_stride
    const double          *inv_u;         // Strided per inv_stride
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_thickness_momentum;

// Per writeups/thicknesses.pdf, when p->inner_cutoff == \delta^\ast + \theta
// this integral evaluated over [0, \infty) should be zero.
static
double integrand_thickness_momentum(
        double y, void * params)
{
    const params_thickness_momentum * const p
            = (params_thickness_momentum *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_u = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_u += p->vis_u[i] * gsl_vector_get(p->Bk, i-istart);
        }

        integrand = - vis_rhou * vis_u;

        if (y >= p->inner_cutoff) {
            const int inv_stride = p->inv_stride;

            double inv_rhou = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_rhou += p->inv_rhou[i*inv_stride]
                          * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_u = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_u += p->inv_u[i*inv_stride]
                       * gsl_vector_get(p->Bk, i-istart);
            }

            integrand += inv_rhou * inv_u;

        }

    }
    return integrand;
}

// Parameters for baseflow-ready integrand_reynolds_momentum
typedef struct {
    const double          *vis_rhou;    // Always stride one
    const double          *vis_u;       // Always stride one
    int                    inv_stride;  // Stride for inv_u
    const double          *inv_u;       // Strided per inv_stride
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_reynolds_momentum;

// Per writeups/thicknesses.pdf, this integral evaluated over
// [0, \infty) should produce \mu * Re_{\theta}
static
double integrand_reynolds_momentum(
        double y, void * params)
{
    const params_reynolds_momentum * const p
            = (params_reynolds_momentum *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_u = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_u += p->vis_u[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_u = 0;
        const int inv_stride = p->inv_stride;
        for (size_t i = istart; i <= iend; ++i) {
            inv_u += p->inv_u[i*inv_stride] * gsl_vector_get(p->Bk, i-istart);
        }

        integrand = 1 - vis_u / inv_u;

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        integrand *= vis_rhou;

    }
    return integrand;
}

// Parameters for baseflow-ready integrand_thickness_energy
typedef struct {
    double                 inner_cutoff;  // Must be first member
    const double          *vis_ke;        // Always stride one
    const double          *vis_rhou;      // Always stride one
    int                    inv_stride;    // Stride for inv_rhou, inv_u, inv_v
    const double          *inv_rhou;      // May have stride zero
    const double          *inv_u;         // May have stride zero
    const double          *inv_v;         // May have stride zero
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_thickness_energy;

// Per writeups/thicknesses.pdf, when p->inner_cutoff == \delta^\ast + \delta_3
// this integral evaluated over [0, \infty) should be zero.
static
double integrand_thickness_energy(
        double y, void * params)
{
    const params_thickness_energy * const p
            = (params_thickness_energy *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_ke = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_ke += p->vis_ke[i] * gsl_vector_get(p->Bk, i-istart);
        }

        integrand = - vis_rhou * (2*vis_ke);

        if (y >= p->inner_cutoff) {
            const int inv_stride = p->inv_stride;

            double inv_rhou = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_rhou += p->inv_rhou[i*inv_stride]
                          * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_u = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_u += p->inv_u[i*inv_stride]
                       * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_v = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_v += p->inv_v[i*inv_stride]
                       * gsl_vector_get(p->Bk, i-istart);
            }

            integrand += inv_rhou * (inv_u*inv_u + inv_v*inv_v);

        }

    }
    return integrand;
}

// Parameters for baseflow-ready integrand_reynolds_energy
typedef struct {
    const double          *vis_ke;       // Always stride one
    const double          *vis_rhou;     // Always stride one
    int                    inv_stride;   // Stride for inv_u
    const double          *inv_u;        // Strided per inv_stride
    const double          *inv_v;        // Strided per inv_stride
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_reynolds_energy;

// Per writeups/thicknesses.pdf, this integral evaluated over
// [0, \infty) should produce \mu * Re_{\delta_3}
static
double integrand_reynolds_energy(
        double y, void * params)
{
    const params_reynolds_energy * const p
            = (params_reynolds_energy *) params;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_ke = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_ke += p->vis_ke[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_u = 0;
        for (size_t i = istart; i <= iend; ++i) {
            inv_u += p->inv_u[i*p->inv_stride]
                   * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_v = 0;
        for (size_t i = istart; i <= iend; ++i) {
            inv_v += p->inv_v[i*p->inv_stride]
                   * gsl_vector_get(p->Bk, i-istart);
        }

        integrand = vis_rhou
                  * (1 - (2*vis_ke) / (inv_u*inv_u + inv_v*inv_v));

    }
    return integrand;
}

// Parameters for baseflow-ready integrand_thickness_enthalpy
typedef struct {
    double                 inner_cutoff;  // Must be first member
    double                 code_Ma;       // Scales kinetic energies
    const double          *vis_H0;        // Always stride one
    const double          *vis_ke;        // Always stride one
    const double          *vis_rhou;      // Always stride one
    int                    inv_stride;    // Stride for inv_rhou, inv_h
    const double          *inv_H0;        // Strided per inv_stride
    const double          *inv_rhou;      // Strided per inv_stride
    const double          *inv_u;         // Strided per inv_stride
    const double          *inv_v;         // Strided per inv_stride
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_thickness_enthalpy;

// Per writeups/thicknesses.pdf, when p->inner_cutoff == \delta^\ast + \delta_h
// this integral evaluated over [0, \infty) should be zero.
static
double integrand_thickness_enthalpy(
        double y, void * params)
{
    const params_thickness_enthalpy * const p
            = (params_thickness_enthalpy *) params;
    const double Ma2 = p->code_Ma*p->code_Ma;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_H0 = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_H0 += p->vis_H0[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_ke = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_ke += p->vis_ke[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        integrand = - vis_rhou * (vis_H0 - Ma2*vis_ke);

        if (y >= p->inner_cutoff) {
            const int inv_stride = p->inv_stride;

            double inv_H0 = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_H0 += p->inv_H0[i*inv_stride]
                        * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_rhou = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_rhou += p->inv_rhou[i*inv_stride]
                          * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_u = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_u += p->inv_u[i*inv_stride]
                       * gsl_vector_get(p->Bk, i-istart);
            }

            double inv_v = 0;
            for (size_t i = istart; i <= iend; ++i) {
                inv_v += p->inv_v[i*inv_stride]
                       * gsl_vector_get(p->Bk, i-istart);
            }

            const double inv_h = inv_H0 - Ma2*(inv_u*inv_u + inv_v*inv_v)/2;
            integrand += inv_rhou * inv_h;

        }

    }
    return integrand;
}

// Parameters for baseflow-ready integrand_reynolds_enthalpy
typedef struct {
    double                 code_Ma;     // Scales kinetic energies
    const double          *vis_H0;      // Always stride one
    const double          *vis_ke;      // Always stride one
    const double          *vis_rhou;    // Always stride one
    int                    inv_stride;  // Stride for inv_rhou
    const double          *inv_H0;      // Strided per inv_stride
    const double          *inv_u;       // Strided per inv_stride
    const double          *inv_v;       // Strided per inv_stride
    gsl_vector            *Bk;
    gsl_bspline_workspace *w;
} params_reynolds_enthalpy;

// Per writeups/thicknesses.pdf, this integral evaluated over
// [0, \infty) should produce \mu * Re_{\delta_h}
static
double integrand_reynolds_enthalpy(
        double y, void * params)
{
    const params_reynolds_enthalpy * const p
            = (params_reynolds_enthalpy *) params;
    const double Ma2 = p->code_Ma*p->code_Ma;
    size_t istart, iend;
    double integrand = GSL_NAN;
    if (!gsl_bspline_eval_nonzero(y, p->Bk, &istart, &iend, p->w)) {

        double vis_H0 = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_H0 += p->vis_H0[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_ke = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_ke += p->vis_ke[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double vis_rhou = 0;
        for (size_t i = istart; i <= iend; ++i) {
            vis_rhou += p->vis_rhou[i] * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_H0 = 0;
        for (size_t i = istart; i <= iend; ++i) {
            inv_H0 += p->inv_H0[i*p->inv_stride]
                    * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_u = 0;
        for (size_t i = istart; i <= iend; ++i) {
            inv_u += p->inv_u[i*p->inv_stride]
                    * gsl_vector_get(p->Bk, i-istart);
        }

        double inv_v = 0;
        for (size_t i = istart; i <= iend; ++i) {
            inv_v += p->inv_v[i*p->inv_stride]
                    * gsl_vector_get(p->Bk, i-istart);
        }

        const double vis_h = vis_H0 - Ma2*vis_ke;
        const double inv_h = inv_H0 - Ma2*(inv_u*inv_u + inv_v*inv_v)/2;
        integrand = vis_rhou*(1 - vis_h / inv_h);

    }
    return integrand;
}

int
suzerain_bl_compute_reynolds_baseflow(
    const double                    code_Ma,
    const double                    code_Re,
    const double            * const coeffs_vis_H0,
    const double            * const coeffs_vis_ke,
    const double            * const coeffs_vis_rhou,
    const double            * const coeffs_vis_u,
    const int                       inv_stride,
    const double            * const coeffs_inv_H0,
    const double            * const coeffs_inv_rhou,
    const double            * const coeffs_inv_u,
    const double            * const coeffs_inv_v,
    const suzerain_bl_local * const edge,
    suzerain_bl_reynolds    * const reynolds,
    gsl_bspline_workspace   * const w)
{
    FILL_WITH_NANS(reynolds);

    // Tracks status of first failure for return from routine.
    // Where sensible, processing continues on a "best effort" basis.
    int status = SUZERAIN_SUCCESS;

    // See suzerain_bl_compute_reynolds for the constant baseflow version of
    // these quantities, including commentary on the code unit scaling.  See
    // writeups/thicknesses.pdf for how the integrals in this function are
    // defined and used.

    // Prepare repeatedly used inputs to gsl_integration_qagp(...)
    assert(w->knots->stride == 1);                // Breakpoints are contiguous
    double * const pts  = w->knots->data+w->k-1;  // Per gsl_bspline_nbreak
    const size_t npts   = w->nbreak;              // "Singular" at breakpoints
    const size_t limit  = npts * (1+w->k/2);      // Hopefully sufficient for..
    const double epsabs = GSL_DBL_EPSILON;        // ..either this...
    const double epsrel = GSL_SQRT_DBL_EPSILON;   // ..or this tolerance
    double abserr       = GSL_NAN;                // Repeatedly used in calls

    gsl_integration_workspace * iw = NULL;
    gsl_vector                * Bk = NULL;
    if (NULL == (iw = gsl_integration_workspace_alloc(limit))) {
        SUZERAIN_ERROR_REPORT("failed to allocate iw",
                              (status = SUZERAIN_ENOMEM));
    }
    if (NULL == (Bk = gsl_vector_alloc(w->k))) {
        SUZERAIN_ERROR_REPORT("failed to allocate Bk",
                              (status = SUZERAIN_ENOMEM));
    }
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    // Re_delta from edge->rho, edge->u, thick->delta, edge->mu, and code_Re
    {
        reynolds->delta = edge->rho * edge->u * edge->y / edge->mu
                        * code_Re;
    }

    // Nondimensional quantities are computed with the first line being the
    // quantity and the final line being any needed "code unit" correction.

    // Re_delta1 from integrand_reynolds_displacement, edge->mu, code_Re
    {
        params_reynolds_displacement params = {
            .vis_rhou   = coeffs_vis_rhou,
            .inv_stride = inv_stride,
            .inv_rhou   = coeffs_inv_rhou,
            .Bk         = Bk,
            .w          = w
        };
        gsl_function F = {
            &integrand_reynolds_displacement,
            &params
        };
        int tmp = gsl_integration_qagp(&F, pts, npts, epsabs, epsrel,
                                       limit, iw, &reynolds->delta1, &abserr);
        reynolds->delta1 *= code_Re / edge->mu;
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

    // Re_delta2 from integrand_reynolds_momentum, edge->mu, code_Re
    {
        params_reynolds_momentum params = {
            .vis_rhou    = coeffs_vis_rhou,
            .vis_u       = coeffs_vis_u,
            .inv_stride  = inv_stride,
            .inv_u       = coeffs_inv_u,
            .Bk          = Bk,
            .w           = w
        };
        gsl_function F = {
            .function = &integrand_reynolds_momentum,
            .params   = &params
        };
        int tmp = gsl_integration_qagp(&F, pts, npts, epsabs, epsrel,
                                       limit, iw, &reynolds->delta2, &abserr);
        reynolds->delta2 *= code_Re / edge->mu;
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

    // Re_delta3 from integrand_reynolds_energy, edge->mu, code_Re
    {
        params_reynolds_energy params = {
            .vis_ke      = coeffs_vis_ke,
            .vis_rhou    = coeffs_vis_rhou,
            .inv_stride  = inv_stride,
            .inv_u       = coeffs_inv_u,
            .inv_v       = coeffs_inv_v,
            .Bk          = Bk,
            .w           = w
        };
        gsl_function F = {
            .function = &integrand_reynolds_energy,
            .params   = &params
        };
        int tmp = gsl_integration_qagp(&F, pts, npts, epsabs, epsrel,
                                       limit, iw, &reynolds->delta3, &abserr);
        reynolds->delta3 *= code_Re / edge->mu;
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

    // Re_delta3 from integrand_reynolds_enthalpy, edge->mu, code_Re
    {
        params_reynolds_enthalpy params = {
            .code_Ma     = code_Ma,
            .vis_H0      = coeffs_vis_H0,
            .vis_ke      = coeffs_vis_ke,
            .vis_rhou    = coeffs_vis_rhou,
            .inv_stride  = inv_stride,
            .inv_H0      = coeffs_inv_H0,
            .inv_u       = coeffs_inv_u,
            .inv_v       = coeffs_inv_v,
            .Bk          = Bk,
            .w           = w
        };
        gsl_function F = {
            .function = &integrand_reynolds_enthalpy,
            .params   = &params
        };
        int tmp = gsl_integration_qagp(&F, pts, npts, epsabs, epsrel,
                                       limit, iw, &reynolds->deltah, &abserr);
        reynolds->deltah *= code_Re / edge->mu;
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

done:
    gsl_integration_workspace_free(iw);
    gsl_vector_free(Bk);

    return status;
}

// Parameters for baseflow-ready integral_thickness_residual requiring a
// function for root-finding wrapping a function for the integration integrand
typedef struct params_integral_thickness_residual {
    gsl_function   integrand;  // Using, e.g. integrand_thickness_displacement
    size_t         ndis;       // # of discontinuities (breakpoints) in basis
    double       * dis;        // Increasing list of discontinuities
    size_t         npts;       // Buffer size for QAGP discontinuities >= ndis+1
    double       * pts;        // Buffer for QAGP discontinuities
    double         epsabs;     // Request relative...
    double         epsrel;     // ...or absolute tolerance to satistfy.
    int            status;     // Stores abserr from gsl_integration_qagp
    double         abserr;     // Stores return value from gsl_integration_qagp
    size_t         limit;      // Maximum number of subdivisions to perform
    gsl_integration_workspace * iw;
} params_integral_thickness_residual;

// Per writeups/thicknesses.pdf, evaluate the residual of an implicit thickness
// equation (e.g. \delta^\ast as defined by equation (*) on page 2) given that
// thickness.
static
double integral_thickness_residual(
        const double thick,
        void * params)
{
    params_integral_thickness_residual * const p =
        (params_integral_thickness_residual *) params;

    // Statically assert that inner_cutoff is the first member for
    // the variety of integrand.params types that may be encountered...
    enum {
        assert1 = 1/(0==offsetof(params_thickness_displacement, inner_cutoff)),
        assert2 = 1/(0==offsetof(params_thickness_energy,       inner_cutoff)),
        assert3 = 1/(0==offsetof(params_thickness_enthalpy,     inner_cutoff)),
        assert4 = 1/(0==offsetof(params_thickness_momentum,     inner_cutoff))
    };
    *((double *) p->integrand.params) = thick; // ...so it can be set thusly.

    // Insert thick into a copied list of discontinuities/singularities
    {
        assert(p->npts >= p->ndis + 1);
        size_t i = 0;
        for (; i < p->ndis && p->dis[i] <= thick; ++i)
            p->pts[i] = p->dis[i];
        p->pts[i++] = thick;
        for (; i < p->ndis; ++i)
            p->pts[i] = p->dis[i-1];
    }

    // Uniquely insert thick into a copied list of discontinuities:
    //   (a) When thick is on a discontinuity, we have ndis points.
    //   (b) When thick is not, we have ndis+1 points.
    // Be careful to distinguish the two.
    assert(p->npts >= p->ndis + 1);
    size_t i = 0;
    for (; i < p->ndis && p->dis[i] < thick; ++i)
        p->pts[i] = p->dis[i];
    size_t j = i;
#pragma warning(push,disable:1572)
    if (i < p->ndis && thick != p->dis[i])
        p->pts[i++] = thick;
#pragma warning(pop)
    for (; j < p->ndis; ++i, ++j)
        p->pts[i] = p->dis[j];
    assert(i == p->ndis || i == p->ndis + 1);

    // Adaptively evaluate the integrand given p->pts[0], ..., p->pts[i-1]
    double result = GSL_NAN;
    p->status = gsl_integration_qagp(&p->integrand,
                                     p->pts, i,
                                     p->epsabs, p->epsrel,
                                     p->limit, p->iw,
                                     &result, &p->abserr);
    return result;
}

static
int fsolver_solve(
        gsl_root_fsolver * const s,
        gsl_function     * const f,
        const size_t             maxiter,
        const double             epsabs,
        const double             epsrel,
        double           * const lower,
        double           * const upper,
        double           * const root)
{
    // Initialize fsolver on [lower, upper]
    gsl_error_handler_t * h = gsl_set_error_handler_off();    // Push handler
    int status = gsl_root_fsolver_set(s, f, *lower, *upper);
    status = (GSL_SUCCESS == status) ? GSL_CONTINUE : status;

    // Proceed until success, failure, or maximum iterations reached
    // On success, overwrite *location as described in API documentation.
    *root = GSL_NAN;
    for (size_t iter = 1; iter <= maxiter && status == GSL_CONTINUE; ++iter) {
        if (GSL_SUCCESS != (status = gsl_root_fsolver_iterate(s))) break;
        *lower = gsl_root_fsolver_x_lower(s);
        *upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(*lower, *upper, epsabs, epsrel);
    }
    if (status == GSL_SUCCESS) {
        *root = gsl_root_fsolver_root(s);
    }
    gsl_set_error_handler(h);                                 // Pop handler
    return status;
}

int
suzerain_bl_compute_thicknesses_baseflow(
    const double                        code_Ma,
    const double                * const coeffs_vis_H0,
    const double                * const coeffs_vis_ke,
    const double                * const coeffs_vis_rhou,
    const double                * const coeffs_vis_u,
    const int                           inv_stride,
    const double                * const coeffs_inv_H0,
    const double                * const coeffs_inv_rhou,
    const double                * const coeffs_inv_u,
    const double                * const coeffs_inv_v,
    suzerain_bl_thicknesses     * const thick,
    gsl_bspline_workspace       * const w,
    gsl_bspline_deriv_workspace * const dw)
{
    FILL_WITH_NANS(thick);

    // Tracks status of first failure for return from routine.
    // Where sensible, processing continues on a "best effort" basis.
    int status = SUZERAIN_SUCCESS;

    // Prepare inputs and allocate resources for subroutine calls
    enum { maxiter = 255 };                     // How long can we wait?
    params_integral_thickness_residual params = { .pts = NULL, .iw = NULL };
    assert(w->knots->stride == 1);              // Breakpoints are contiguous
    params.ndis   = w->nbreak;                  // Smoothness drops @ breakpoint
    params.dis    = w->knots->data + w->k - 1;  // Per gsl_bspline_nbreak
    params.npts   = params.ndis + 1;            // Algorithmic requirement
    params.limit  = params.npts * (1+w->k/2);   // Hopefully sufficient for..
    params.epsabs = GSL_DBL_EPSILON;            // ..either this...
    params.epsrel = GSL_SQRT_DBL_EPSILON;       // ..or this tolerance
    if (NULL == (params.iw = gsl_integration_workspace_alloc(params.limit))) {
        SUZERAIN_ERROR_REPORT("failed to allocate params.iw",
                              (status = SUZERAIN_ENOMEM));
    }
    if (NULL == (params.pts = malloc(params.npts * sizeof(params.pts[0])))) {
        SUZERAIN_ERROR_REPORT("failed to allocate params.pts",
                              (status = SUZERAIN_ENOMEM));
    }
    gsl_matrix * const dB = gsl_matrix_alloc(w->k, 3);
    if (NULL == dB) {
        SUZERAIN_ERROR_REPORT("failed to allocate dB",
                              (status = SUZERAIN_ENOMEM));
    }
    gsl_vector * const Bk = gsl_vector_alloc(w->k);
    if (NULL == Bk) {
        SUZERAIN_ERROR_REPORT("failed to allocate Bk",
                              (status = SUZERAIN_ENOMEM));
    }
    gsl_root_fsolver * const s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    if (NULL == s) {
        SUZERAIN_ERROR_REPORT("failed to allocate s",
                              (status = SUZERAIN_ENOMEM));
    }
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    // See suzerain_bl_compute_thicknesses for the constant baseflow version of
    // these quantities.  See writeups/thicknesses.pdf for how the integrals in
    // this function are defined and used.

    // Find implicitly-defined displacement thickness: thick->delta1
    // Displacement thickness more robust than delta99 or delta{2,3,H}
    {
        params_thickness_displacement integrand_params = {
            .vis_rhou   = coeffs_vis_rhou,
            .inv_stride = inv_stride,
            .inv_rhou   = coeffs_inv_rhou,
            .Bk         = Bk,
            .w          = w
        };
        params.integrand.function = &integrand_thickness_displacement;
        params.integrand.params   = &integrand_params;
        gsl_function f            = { &integral_thickness_residual, &params };

        // Look across the whole domain, usually [0, oo), for the thickness
        double lower = params.dis[0];
        double upper = params.dis[params.ndis-1];
        int tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->delta1);
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

    // As \delta_{99} must always be outside \delta_1 and the latter
    // is more robust than the former, use \delta_1 as a lower bound on
    // where we might find \delta_{99}.
    if (!gsl_isnan(thick->delta1)) {
        int tmp = suzerain_bl_find_edge(
                coeffs_vis_H0, thick->delta1,
                gsl_bspline_breakpoint(w->nbreak-1, w),
                &thick->delta, dB, w, dw);
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }

    // Find implicitly-defined thickness: thick->delta1 + thick->delta2
    if (!gsl_isnan(thick->delta1)) {
        params_thickness_momentum integrand_params = {
            .vis_rhou    = coeffs_vis_rhou,
            .vis_u       = coeffs_vis_u,
            .inv_stride  = inv_stride,
            .inv_rhou    = coeffs_inv_rhou,
            .inv_u       = coeffs_inv_u,
            .Bk          = Bk,
            .w           = w
        };
        params.integrand.function = &integrand_thickness_momentum;
        params.integrand.params   = &integrand_params;
        gsl_function f            = { &integral_thickness_residual, &params };

        // Try interval [delta1, oo) which may yield nonnegative delta2...
        // (thus increasing our robustness/speed in the common delta2>0 case)
        double lower = thick->delta1;
        double upper = params.dis[params.ndis-1];
        int tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->delta2);
        if (SUZERAIN_UNLIKELY(tmp != SUZERAIN_SUCCESS)) {
            // ...but accept a negative thick->delta2 if observed.  So it goes.
            lower = params.dis[0];
            upper = thick->delta1;
            tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->delta2);
        }
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }
    // Adjust to obtain momentum thickness, thick->delta2
    thick->delta2 -= thick->delta1;

    // Find implicitly-defined energy thickness: thick->delta1 + thick->delta3
    if (!gsl_isnan(thick->delta1)) {
        params_thickness_energy integrand_params = {
            .vis_ke      = coeffs_vis_ke,
            .vis_rhou    = coeffs_vis_rhou,
            .inv_stride  = inv_stride,
            .inv_rhou    = coeffs_inv_rhou,
            .inv_u       = coeffs_inv_u,
            .inv_v       = coeffs_inv_v,
            .Bk          = Bk,
            .w           = w
        };
        params.integrand.function = &integrand_thickness_energy;
        params.integrand.params   = &integrand_params;
        gsl_function f            = { &integral_thickness_residual, &params };

        // Try interval [delta1, oo) which may yield nonnegative delta3...
        // (thus increasing our robustness/speed in the common delta3>0 case)
        double lower = thick->delta1;
        double upper = params.dis[params.ndis-1];
        int tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->delta3);
        if (SUZERAIN_UNLIKELY(tmp != SUZERAIN_SUCCESS)) {
            // ...but accept a negative thick->delta3 if observed.  So it goes.
            lower = params.dis[0];
            upper = thick->delta1;
            tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->delta3);
        }
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }
    // Adjust to obtain energy thickness, thick->delta3
    thick->delta3 -= thick->delta1;

    // Find implicitly-defined thickness: thick->delta1 + thick->deltah
    if (!gsl_isnan(thick->delta1)) {
        params_thickness_enthalpy integrand_params = {
            .code_Ma     = code_Ma,
            .vis_H0      = coeffs_vis_H0,
            .vis_ke      = coeffs_vis_ke,
            .vis_rhou    = coeffs_vis_rhou,
            .inv_stride  = inv_stride,
            .inv_H0      = coeffs_inv_H0,
            .inv_rhou    = coeffs_inv_rhou,
            .inv_u       = coeffs_inv_u,
            .inv_v       = coeffs_inv_v,
            .Bk          = Bk,
            .w           = w
        };
        params.integrand.function = &integrand_thickness_enthalpy;
        params.integrand.params   = &integrand_params;
        gsl_function f            = { &integral_thickness_residual, &params };

        // Try interval [delta1, oo) which may yield nonnegative deltah...
        // (thus increasing our robustness/speed in the common deltah>0 case)
        double lower = thick->delta1;
        double upper = params.dis[params.ndis-1];
        int tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->deltah);
        if (SUZERAIN_UNLIKELY(tmp != SUZERAIN_SUCCESS)) {
            // ...but accept a negative thick->deltah if observed.  So it goes.
            lower = params.dis[0];
            upper = thick->delta1;
            tmp = fsolver_solve(s, &f, maxiter, params.epsabs, params.epsrel,
                                &lower, &upper, &thick->deltah);
        }
        if (status == SUZERAIN_SUCCESS) status = tmp;
    }
    // Adjust to obtain enthalpy thickness, thick->deltah
    thick->deltah -= thick->delta1;

done:
    gsl_integration_workspace_free(params.iw);
    free(params.pts);
    gsl_matrix_free(dB);
    gsl_vector_free(Bk);
    gsl_root_fsolver_free(s);

    return status;
}
