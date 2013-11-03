//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
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

#ifndef SUZERAIN_BSPLINE_HPP
#define SUZERAIN_BSPLINE_HPP

/** @file
 * C++ wrappers for the C-based GSL and Suzerain B-spline APIs
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <suzerain/bspline.h>
#include <suzerain/bsplineop.h>
#include <suzerain/complex.hpp>

/** @file
 *
 * Provides C++ RAII wrappers for the C-based <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html">
 * B-spline API</a>.
 */

namespace suzerain {

/**
 * Provides a thin RAII wrapper for a combined
 * <code>gsl_bspline_workspace</code> and
 * <code>gsl_bspline_deriv_workspace</code>.  The wrapper includes convenience
 * methods to simplify working with the GSL API.
 * Non-constant methods are not
 * thread-safe.
 */
class bspline : public boost::noncopyable {
public:

    /**
     * Constructor selection tags.  Constructors are selected by tag dispatch
     * because the usual named constructor idioms do not translate well to
     * noncopyable classes.
     * @{
     */

    struct from_breakpoints {};
    struct from_abscissae   {};

    /** @} */

    /**
     * Create a B-spline workspace capable of evaluating piecewise polynomials
     * of degree \c k - 1 on \c nbreak breakpoints located at \c breakpoints.
     * The resulting basis will have <tt>nbreak + k - 2</tt> degrees of
     * freedom.
     *
     * @param k        B-spline order per GSL/PPPACK conventions.
     *                 For example, 4 denotes piecewise cubics.
     * @param tag      Dispatch tag to identify this constructor.
     * @param nbreak   Number of breakpoints.
     * @param breakpts Strictly increasing breakpoint locations.
     *
     * @see gsl_bspline_alloc(), gsl_bspline_deriv_alloc(), and
     *      gsl_bspline_knots()
     */
    bspline(int k, const from_breakpoints &tag,
            int nbreak, const real_t * breakpts)
        : bw(gsl_bspline_alloc(k, nbreak)),
          dbw(gsl_bspline_deriv_alloc(k)),
          db_(gsl_matrix_alloc(k, k))
    {
        SUZERAIN_UNUSED(tag);
        gsl_vector_const_view view
            = gsl_vector_const_view_array(breakpts, nbreak);
        gsl_bspline_knots(&view.vector, bw);
    }

    /**
     * Create B-spline workspace capable of evaluating piecewise polynomials of
     * degree \c k - 1 on a knot vector created to best approximate the \c
     * nabscissae provided Greville \c abscissae.  The resulting basis will
     * have \c nabscissae degrees of freedom.
     *
     * @param[in]  k          B-spline order per GSL/PPPACK conventions.
     *       [in]             For example, 4 denotes piecewise cubics.
     * @param[in]  tag        Dispatch tag to identify this constructor.
     * @param[in]  nabscissae Number of abscissae.
     * @param[in]  abscissae  Strictly increasing Greville abscissae locations.
     * @param[out] abserr     Absolute error over all approximated \c abscissae
     *                        when not NULL.
     *
     * @see gsl_bspline_alloc(), gsl_bspline_deriv_alloc(), and
     *      gsl_bspline_knots_greville() with the latter detailing the best
     *      approximation requirements.
     */
    bspline(int k, const from_abscissae &tag,
            int nabscissae, const real_t * abscissae, real_t * abserr = NULL)
        : bw(gsl_bspline_alloc(k, /* nbreak == */ nabscissae - k + 2)),
          dbw(gsl_bspline_deriv_alloc(k)),
          db_(gsl_matrix_alloc(k, k))
    {
        SUZERAIN_UNUSED(tag);
        gsl_vector_const_view view
            = gsl_vector_const_view_array(abscissae, nabscissae);
        gsl_bspline_knots_greville(&view.vector, bw, abserr);
    }

    /** @see gsl_bspline_free */
    ~bspline() {
        gsl_matrix_free(db_);
        gsl_bspline_deriv_free(dbw);
        gsl_bspline_free(bw);
    }

/** @name General inquiry */
/**@{*/

    /** The wrapped GSL B-spline workspace. */
    gsl_bspline_workspace * const bw;

    /** The wrapped GSL B-spline derivative workspace. */
    gsl_bspline_deriv_workspace * const dbw;

    /**
     * @copybrief  suzerain_bsplineop_workspace#k
     * @see        suzerain_bsplineop_workspace#k
     */
    int k() const { return bw->k; }

    /** The number of breakpoints for the basis */
    int nbreak() const { return bw->nbreak; }

    /** The number of knots for the basis */
    int nknot() const { return bw->knots->size; }

    /**
     * @copybrief suzerain_bsplineop_workspace#n
     * @see       suzerain_bsplineop_workspace#n
     */
    int n() const { return bw->n; }

    /** Retrieve the <tt>i</tt>th breakpoint. */
    real_t breakpoint(std::size_t i) const
    {
        return gsl_bspline_breakpoint(i, bw);
    }

    /** Retrieve the <tt>i</tt>th knot. */
    real_t knot(std::size_t i) const
    {
        return gsl_vector_get(bw->knots, i);
    }

    /**
     * Retrieve the <tt>i</tt>th collocation point.  These are the Greville
     * abscissae per <code>gsl_bspline_greville_abscissa()</code>.
     */
    real_t collocation_point(int i) const
    {
        return gsl_bspline_greville_abscissa(i, bw);
    }

    /**
     * Retrieve the minimum separation between collocation points
     * <tt>[i-1,i]</tt> and <tt>[i,i+1]</tt>.
     */
    real_t spacing_collocation_point(int i) const
    {
        return suzerain_bspline_spacing_greville_abscissae(i, bw);
    }

    /**
     * Retrieve the minimum separation between the breakpoints
     * providing the support for the <tt>i</tt>th basis function.
     */
    real_t spacing_breakpoint(int i) const
    {
        return suzerain_bspline_spacing_breakpoints(i, bw);
    }

/**@}*/

/** @name Operations */
/**@{*/

    /**
     * @copybrief suzerain_bspline_linear_combination
     * @see       suzerain_bspline_linear_combination
     */
    int linear_combination(const std::size_t nderiv,
                           const real_t * coeffs,
                           const std::size_t npoints,
                           const real_t * points,
                           real_t * values,
                           const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination(nderiv, coeffs, npoints,
                points, values, ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_linear_combination
     * @see       suzerain_bspline_linear_combination
     */
    int linear_combination(const std::size_t nderiv,
                           const real_t * coeffs,
                           const real_t   point,
                           real_t * values,
                           const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination(nderiv, coeffs,
                1, &point, // Evaluate at one single point
                values, ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_linear_combination_complex
     * @see       suzerain_bspline_linear_combination_complex
     */
    template< typename Complex1,
              typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_t<Complex1>,
        suzerain::complex::traits::is_complex_t<Complex2>
    >, int>::type linear_combination(const std::size_t nderiv,
                                     const Complex1 *coeffs,
                                     const std::size_t npoints,
                                     const real_t * points,
                                     Complex2 *values,
                                     const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination_complex(
                nderiv,
                reinterpret_cast<const complex_t *>(coeffs),
                npoints, points,
                reinterpret_cast<complex_t *>(values),
                ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_linear_combination_complex
     * @see       suzerain_bspline_linear_combination_complex
     */
    template< typename Complex1,
              typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_t<Complex1>,
        suzerain::complex::traits::is_complex_t<Complex2>
    >, int>::type linear_combination(const std::size_t nderiv,
                                     const Complex1 *coeffs,
                                     const real_t point,
                                     Complex2 *values,
                                     const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination_complex(
                nderiv,
                reinterpret_cast<const complex_t *>(coeffs),
                1, &point, // Evaluate at one single point
                reinterpret_cast<complex_t *>(values),
                ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_crossing
     * @see       suzerain_bspline_crossing
     */
    int crossing(const size_t nderiv,
                 const double * coeffs,
                 const double value,
                 double * lower,
                 double * upper,
                 const size_t maxiter,
                 const double epsabs,
                 const double epsrel,
                 double * location)
    {
        return suzerain_bspline_crossing(
                nderiv, coeffs, value, lower, upper,
                maxiter, epsabs, epsrel, location, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_integration_coefficients
     * @see       suzerain_bspline_integration_coefficients
     */
    int integration_coefficients(
            const std::size_t nderiv,
            real_t * coeffs,
            const real_t lo = -std::numeric_limits<real_t>::max(),
            const real_t hi =  std::numeric_limits<real_t>::max())
    {
        return suzerain_bspline_integration_coefficients(
                nderiv, coeffs, 1, lo, hi, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_integration_coefficients
     * @see       suzerain_bspline_integration_coefficients
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type integration_coefficients(
            const std::size_t nderiv,
            Complex *coeffs,
            const real_t lo = -std::numeric_limits<real_t>::max(),
            const real_t hi =  std::numeric_limits<real_t>::max())
    {
        // Zero real and imaginary components
        for (std::size_t i = 0; i < bw->n; ++i) {
            suzerain::complex::assign_components(coeffs[i], 0, 0);
        }

        // Populate the real components
        return suzerain_bspline_integration_coefficients(
                nderiv, reinterpret_cast<real_t *>(coeffs),
                sizeof(complex_t)/sizeof(real_t), lo, hi, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_distance
     * @see       suzerain_bspline_distance
     */
    real_t distance_to(const bspline& other) const
    {
        return suzerain_bspline_distance(this->bw, other.bw);
    }

/**@}*/

private:
    gsl_matrix *db_; /**< Scratch storage for evaluation routines */
};

/**
 * Provides a thin RAII wrapper for ::suzerain_bsplineop_workspace.
 * Non-constant methods are not thread-safe.
 */
class bsplineop : public boost::noncopyable {
public:

    /** @see suzerain_bsplineop_alloc */
    bsplineop(bspline &b,
              int nderiv,
              enum suzerain_bsplineop_method method
                    = SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE)
        : w_(suzerain_bsplineop_alloc(b.bw, b.dbw, nderiv, method)) {}

    ~bsplineop() { suzerain_bsplineop_free(w_); }

/** @name General inquiry */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_workspace#k
     * @see       suzerain_bsplineop_workspace#k
     */
    int k() const { return w_->k; }

    /**
     * @copybrief suzerain_bsplineop_workspace#n
     * @see       suzerain_bsplineop_workspace#n
     */
    int n() const { return w_->n; }

    /**
     * @copybrief suzerain_bsplineop_workspace#nderiv
     * @see       suzerain_bsplineop_workspace#nderiv
     */
    int nderiv() const { return w_->nderiv; }

    /**
     * @copybrief suzerain_bsplineop_workspace#kl
     * @see       suzerain_bsplineop_workspace#kl
     */
    int kl(int i) const { return w_->kl[i]; }

    /**
     * @copybrief suzerain_bsplineop_workspace#ku
     * @see       suzerain_bsplineop_workspace#ku
     */
    int ku(int i) const { return w_->ku[i]; }

    /**
     * @copybrief suzerain_bsplineop_workspace#max_kl
     * @see       suzerain_bsplineop_workspace#max_kl
     */
    int max_kl() const { return w_->max_kl; }

    /**
     * @copybrief suzerain_bsplineop_workspace#max_ku
     * @see       suzerain_bsplineop_workspace#max_ku
     */
    int max_ku() const { return w_->max_ku; }

    /**
     * @copybrief suzerain_bsplineop_workspace#ld
     * @see       suzerain_bsplineop_workspace#ld
     */
    int ld() const { return w_->ld; }

    /**
     * @copybrief suzerain_bsplineop_workspace#D_T
     * @see       suzerain_bsplineop_workspace#D_T
     */
    const real_t * D_T(int i) const { return w_->D_T[i]; }

    /**
     * @copybrief suzerain_bsplineop_workspace#D_T
     * @see       suzerain_bsplineop_workspace#D_T
     */
    real_t * D_T(int i) { return w_->D_T[i]; }

    /** @return The wrapped suzerain_bsplineop_workspace pointer. */
    const suzerain_bsplineop_workspace* get() const { return w_; }

/**@}*/

/** @name Real-valued operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_accumulate
     * @see       suzerain_bsplineop_accumulate
     */
    int accumulate(int nderiv, int nrhs,
                   real_t alpha, const real_t *x, int incx, int ldx,
                   real_t beta, real_t *y, int incy, int ldy) const
    {
        return suzerain_bsplineop_accumulate(nderiv, nrhs,
                                             alpha, x, incx, ldx,
                                             beta, y, incy, ldy, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_accumulate
     * @see       suzerain_bsplineop_accumulate
     */
    int accumulate(int nderiv,
                   real_t alpha, const real_t *x, int incx,
                   real_t beta,        real_t *y, int incy) const
    {
        return this->accumulate(
                nderiv, 1, alpha, x, incx, 0, beta,  y, incy, 0);
    }

    /**
     * @copybrief suzerain_bsplineop_apply
     * @see       suzerain_bsplineop_apply
     */
    int apply(int nderiv, int nrhs, real_t alpha,
              real_t *x, int incx, int ldx) const
    {
        return suzerain_bsplineop_apply(nderiv, nrhs, alpha, x, incx, ldx, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply
     * @see       suzerain_bsplineop_apply
     */
    int apply(int nderiv, real_t alpha, real_t *x, int incx) const
    {
        return this->apply(nderiv, 1, alpha, x, incx, 0);
    }

    /**
     * @copybrief suzerain_bsplineop_interpolation_rhs
     * @see       suzerain_bsplineop_interpolation_rhs
     */
    int interpolation_rhs(const suzerain_function * function,
                          real_t * rhs,
                          bspline &b) const
    {
        return suzerain_bsplineop_interpolation_rhs(function, rhs, b.bw, w_);
    }

/**@}*/

/** @name Complex-valued operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_accumulate_complex
     * @see       suzerain_bsplineop_accumulate_complex
     */
    template< typename AlphaType, typename BetaType,
              typename Complex1,  typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_t<Complex1>,
        suzerain::complex::traits::is_complex_t<Complex2>
    >, int>::type accumulate(
            int nderiv, int nrhs,
            const AlphaType &alpha, const Complex1 *x, int incx, int ldx,
            const BetaType  &beta,        Complex2 *y, int incy, int ldy) const
    {
        complex_t alpha_complex;
        suzerain::complex::assign_complex(alpha_complex, alpha);
        complex_t beta_complex;
        suzerain::complex::assign_complex(beta_complex, beta);
        return suzerain_bsplineop_accumulate_complex(
                nderiv, nrhs,
                alpha_complex,
                reinterpret_cast<const complex_t *>(x),
                incx, ldx,
                beta_complex,
                reinterpret_cast<complex_t *>(y),
                incy, ldy, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_accumulate_complex
     * @see       suzerain_bsplineop_accumulate_complex
     */
    template< typename AlphaType, typename BetaType,
              typename Complex1,  typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_t<Complex1>,
        suzerain::complex::traits::is_complex_t<Complex2>
    >, int>::type accumulate(
                int nderiv,
                const AlphaType &alpha, const Complex1 *x, int incx,
                const BetaType  &beta,        Complex2 *y, int incy) const
    {
        return this->accumulate(nderiv, 1, alpha, x, incx, 0, beta, y, incy, 0);
    }

    /**
     * @copybrief suzerain_bsplineop_apply_complex
     * @see       suzerain_bsplineop_apply_complex
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type apply(int nderiv, int nrhs, real_t alpha,
                  Complex *x, int incx, int ldx) const
    {
        return suzerain_bsplineop_apply_complex(
                nderiv, nrhs, alpha,
                reinterpret_cast<complex_t *>(x),
                incx, ldx, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply_complex
     * @see       suzerain_bsplineop_apply_complex
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type apply(int nderiv, real_t alpha,
                  Complex *x, int incx) const
    {
        return this->apply(nderiv, 1, alpha, x, incx, 0);
    }

    /**
     * @copybrief suzerain_bsplineop_interpolation_rhs
     * @see       suzerain_bsplineop_interpolation_rhs
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type interpolation_rhs(const suzerain_zfunction * zfunction,
                              Complex * rhs,
                              bspline &b) const
    {
        return suzerain_bsplineop_interpolation_rhs_complex(
                zfunction, reinterpret_cast<complex_t *>(rhs),
                b.bw, w_);
    }

/**@}*/

private:
    suzerain_bsplineop_workspace *w_; /**< The wrapped instance */
};


/**
 * Provides a thin RAII wrapper for suzerain_bsplineop_lu_workspace.
 * @see suzerain_bsplineop_lu_workspace.
 */
class bsplineop_lu : public boost::noncopyable {
public:

    /** @see suzerain_bsplineop_lu_alloc */
    bsplineop_lu(const bsplineop &op)
        : luw_(suzerain_bsplineop_lu_alloc(op.get())) {}

    /** @see suzerain_bsplineop_lu_free */
    ~bsplineop_lu() { suzerain_bsplineop_lu_free(luw_); }

    /** @return The wrapped suzerain_bsplineop_lu_workspace pointer. */
    const suzerain_bsplineop_lu_workspace* get() const { return luw_; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#n
     * @see       suzerain_bsplineop_lu_workspace#n
     */
    int n() const { return luw_->n; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#kl
     * @see       suzerain_bsplineop_lu_workspace#kl
     */
    int kl() const { return luw_->kl; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#ku
     * @see       suzerain_bsplineop_lu_workspace#ku
     */
    int ku() const { return luw_->ku; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#ld
     * @see       suzerain_bsplineop_lu_workspace#ld
     */
    int ld() const { return luw_->ld; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#ipiv
     * @see       suzerain_bsplineop_lu_workspace#ipiv
     */
    const int * ipiv() const { return luw_->ipiv; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#ipiv
     * @see       suzerain_bsplineop_lu_workspace#ipiv
     */
    int * ipiv() { return luw_->ipiv; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#A_T
     * @see       suzerain_bsplineop_lu_workspace#A_T
     */
    const real_t * A_T() const { return luw_->A_T; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#A_T
     * @see       suzerain_bsplineop_lu_workspace#A_T
     */
    real_t * A_T() { return luw_->A_T; }

/** @name Operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_lu_opaccumulate
     * @see       suzerain_bsplineop_lu_opaccumulate
     */
    int opaccumulate(int ncoefficients,
                     const real_t * coefficients,
                     const bsplineop &op,
                     const real_t scale) const
    {
        return suzerain_bsplineop_lu_opaccumulate(
                ncoefficients, coefficients, op.get(), scale, luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_opnorm
     * @see       suzerain_bsplineop_lu_opnorm
     */
    int opnorm(real_t &norm) const
    {
        return suzerain_bsplineop_lu_opnorm(luw_, &norm);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_factor
     * @see       suzerain_bsplineop_lu_factor
     */
    int factor()
    {
        return suzerain_bsplineop_lu_factor(luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_rcond
     * @see       suzerain_bsplineop_lu_rcond
     */
    int rcond(const real_t norm, real_t &rcond) const
    {
        return suzerain_bsplineop_lu_rcond(norm, &rcond, luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_solve
     * @see       suzerain_bsplineop_lu_solve
     */
    int solve(int nrhs, real_t *B, int incb, int ldb) const
    {
        return suzerain_bsplineop_lu_solve(nrhs, B, incb, ldb, luw_);
    }

/**@}*/

/** @name Convenience operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_lu_opform
     * @see       suzerain_bsplineop_lu_opform
     */
    int opform(int ncoefficients,
               const real_t * coefficients,
               const bsplineop &op)
    {
        return suzerain_bsplineop_lu_opform(
                ncoefficients, coefficients, op.get(), luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_opform_mass
     * @see       suzerain_bsplineop_lu_opform_mass
     */
    int opform_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_lu_opform_mass(op.get(), luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_factor_mass
     * @see       suzerain_bsplineop_lu_factor_mass
     */
    int factor_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_lu_factor_mass(op.get(), luw_);
    }

/**@}*/

private:
    suzerain_bsplineop_lu_workspace *luw_; /**< The wrapped instance */
};


/**
 * Provides a thin RAII wrapper for suzerain_bsplineop_luz_workspace.
 * @see suzerain_bsplineop_luz_workspace.
 */
class bsplineop_luz : public boost::noncopyable {
public:

    /** @see suzerain_bsplineop_luz_alloc */
    bsplineop_luz(const bsplineop &op)
        : luzw_(suzerain_bsplineop_luz_alloc(op.get())) {}

    /** @see suzerain_bsplineop_luz_free */
    ~bsplineop_luz() { suzerain_bsplineop_luz_free(luzw_); }

/** @name General inquiry */
/**@{*/

    /** @return The wrapped suzerain_bsplineop_luz_workspace pointer. */
    const suzerain_bsplineop_luz_workspace* get() const { return luzw_; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#n
     * @see       suzerain_bsplineop_luz_workspace#n
     */
    int n() const { return luzw_->n; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#kl
     * @see       suzerain_bsplineop_luz_workspace#kl
     */
    int kl() const { return luzw_->kl; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#ku
     * @see       suzerain_bsplineop_luz_workspace#ku
     */
    int ku() const { return luzw_->ku; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#ld
     * @see       suzerain_bsplineop_luz_workspace#ld
     */
    int ld() const { return luzw_->ld; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#ipiv
     * @see       suzerain_bsplineop_luz_workspace#ipiv
     */
    const int * ipiv() const { return luzw_->ipiv; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#ipiv
     * @see       suzerain_bsplineop_luz_workspace#ipiv
     */
    int * ipiv() { return luzw_->ipiv; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#A_T
     * @see       suzerain_bsplineop_luz_workspace#A_T
     */
    const complex_t * A_T() const { return luzw_->A_T; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#A_T
     * @see       suzerain_bsplineop_luz_workspace#A_T
     */
    complex_t * A_T() { return luzw_->A_T; }

/**@}*/

/** @name Operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_luz_opaccumulate
     * @see       suzerain_bsplineop_luz_opaccumulate
     */
    template< typename ScaleType, typename Complex1 >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex1>, int
    >::type opaccumulate(int ncoefficients,
                         const Complex1* coefficients,
                         const bsplineop &op,
                         const ScaleType &scale)
    {
        complex_t scale_complex;
        suzerain::complex::assign_complex(scale_complex, scale);
        return suzerain_bsplineop_luz_opaccumulate(
                ncoefficients,
                reinterpret_cast<const complex_t *>(coefficients),
                op.get(), scale_complex, luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_opnorm
     * @see       suzerain_bsplineop_luz_opnorm
     */
    int opnorm(real_t &norm) const
    {
        return suzerain_bsplineop_luz_opnorm(luzw_, &norm);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_factor
     * @see       suzerain_bsplineop_luz_factor
     */
    int factor()
    {
        return suzerain_bsplineop_luz_factor(luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_rcond
     * @see       suzerain_bsplineop_luz_rcond
     */
    int rcond(const real_t norm, real_t &rcond) const
    {
        return suzerain_bsplineop_luz_rcond(norm, &rcond, luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_solve
     * @see       suzerain_bsplineop_luz_solve
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type solve(int nrhs, Complex *B, int incb, int ldb) const
    {
        return suzerain_bsplineop_luz_solve(
                nrhs, reinterpret_cast<complex_t *>(B),
                incb, ldb, luzw_);
    }

/**@}*/

/** @name Convenience operations */
/**@{*/

    /**
     * @copybrief suzerain_bsplineop_luz_opform
     * @see       suzerain_bsplineop_luz_opform
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_t<Complex>, int
    >::type opform(int ncoefficients,
                   const Complex* coefficients,
                   const bsplineop &op)
    {
        return suzerain_bsplineop_luz_opform(
                ncoefficients,
                reinterpret_cast<const complex_t *>(coefficients),
                op.get(), luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_opform_mass
     * @see       suzerain_bsplineop_luz_opform_mass
     */
    int opform_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_luz_opform_mass(op.get(), luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_factor_mass
     * @see       suzerain_bsplineop_luz_factor_mass
     */
    int factor_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_luz_factor_mass(op.get(), luzw_);
    }

/**@}*/

private:
    suzerain_bsplineop_luz_workspace *luzw_; /**< The wrapped instance */
};

} // namespace suzerain

#endif // SUZERAIN_BSPLINE_HPP
