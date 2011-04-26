/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * bspline.hpp: C++ wrappers for the C-based GSL and Suzerain B-spline APIs
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BSPLINE_HPP
#define __SUZERAIN_BSPLINE_HPP

#include <gsl/gsl_bspline.h>
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
 * methods to simplify working with the GSL API.  Non-constant methods are not
 * thread-safe.
 */
class bspline : public boost::noncopyable {
public:

    /**
     * Create B-spline workspaces capable of evaluating piecewise polynomials
     * of degree \c k - 1 on \c nbreak breakpoints located at \c breakpoints.
     * The resulting basis will have <tt>nbreak + k - 2</tt> degrees of
     * freedom.
     *
     * @param k B-spline order per GSL/PPPACK conventions.
     *          For example, 4 denotes piecewise cubics.
     * @param nbreak   Number of breakpoints.
     * @param breakpts Strictly increasing breakpoint locations.
     *
     * @see gsl_bspline_alloc(), gsl_bspline_deriv_alloc(), gsl_bspline_knots()
     */
    bspline(int k, int nbreak, const double * breakpts)
        : bw(gsl_bspline_alloc(k, nbreak)),
          dbw(gsl_bspline_deriv_alloc(k)),
          db_(gsl_matrix_alloc(k, k))
    {
        gsl_vector_const_view view
            = gsl_vector_const_view_array(breakpts, nbreak);
        gsl_bspline_knots(&view.vector, bw);
    }

    /** @see gsl_bspline_free */
    ~bspline() {
        gsl_matrix_free(db_);
        gsl_bspline_deriv_free(dbw);
        gsl_bspline_free(bw);
    }

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
    double breakpoint(std::size_t i) const
    {
        return gsl_bspline_breakpoint(i, bw);
    }

    /** Retrieve the <tt>i</tt>th knot. */
    double knot(std::size_t i) const
    {
        return gsl_vector_get(bw->knots, i);
    }

    /**
     * Retrieve the <tt>i</tt>th collocation point.  These are the Greville
     * abscissae per <code>gsl_bspline_greville_abscissa()</code>.
     */
    double collocation_point(int i) const
    {
        return gsl_bspline_greville_abscissa(i, bw);
    }

    /**
     * @copybrief suzerain_bspline_linear_combination
     * @see       suzerain_bspline_linear_combination
     */
    int linear_combination(const std::size_t nderiv,
                           const double * coeffs,
                           const std::size_t npoints,
                           const double * points,
                           double * values,
                           const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination(nderiv, coeffs, npoints,
                points, values, ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_linear_combination_complex
     * @see       suzerain_bspline_linear_combination_complex
     */
    template< typename Complex1,
              typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>
    >, int>::type linear_combination(const std::size_t nderiv,
                                     const Complex1 *coeffs,
                                     const std::size_t npoints,
                                     const double * points,
                                     Complex2 *values,
                                     const std::size_t ldvalues = 0)
    {
        return suzerain_bspline_linear_combination_complex(
                nderiv,
                reinterpret_cast<const gsl_complex *>(coeffs),
                npoints, points,
                reinterpret_cast<gsl_complex *>(values),
                ldvalues, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_integration_coefficients
     * @see       suzerain_bspline_integration_coefficients
     */
    int integration_coefficients(const std::size_t nderiv,
                                 double * coeffs)
    {
        return suzerain_bspline_integration_coefficients(
                nderiv, coeffs, 1, db_, bw, dbw);
    }

    /**
     * @copybrief suzerain_bspline_integration_coefficients
     * @see       suzerain_bspline_integration_coefficients
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type integration_coefficients(const std::size_t nderiv,
                                     Complex *coeffs)
    {
        // Zero real and imaginary components
        for (std::size_t i = 0; i < bw->n; ++i) {
            suzerain::complex::assign_components(coeffs[i], 0, 0);
        }

        // Populate the real components
        return suzerain_bspline_integration_coefficients(
                nderiv, reinterpret_cast<double *>(coeffs),
                sizeof(Complex)/sizeof(double), db_, bw, dbw);
    }

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
    bsplineop(suzerain::bspline &b,
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
     * @copybrief suzerain_bsplineop_workspace#D
     * @see       suzerain_bsplineop_workspace#D
     */
    const double * D(int i) const { return w_->D[i]; }

    /**
     * @copybrief suzerain_bsplineop_workspace#D
     * @see       suzerain_bsplineop_workspace#D
     */
    double * D(int i) { return w_->D[i]; }

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
                   double alpha, const double *x, int incx, int ldx,
                   double beta, double *y, int incy, int ldy) const
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
                   double alpha, const double *x, int incx,
                   double beta, double *y, int incy) const
    {
        return suzerain_bsplineop_accumulate(nderiv, 1,
                                             alpha, x, incx, 0,
                                             beta,  y, incy, 0, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply
     * @see       suzerain_bsplineop_apply
     */
    int apply(int nderiv, int nrhs, double alpha,
              double *x, int incx, int ldx) const
    {
        return suzerain_bsplineop_apply(nderiv, nrhs, alpha, x, incx, ldx, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply
     * @see       suzerain_bsplineop_apply
     */
    int apply(int nderiv, double alpha,
              double *x, int incx) const
    {
        return suzerain_bsplineop_apply(nderiv, 1, alpha, x, incx, 0, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_interpolation_rhs
     * @see       suzerain_bsplineop_interpolation_rhs
     */
    int interpolation_rhs(const suzerain_function * function,
                          double * rhs,
                          suzerain::bspline &b) const
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
    template< typename Complex1,
              typename Complex2,
              typename Complex3,
              typename Complex4 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>,
        suzerain::complex::traits::is_complex_double<Complex3>,
        suzerain::complex::traits::is_complex_double<Complex4>
    >, int>::type accumulate(
                int nderiv, int nrhs,
                const Complex1 &alpha, const Complex2 *x, int incx, int ldx,
                const Complex3 &beta, Complex4 *y, int incy, int ldy) const
    {
        return suzerain_bsplineop_accumulate_complex(
                nderiv, nrhs,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double (*)[2]>(x),
                incx, ldx,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<double (*)[2]>(y),
                incy, ldy, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_accumulate_complex
     * @see       suzerain_bsplineop_accumulate_complex
     */
    template< typename Complex1,
              typename Complex2,
              typename Complex3,
              typename Complex4 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>,
        suzerain::complex::traits::is_complex_double<Complex3>,
        suzerain::complex::traits::is_complex_double<Complex4>
    >, int>::type accumulate(
                int nderiv,
                const Complex1 &alpha, const Complex2 *x, int incx,
                const Complex3 &beta, Complex4 *y, int incy) const
    {
        return suzerain_bsplineop_accumulate_complex(
                nderiv, 1,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double (*)[2]>(x),
                incx, 0,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<double (*)[2]>(y),
                incy, 0, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_accumulate_complex
     * @see       suzerain_bsplineop_accumulate_complex
     */
    template< typename Complex1, typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>
    >, int>::type accumulate(
                int nderiv, int nrhs,
                const double alpha, const Complex1 *x, int incx, int ldx,
                const double beta, Complex2 *y, int incy, int ldy) const
    {
        const double alpha_complex[2] = { alpha, 0 };
        const double beta_complex[2]  = { beta,  0 };
        return suzerain_bsplineop_accumulate_complex(
                nderiv, nrhs,
                alpha_complex,
                reinterpret_cast<const double (*)[2]>(x),
                incx, ldx,
                beta_complex,
                reinterpret_cast<double (*)[2]>(y),
                incy, ldy, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_accumulate_complex
     * @see       suzerain_bsplineop_accumulate_complex
     */
    template< typename Complex1, typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>
    >, int>::type accumulate(
                int nderiv,
                const double alpha, const Complex1 *x, int incx,
                const double beta, Complex2 *y, int incy) const
    {
        const double alpha_complex[2] = { alpha, 0 };
        const double beta_complex[2]  = { beta, 0 };
        return suzerain_bsplineop_accumulate_complex(
                nderiv, 1,
                alpha_complex,
                reinterpret_cast<const double (*)[2]>(x),
                incx, 0,
                beta_complex,
                reinterpret_cast<double (*)[2]>(y),
                incy, 0, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply_complex
     * @see       suzerain_bsplineop_apply_complex
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type apply(int nderiv, int nrhs, double alpha,
                  Complex *x, int incx, int ldx) const
    {
        return suzerain_bsplineop_apply_complex(
                nderiv, nrhs, alpha,
                reinterpret_cast<double (*)[2]>(x),
                incx, ldx, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_apply_complex
     * @see       suzerain_bsplineop_apply_complex
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type apply(int nderiv, double alpha,
                  Complex *x, int incx) const
    {
        return suzerain_bsplineop_apply_complex(
                nderiv, 1, alpha,
                reinterpret_cast<double (*)[2]>(x),
                incx, 0, w_);
    }

    /**
     * @copybrief suzerain_bsplineop_interpolation_rhs
     * @see       suzerain_bsplineop_interpolation_rhs
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type interpolation_rhs(const suzerain_zfunction * zfunction,
                              Complex * rhs,
                              suzerain::bspline &b) const
    {
        return suzerain_bsplineop_interpolation_rhs_complex(
                zfunction, reinterpret_cast<double (*)[2]>(rhs), b.bw, w_);
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
     * @copybrief suzerain_bsplineop_lu_workspace#A
     * @see       suzerain_bsplineop_lu_workspace#A
     */
    const double * A() const { return luw_->A; }

    /**
     * @copybrief suzerain_bsplineop_lu_workspace#A
     * @see       suzerain_bsplineop_lu_workspace#A
     */
    double * A() { return luw_->A; }

    /**
     * @copybrief suzerain_bsplineop_lu_form
     * @see       suzerain_bsplineop_lu_form
     */
    int form(int ncoefficients,
             const double * coefficients,
             const bsplineop &op)
    {
        return suzerain_bsplineop_lu_form(
                ncoefficients, coefficients, op.get(), luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_form_mass
     * @see       suzerain_bsplineop_lu_form_mass
     */
    int form_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_lu_form_mass(op.get(), luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_solve
     * @see       suzerain_bsplineop_lu_solve
     */
    int solve(int nrhs, double *b, int incb, int ldb) const
    {
        return suzerain_bsplineop_lu_solve(nrhs, b, incb, ldb, luw_);
    }

    /**
     * @copybrief suzerain_bsplineop_lu_rcond
     * @see       suzerain_bsplineop_lu_rcond
     */
    int rcond(double *rcond) const
    {
        return suzerain_bsplineop_lu_rcond(rcond, luw_);
    }

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
     * @copybrief suzerain_bsplineop_luz_workspace#A
     * @see       suzerain_bsplineop_luz_workspace#A
     */
    const double (*A() const)[2] { return luzw_->A; }

    /**
     * @copybrief suzerain_bsplineop_luz_workspace#A
     * @see       suzerain_bsplineop_luz_workspace#A
     */
    double (*A())[2] { return luzw_->A; }

    /**
     * @copybrief suzerain_bsplineop_luz_form
     * @see       suzerain_bsplineop_luz_form
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type form(int ncoefficients,
                 const Complex* coefficients,
                 const bsplineop &op)
    {
        return suzerain_bsplineop_luz_form(
                ncoefficients,
                reinterpret_cast<const double (*)[2]>(coefficients),
                op.get(), luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_form_mass
     * @see       suzerain_bsplineop_luz_form_mass
     */
    int form_mass(const bsplineop &op)
    {
        return suzerain_bsplineop_luz_form_mass(op.get(), luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_solve
     * @see       suzerain_bsplineop_luz_solve
     */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type solve(int nrhs, Complex *b, int incb, int ldb) const
    {
        return suzerain_bsplineop_luz_solve(
                nrhs, reinterpret_cast<double (*)[2]>(b), incb, ldb, luzw_);
    }

    /**
     * @copybrief suzerain_bsplineop_luz_rcond
     * @see       suzerain_bsplineop_luz_rcond
     */
    int rcond(double *rcond) const
    {
        return suzerain_bsplineop_luz_rcond(rcond, luzw_);
    }

private:
    suzerain_bsplineop_luz_workspace *luzw_; /**< The wrapped instance */
};

} // namespace suzerain

#endif // __SUZERAIN_BSPLINE_HPP
