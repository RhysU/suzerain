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
 * bspline.hpp: C++ wrappers for the C-based suzerain_bspline API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BSPLINE_HPP
#define __SUZERAIN_BSPLINE_HPP

#include <suzerain/bspline.h>
#include <suzerain/complex.hpp>

/** @file
 * Provides C++ wrappers for the C-based API in bspline.h.  In particular,
 * provides RAII semantics for opaque types.
 */

namespace suzerain {

/**
 * Provides a thin RAII wrapper for suzerain_bspline_workspace.
 * @see suzerain_bspline_workspace.
 */
class bspline : public boost::noncopyable {
public:

    /** @see suzerain_bspline_alloc */
    bspline(int order, int nderivatives, int nbreakpoints,
            const double * breakpoints,
            enum suzerain_bspline_method method
                = SUZERAIN_BSPLINE_COLLOCATION_GREVILLE)
        : w_(suzerain_bspline_alloc(order,
                                    nderivatives,
                                    nbreakpoints,
                                    breakpoints,
                                    method)) {}

    /** @see suzerain_bspline_free */
    ~bspline() { suzerain_bspline_free(w_); }

    /** @return The wrapped suzerain_bspline_workspace pointer. */
    const suzerain_bspline_workspace* get() const { return w_; }

    /** @see suzerain_bspline_workspace#ndof */
    int ndof() const { return w_->ndof; }

    /** @see suzerain_bspline_accumulate_operator */
    int accumulate_operator(
            int nderivative, int nrhs,
            double alpha, const double *x, int incx, int ldx,
            double beta, double *y, int incy, int ldy) const
    {
        return suzerain_bspline_accumulate_operator(nderivative, nrhs,
                                                    alpha, x, incx, ldx,
                                                    beta, y, incy, ldy, w_);
    }

    /** @see suzerain_bspline_zaccumulate_operator */
    template< typename Complex1,
              typename Complex2,
              typename Complex3,
              typename Complex4 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>,
        suzerain::complex::traits::is_complex_double<Complex3>,
        suzerain::complex::traits::is_complex_double<Complex4>
    >, int>::type accumulate_operator(
                int nderivative, int nrhs,
                const Complex1 &alpha, const Complex2 *x, int incx, int ldx,
                const Complex3 &beta, Complex4 *y, int incy, int ldy) const
    {
        return suzerain_bspline_zaccumulate_operator(
                nderivative, nrhs,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double (*)[2]>(x),
                incx, ldx,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<double (*)[2]>(y),
                incy, ldy, w_);
    }

    /** @see suzerain_bspline_zaccumulate_operator */
    template< typename Complex1, typename Complex2 >
    typename boost::enable_if<boost::mpl::and_<
        suzerain::complex::traits::is_complex_double<Complex1>,
        suzerain::complex::traits::is_complex_double<Complex2>
    >, int>::type accumulate_operator(
                int nderivative, int nrhs,
                const double alpha, const Complex1 *x, int incx, int ldx,
                const double beta, Complex2 *y, int incy, int ldy) const
    {
        const double alpha_complex[2] = { alpha, 0 };
        const double beta_complex[2]  = { beta, 0 };
        return suzerain_bspline_zaccumulate_operator(
                nderivative, nrhs,
                alpha_complex,
                reinterpret_cast<const double (*)[2]>(x),
                incx, ldx,
                beta_complex,
                reinterpret_cast<double (*)[2]>(y),
                incy, ldy, w_);
    }

    /** @see suzerain_bspline_apply_operator */
    int apply_operator(
            int nderivative, int nrhs, double *b, int incb, int ldb) const
    {
        return suzerain_bspline_apply_operator(nderivative, nrhs,
                                               b, incb, ldb, w_);
    }

    /** @see suzerain_bspline_zapply_operator */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type apply_operator(int nderivative, int nrhs,
                           Complex *b, int incb, int ldb) const
    {
        return suzerain_bspline_zapply_operator(
                nderivative, nrhs,
                reinterpret_cast<double (*)[2]>(b),
                incb, ldb, w_);
    }

    /** @see suzerain_bspline_evaluate */
    int evaluate(int nderivative,
            const double * coefficients, int npoints, const double * points,
            double * values, int ldvalues) const
    {
        return suzerain_bspline_evaluate(nderivative,
                coefficients, npoints, points,
                values, ldvalues, w_);
    }


    /** @see suzerain_bspline_find_interpolation_problem_rhs */
    int find_interpolation_problem_rhs(const suzerain_function * function,
                                       double * rhs) const
    {
        return suzerain_bspline_find_interpolation_problem_rhs(
                function, rhs, w_);
    }

    /** @see suzerain_bspline_collocation_point */
    int collocation_point(int j, double *x_j) const
    {
        return suzerain_bspline_collocation_point(j, x_j, w_);
    }

    /** @see suzerain_bspline_collocation_points */
    int collocation_points(double *x, int incx) const
    {
        return suzerain_bspline_collocation_points(x, incx, w_);
    }

private:
    suzerain_bspline_workspace *w_; /**< The wrapped instance */
};


/**
 * Provides a thin RAII wrapper for suzerain_bspline_lu_workspace.
 * @see suzerain_bspline_lu_workspace.
 */
class bspline_lu : public boost::noncopyable {
public:

    /** @see suzerain_bspline_lu_alloc */
    bspline_lu(const bspline &w)
        : luw_(suzerain_bspline_lu_alloc(w.get())) {}

    /** @see suzerain_bspline_lu_free */
    ~bspline_lu() { suzerain_bspline_lu_free(luw_); }

    /** @return The wrapped suzerain_bspline_lu_workspace pointer. */
    const suzerain_bspline_lu_workspace* get() const { return luw_; }

    /** @see suzerain_bspline_lu_workspace#ndof */
    int ndof() const { return luw_->ndof; }

    /** @see suzerain_bspline_lu_form_general */
    int form_general(int ncoefficients,
                     const double * coefficients,
                     const bspline &w)
    {
        return suzerain_bspline_lu_form_general(
                ncoefficients, coefficients, w.get(), luw_);
    }

    /** @see suzerain_bspline_lu_form_mass */
    int form_mass(const bspline &w)
    {
        return suzerain_bspline_lu_form_mass(w.get(), luw_);
    }

    /** @see suzerain_bspline_lu_solve */
    int solve(int nrhs, double *b, int incb, int ldb) const
    {
        return suzerain_bspline_lu_solve(nrhs, b, incb, ldb, luw_);
    }

private:
    suzerain_bspline_lu_workspace *luw_; /**< The wrapped instance */
};


/**
 * Provides a thin RAII wrapper for suzerain_bspline_luz_workspace.
 * @see suzerain_bspline_luz_workspace.
 */
class bspline_luz : public boost::noncopyable {
public:

    /** @see suzerain_bspline_luz_alloc */
    bspline_luz(const bspline &w)
        : luzw_(suzerain_bspline_luz_alloc(w.get())) {}

    /** @see suzerain_bspline_luz_free */
    ~bspline_luz() { suzerain_bspline_luz_free(luzw_); }

    /** @return The wrapped suzerain_bspline_luz_workspace pointer. */
    const suzerain_bspline_luz_workspace* get() const { return luzw_; }

    /** @see suzerain_bspline_luz_workspace#ndof */
    int ndof() const { return luzw_->ndof; }

    /** @see suzerain_bspline_luz_form_general */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type form_general(int ncoefficients,
                         const Complex* coefficients,
                         const bspline &w)
    {
        return suzerain_bspline_luz_form_general(
                ncoefficients,
                reinterpret_cast<const double (*)[2]>(coefficients),
                w.get(), luzw_);
    }

    /** @see suzerain_bspline_luz_form_mass */
    int form_mass(const bspline &w)
    {
        return suzerain_bspline_luz_form_mass(w.get(), luzw_);
    }

    /** @see suzerain_bspline_luz_solve */
    template< typename Complex >
    typename boost::enable_if<
        suzerain::complex::traits::is_complex_double<Complex>, int
    >::type solve(int nrhs, Complex *b, int ldb) const
    {
        return suzerain_bspline_luz_solve(
                nrhs, reinterpret_cast<const double (*)[2]>(b), ldb, luzw_);
    }

private:
    suzerain_bspline_luz_workspace *luzw_; /**< The wrapped instance */
};

} // namespace suzerain

#endif // __SUZERAIN_BSPLINE_HPP
