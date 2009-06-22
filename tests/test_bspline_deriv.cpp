#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <log4cxx/logger.h>

#include <boost/test/included/unit_test.hpp>

using namespace log4cxx;
using boost::shared_ptr;

LoggerPtr logger = Logger::getRootLogger();

void print_gsl_matrix(FILE * stream, gsl_matrix * m)
{
    for (size_t i = 0; i < m->size1; ++i) {
        for (size_t j = 0; j < m->size2; ++j) {
            fprintf(stream, " %14g", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

BOOST_AUTO_TEST_CASE( main_test )
{
/*
    const size_t k = 3;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
*/
/*
    const size_t k = 3;
    const double bpoint_data[]    = { 0.0, 2.0, 4.0 };
*/
    const size_t k = 2;
    const double bpoint_data[]    = { 0.0, 1.0, 2.0, 3.0 };

    const size_t nbreak = sizeof(bpoint_data)/sizeof(bpoint_data[0]);
    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(k);
    gsl_bspline_knots((const gsl_vector *) &bpoints, bw);
    const size_t ncoeff = gsl_bspline_ncoeffs(bw);

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < bw->knots->size; ++i) {
            LOG4CXX_DEBUG(logger, "knot[" << i << "] = "
                                  << gsl_vector_get(bw->knots, i));
        }
    }

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < ncoeff; ++i) {
            LOG4CXX_DEBUG(logger, "abscissa[" << i << "] = "
                                  << gsl_bspline_greville_abscissa(i, bw));
        }
    }

    const size_t nderiv = 2;
    gsl_matrix *M  = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *D1 = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *D2 = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *dB = gsl_matrix_alloc(ncoeff, nderiv+1);

    for (size_t j = 0; j < M->size1; ++j) {
        double x = gsl_bspline_greville_abscissa(j, bw);
        gsl_bspline_deriv_eval(x, nderiv, dB, bw, bdw);
        gsl_matrix_set_row(M,  j, (const gsl_vector *) &gsl_matrix_column(dB, 0));
        gsl_matrix_set_row(D1, j, (const gsl_vector *) &gsl_matrix_column(dB, 1));
        gsl_matrix_set_row(D2, j, (const gsl_vector *) &gsl_matrix_column(dB, 2));
    }

    LOG4CXX_INFO(logger, "M follows... ");
    print_gsl_matrix(stdout, M);

    LOG4CXX_INFO(logger, "D1 follows... ");
    print_gsl_matrix(stdout, D1);

    int signum;
    gsl_permutation *p = gsl_permutation_calloc(M->size1);
    gsl_matrix *MLU    = gsl_matrix_alloc(M->size1, M->size2);
    gsl_matrix_memcpy(MLU, M);
    gsl_linalg_LU_decomp(MLU, p, &signum);
    LOG4CXX_INFO(logger, "MLU follows... ");
    print_gsl_matrix(stdout, MLU);

    {
        // Vector of function coefficients for f(x) = x+1
        gsl_vector * c  = gsl_vector_calloc(ncoeff);
        for (size_t i = 0; i < ncoeff; ++i) {
            gsl_vector_set(c, i, i+1);
        }

        // Form rhs = D1 c
        gsl_vector * rhs  = gsl_vector_calloc(ncoeff);
        gsl_blas_dgemv(CblasNoTrans, 1.0, D1, c, 0.0, rhs);
        // Solve MLU dc = rhs
        gsl_vector *dc = gsl_vector_alloc(ncoeff);
        gsl_linalg_LU_solve(MLU, p, rhs, dc);

        if (logger->isDebugEnabled()) {
            for (size_t i = 0; i < dc->size; ++i) {
                LOG4CXX_DEBUG(logger, "dc[" << i << "] = "
                                      << gsl_vector_get(dc, i));
            }
        }

        // Check residual for M dc - D1 c = 0
        gsl_vector * residual = gsl_vector_alloc(ncoeff);
        gsl_blas_dgemv(CblasNoTrans, -1.0, D1, c, 0.0, residual);
        gsl_blas_dgemv(CblasNoTrans,  1.0, M, dc, 1.0, residual);
        if (logger->isTraceEnabled()) {
            for (size_t i = 0; i < residual->size; ++i) {
                LOG4CXX_TRACE(logger, "residual[" << i << "] = "
                                      << gsl_vector_get(residual, i));
            }
        }
        gsl_vector_free(residual);

        const int nsample = 5;
        const size_t nderiv = 2;
        for (size_t i = 0; i < ncoeff-1; ++i) {
            const double gai   = gsl_bspline_greville_abscissa(i, bw);
            const double gaip1 = gsl_bspline_greville_abscissa(i+1, bw);
            for (size_t j = 0; j < nsample; ++j) {
                double direct_d1, collocation_d1;
                const double x = gai + j/((double)nsample)*(gaip1-gai);
                gsl_bspline_deriv_eval(x, nderiv, dB, bw, bdw);
                gsl_blas_ddot(c,  (const gsl_vector *) &gsl_matrix_column(dB, 1),
                              &direct_d1);
                gsl_blas_ddot(dc, (const gsl_vector *) &gsl_matrix_column(dB, 0),
                              &collocation_d1);
                LOG4CXX_DEBUG(logger, "At x = " << x << " direct (" << direct_d1 << ") versus collocation (" << collocation_d1 << ") gives error " << direct_d1 - collocation_d1);
            }
        }

        gsl_vector_free(c);
        gsl_vector_free(dc);
        gsl_vector_free(rhs);
    }

    gsl_matrix_free(dB);
    gsl_matrix_free(MLU);
    gsl_permutation_free(p);
    gsl_matrix_free(M);
    gsl_matrix_free(D1);
    gsl_matrix_free(D2);
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);
}
