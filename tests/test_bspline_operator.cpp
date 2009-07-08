#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <gsl/gsl_poly.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline_operator.h>
#include <suzerain/function.h>

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int order  = 4;
    const int nderiv = 2;
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);

    suzerain_bspline_operator_workspace *w
        = suzerain_bspline_operator_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    suzerain_bspline_operator_lu_workspace *luw
        = suzerain_bspline_operator_lu_alloc(w);
    suzerain_bspline_operator_lu_free(luw);

    suzerain_bspline_operator_free(w);
}

// Check a simple piecewise linear case's general banded storage
BOOST_AUTO_TEST_CASE( piecewise_linear_memory_application_solution )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 2;

    suzerain_bspline_operator_workspace *w
        = suzerain_bspline_operator_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *   1   0   0   0
         *   0   1   0   0
         *   0   0   1   0
         *   0   0   0   1
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = { /*DK*/   1,     0,
                                       0,   1,     0,
                                       0,   1,     0,
                                       0,   1  /*DK*/ };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
            w->D[0] + w->ku, w->D[0] + w->storagesize - w->kl);

        /* Check w->D[0] application against multiple vectors */
        const int nrhs = 2;
        double vector[] = { 1, 2, 3, 4,
                            5, 6, 7, 8 };
        const double good_result[] = { 1, 2, 3, 4,
                                       5, 6, 7, 8 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_operator_apply(0, nrhs, vector, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

    {
        /* Check w->D[1], the first derivative matrix, against known good:
         *  -1   1   0   0
         *   0  -1   1   0
         *   0   0  -1   1
         *   0   0  -1   1
         * Known good is in general banded matrix column-major order.
         */
        const double good_D1[] = { /*DK*/  -1,     0,
                                       1,  -1,     0,
                                       1,  -1,    -1,
                                       1,   1  /*DK*/ };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
            w->D[1] + w->ku, w->D[1] + w->storagesize - w->kl);

        /* Check w->D[0] application against multiple vectors */
        const int nrhs = 2;
        double vector[] = { 1, 3, 2, 4,
                            7, 6, 5, 8 };
        const double good_result[] = {  2, -1, 2, 2,
                                       -1, -1, 3, 3 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_operator_apply(1, nrhs, vector, ldb, w);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
            vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

    {
        /* Check w->D[2], the second derivative matrix, against zero result.
         */
        const double good_D2[] = { /*DK*/  0,    0,
                                       0,  0,    0,
                                       0,  0,    0,
                                       0,  0  /*DK*/ };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D2, good_D2 + sizeof(good_D2)/sizeof(good_D2[0]),
            w->D[2] + w->ku, w->D[2] + w->storagesize - w->kl);
    }

    suzerain_bspline_operator_lu_workspace *luw
        = suzerain_bspline_operator_lu_alloc(w);

    /* Form 2*D[0] - 3*D[1] operator in LU-ready banded storage.  Answer is
     *   5   -3    0     0
     *   0    5   -3     0
     *   0    0    5    -3
     *   0    0    3    -1
     * which, in LU-form where L has ones on the main diagonal, is
     *   5   -3    0     0
     *   0    5   -3     0
     *   0    0    5    -3
     *   0    0    0.6   0.8
     * The pivot matrix is eye(4).  Check it in octave using [l,u,p] = lu(A).
     * Known good is in general banded matrix column-major order with
     * additional superdiagonal to allow for LU factorization fill-in.
     */
    const double good_A0[] = { /*DK*/   /*DK*/   5,     0,
                               /*DK*/0,    -3,   5,     0,
                                     0,    -3,   5,     0.6,
                                     0,    -3,   0.8  /*DK*/ };
    const double coeff[] = { 2.0, -3.0 };
    suzerain_bspline_operator_lu_form(
        sizeof(coeff)/sizeof(coeff[0]), coeff, w, luw);
    {
        // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
        const double *expected, *actual;
        for (expected = good_A0, actual = luw->A + luw->ku;
             expected < good_A0 + sizeof(good_A0)/sizeof(good_A0[0]);
             ++expected, ++actual) {
            BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
        }
    }

    /* Check that multiple rhs solution works for operator found just above */
    {
        const int nrhs = 2;
        double vector[] = { 1,  2, 3, 4,
                           -4, -1, 1, 3};
        const double good_result[] = {  1.25, 1.75, 2.25, 2.75,
                                       -0.2,  1,    2,    3 };
        const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
        suzerain_bspline_operator_lu_solve(nrhs, vector, ldb, luw);
        {
            // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
            const double *expected, *actual;
            for (expected = good_result, actual = vector;
                 expected < good_result + sizeof(good_result)/sizeof(good_result[0]);
                 ++expected, ++actual) {
                BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
            }
        }
    }

    suzerain_bspline_operator_lu_free(luw);
    suzerain_bspline_operator_free(w);
}

// Check a simple piecewise quadratic case's general banded storage
BOOST_AUTO_TEST_CASE( piecewise_quadratic_memory_application_solution )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 3;
    const int nderiv = 2;

    suzerain_bspline_operator_workspace *w
        = suzerain_bspline_operator_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    {
        /* Check w->D[0], the mass matrix, against known good solution:
         *  1.    0.    0.    0.    0.
         *  1./4. 5./8. 1./8. 0.    0.
         *  0.    1./8. 3./4. 1./8. 0.
         *  0.    0.    1./8. 5./8. 1./4.
         *  0.    0.    0.    0.    1.
         * Known good is in general banded matrix column-major order.
         */
        const double good_D0[] = { /*DK*/      1.,       1./4.,
                                      0.,      5./8.,    1./8.,
                                      1./8.,   3./4.,    1./8.,
                                      1./8.,   5./8.,    0.,
                                      1./4.,       1.     /*DK*/ };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
            w->D[0] + w->ku, w->D[0] + w->storagesize - w->kl);

//      /* Check w->D[0] application against multiple vectors */
//      const int nrhs = 2;
//      double vector[] = { 1, 2, 3, 4,
//                          5, 6, 7, 8 };
//      const double good_result[] = { 1, 2, 3, 4,
//                                     5, 6, 7, 8 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_apply(0, nrhs, vector, ldb, w);
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
//          vector, vector + sizeof(vector)/sizeof(vector[0]));
    }

//  {
//      /* Check w->D[1], the first derivative matrix, against known good:
//       *  -1   1   0   0
//       *   0  -1   1   0
//       *   0   0  -1   1
//       *   0   0  -1   1
//       * Known good is in general banded matrix column-major order.
//       */
//      const double good_D1[] = { /*DK*/  -1,     0,
//                                     1,  -1,     0,
//                                     1,  -1,    -1,
//                                     1,   1  /*DK*/ };
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
//          w->D[1] + w->ku, w->D[1] + w->storagesize - w->kl);

//      /* Check w->D[0] application against multiple vectors */
//      const int nrhs = 2;
//      double vector[] = { 1, 3, 2, 4,
//                          7, 6, 5, 8 };
//      const double good_result[] = {  2, -1, 2, 2,
//                                     -1, -1, 3, 3 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_apply(1, nrhs, vector, ldb, w);
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
//          vector, vector + sizeof(vector)/sizeof(vector[0]));
//  }

//  {
//      /* Check w->D[2], the second derivative matrix, against zero result.
//       */
//      const double good_D2[] = { /*DK*/  0,    0,
//                                     0,  0,    0,
//                                     0,  0,    0,
//                                     0,  0  /*DK*/ };
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_D2, good_D2 + sizeof(good_D2)/sizeof(good_D2[0]),
//          w->D[2] + w->ku, w->D[2] + w->storagesize - w->kl);
//  }

//  suzerain_bspline_operator_lu_workspace *luw
//      = suzerain_bspline_operator_lu_alloc(w);

//  /* Form 2*D[0] - 3*D[1] operator in LU-ready banded storage.  Answer is
//   *   5   -3    0     0
//   *   0    5   -3     0
//   *   0    0    5    -3
//   *   0    0    3    -1
//   * which, in LU-form where L has ones on the main diagonal, is
//   *   5   -3    0     0
//   *   0    5   -3     0
//   *   0    0    5    -3
//   *   0    0    0.6   0.8
//   * The pivot matrix is eye(4).  Check it in octave using [l,u,p] = lu(A).
//   * Known good is in general banded matrix column-major order with
//   * additional superdiagonal to allow for LU factorization fill-in.
//   */
//  const double good_A0[] = { /*DK*/   /*DK*/   5,     0,
//                             /*DK*/0,    -3,   5,     0,
//                                   0,    -3,   5,     0.6,
//                                   0,    -3,   0.8  /*DK*/ };
//  const double coeff[] = { 2.0, -3.0 };
//  suzerain_bspline_operator_lu_form(
//      sizeof(coeff)/sizeof(coeff[0]), coeff, w, luw);
//  {
//      // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
//      const double *expected, *actual;
//      for (expected = good_A0, actual = luw->A + luw->ku;
//           expected < good_A0 + sizeof(good_A0)/sizeof(good_A0[0]);
//           ++expected, ++actual) {
//          BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
//      }
//  }

//  /* Check that multiple rhs solution works for operator found just above */
//  {
//      const int nrhs = 2;
//      double vector[] = { 1,  2, 3, 4,
//                         -4, -1, 1, 3};
//      const double good_result[] = {  1.25, 1.75, 2.25, 2.75,
//                                     -0.2,  1,    2,    3 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_lu_solve(nrhs, vector, ldb, luw);
//      {
//          // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
//          const double *expected, *actual;
//          for (expected = good_result, actual = vector;
//               expected < good_result + sizeof(good_result)/sizeof(good_result[0]);
//               ++expected, ++actual) {
//              BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
//          }
//      }
//  }

//  suzerain_bspline_operator_lu_free(luw);
    suzerain_bspline_operator_free(w);
}


// Check a piecewise cubic case's general banded storage
// BOOST_AUTO_TEST_CASE( piecewise_cubic_memory_application_solution )
// {
//     const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
//     const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
//     const int order  = 4;
//     const int nderiv = 2;
//
//     suzerain_bspline_operator_workspace *w
//         = suzerain_bspline_operator_alloc(order, nderiv, nbreak, breakpoints,
//             SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);
//
//     {
//         /* Check w->D[0], the mass matrix, against known good solution:
//          * 1.       0.        0.       0.        0.       0.
//          * 8./27.  61./108.  43./324.  1./162.   0.       0.
//          * 0.       1./4.     7./12.   1./6.     0.       0.
//          * 0.       0.        1./6.    7./12.    1./4.    0.
//          * 0.       0.        1./162. 43./324.  61./108.  8./27.
//          * 0.       0.        0.       0.        0.       1.
//          * Known good is in general banded matrix column-major order.
//          */
//         const double good_D0[] = {
//             /*DK*/      /*DK*/          1.,       8./27.,         0,
//             /*DK*/0,        0,         61./108.,  1./4.,          0,
//                   0,       43./324.,    7./12.,   1./6.,          1./162.,
//                   1./162.,  1./6.,      7./12.,   43./324.,       0,
//                   0,        1./4.,     61./108.,  0,        /*DK*/0,
//                   0,        8./27.,     1.        /*DK*/    /*DK*/
//         };
//     {
//         // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
//         int i = 0; /* DEBUG */
//         const double *expected, *actual;
//         for (expected = good_D0, actual = w->D[0] + w->ku;
//              expected < good_D0 + sizeof(good_D0)/sizeof(good_D0[0]);
//              ++expected, ++actual, ++i /*DEBUG*/) {
//             std::cout << boost::format("expected[%3d]=%12g, actual[%3d]=%12g\n")
//                 % i % *expected % i % *actual; /*DEBUG*/
//             BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
//         }
//     }
//         BOOST_CHECK_EQUAL_COLLECTIONS(
//             good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
//             w->D[0] + w->ku, w->D[0] + w->storagesize - w->kl);
//
//         /* Check w->D[0] application against multiple vectors */
//      const int nrhs = 2;
//      double vector[] = { 1, 2, 3, 4,
//                          5, 6, 7, 8 };
//      const double good_result[] = { 1, 2, 3, 4,
//                                     5, 6, 7, 8 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_apply(0, nrhs, vector, ldb, w);
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
//          vector, vector + sizeof(vector)/sizeof(vector[0]));
//  }

//  {
//      /* Check w->D[1], the first derivative matrix, against known good:
//       *  -1   1   0   0
//       *   0  -1   1   0
//       *   0   0  -1   1
//       *   0   0  -1   1
//       * Known good is in general banded matrix column-major order.
//       */
//      const double good_D1[] = { /*DK*/  -1,     0,
//                                     1,  -1,     0,
//                                     1,  -1,    -1,
//                                     1,   1  /*DK*/ };
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
//          w->D[1] + w->ku, w->D[1] + w->storagesize - w->kl);

//      /* Check w->D[0] application against multiple vectors */
//      const int nrhs = 2;
//      double vector[] = { 1, 3, 2, 4,
//                          7, 6, 5, 8 };
//      const double good_result[] = {  2, -1, 2, 2,
//                                     -1, -1, 3, 3 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_apply(1, nrhs, vector, ldb, w);
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_result, good_result + sizeof(good_result)/sizeof(good_result[0]),
//          vector, vector + sizeof(vector)/sizeof(vector[0]));
//  }

//  {
//      /* Check w->D[2], the second derivative matrix, against zero result.
//       */
//      const double good_D2[] = { /*DK*/  0,    0,
//                                     0,  0,    0,
//                                     0,  0,    0,
//                                     0,  0  /*DK*/ };
//      BOOST_CHECK_EQUAL_COLLECTIONS(
//          good_D2, good_D2 + sizeof(good_D2)/sizeof(good_D2[0]),
//          w->D[2] + w->ku, w->D[2] + w->storagesize - w->kl);
//  }

//  suzerain_bspline_operator_lu_workspace *luw
//      = suzerain_bspline_operator_lu_alloc(w);

//  /* Form 2*D[0] - 3*D[1] operator in LU-ready banded storage.  Answer is
//   *   5   -3    0     0
//   *   0    5   -3     0
//   *   0    0    5    -3
//   *   0    0    3    -1
//   * which, in LU-form where L has ones on the main diagonal, is
//   *   5   -3    0     0
//   *   0    5   -3     0
//   *   0    0    5    -3
//   *   0    0    0.6   0.8
//   * The pivot matrix is eye(4).  Check it in octave using [l,u,p] = lu(A).
//   * Known good is in general banded matrix column-major order with
//   * additional superdiagonal to allow for LU factorization fill-in.
//   */
//  const double good_A0[] = { /*DK*/   /*DK*/   5,     0,
//                             /*DK*/0,    -3,   5,     0,
//                                   0,    -3,   5,     0.6,
//                                   0,    -3,   0.8  /*DK*/ };
//  const double coeff[] = { 2.0, -3.0 };
//  suzerain_bspline_operator_lu_form(
//      sizeof(coeff)/sizeof(coeff[0]), coeff, w, luw);
//  {
//      // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
//      const double *expected, *actual;
//      for (expected = good_A0, actual = luw->A + luw->ku;
//           expected < good_A0 + sizeof(good_A0)/sizeof(good_A0[0]);
//           ++expected, ++actual) {
//          BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
//      }
//  }

//  /* Check that multiple rhs solution works for operator found just above */
//  {
//      const int nrhs = 2;
//      double vector[] = { 1,  2, 3, 4,
//                         -4, -1, 1, 3};
//      const double good_result[] = {  1.25, 1.75, 2.25, 2.75,
//                                     -0.2,  1,    2,    3 };
//      const int ldb = sizeof(vector)/(sizeof(vector[0]))/nrhs;
//      suzerain_bspline_operator_lu_solve(nrhs, vector, ldb, luw);
//      {
//          // Coarsely emulate BOOST_CHECK_EQUAL_COLLECTIONS with tolerance
//          const double *expected, *actual;
//          for (expected = good_result, actual = vector;
//               expected < good_result + sizeof(good_result)/sizeof(good_result[0]);
//               ++expected, ++actual) {
//              BOOST_CHECK_CLOSE(*expected, *actual, 1.0e-12);
//          }
//      }
//  }

//  suzerain_bspline_operator_lu_free(luw);
//    suzerain_bspline_operator_free(w);
//}

// Polynomial test helpers
typedef struct { int n; double c[]; } poly_params; // Flexible array
double poly_f(double x, void *params)
{
    poly_params *p = (poly_params *) params;
    return gsl_poly_eval(p->c, p->n, x);
}
void poly_params_differentiate(poly_params *params)
{
    for (int i = 1; i < params->n; ++i) params->c[i-1] = params->c[i] *i;
    params->c[params->n-1] = 0;
}

// Sanity check the polynomial test helpers
BOOST_AUTO_TEST_CASE( gsl_poly_eval_and_deriv )
{
    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + 3*sizeof(double));
    p->n    = 3;
    p->c[0] = 1.1; // Constant
    p->c[1] = 2.2; // Linear
    p->c[2] = 3.3; // Quadratic
    suzerain_function f = {poly_f, p};

    const double   dc[] = { 2.2, 6.6, 0.0 }; // (d/dx)(p->c)
    const double  ddc[] = { 6.6, 0.0, 0.0 }; // (d^2/dx^2)(p->c)
    const double dddc[] = { 0.0, 0.0, 0.0 }; // (d^3/dx^3)(p->c)

    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,1.0),  6.6, 0.001);
    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,2.0), 18.7, 0.001);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dc, dc + p->n);

    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, ddc, ddc + p->n);

    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,1.0), 6.6, 0.001);
    BOOST_CHECK_CLOSE(SUZERAIN_FN_EVAL(&f,2.0), 6.6, 0.001);

    // Differentiate a constant
    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dddc, dddc + p->n);

    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,1.0), 0.0);
    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,2.0), 0.0);

    // Differentiate zero
    poly_params_differentiate(p);
    BOOST_CHECK_EQUAL_COLLECTIONS(p->c, p->c + p->n, dddc, dddc + p->n);

    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,1.0), 0.0);
    BOOST_CHECK_EQUAL(SUZERAIN_FN_EVAL(&f,2.0), 0.0);

    free(p);
}

BOOST_AUTO_TEST_CASE( functioncoefficients )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 2;

    poly_params *p = (poly_params *)
                      malloc(sizeof(poly_params) + 3*sizeof(double));
    p->n = 3;
    suzerain_function f = {poly_f, p};

    suzerain_bspline_operator_workspace *w
        = suzerain_bspline_operator_alloc(order, nderiv, nbreak, breakpoints,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);
    const int ncoeff = suzerain_bspline_operator_ncoefficients(w);

    // Form the mass matrix M
    suzerain_bspline_operator_lu_workspace * const mass
        = suzerain_bspline_operator_lu_alloc(w);
    const double d_one = 1.0;
    int    i_one = 1;
    suzerain_bspline_operator_lu_form(i_one, &d_one, w, mass);

    {
        const int derivative = 1;

        p->c[0] = 1.2; // Constant
        p->c[1] = 3.4; // Linear
        p->c[2] = 0.0; // Quadratic

        // Compute the right hand side coefficients for M x = b
        double * coefficient = (double *) malloc(ncoeff * sizeof(double));
        suzerain_bspline_operator_functioncoefficient_rhs(&f, coefficient, w);

        // Solve for function coefficients using the mass matrix
        suzerain_bspline_operator_lu_solve(1, coefficient, ncoeff, mass);

        // Take the n-th derivative of the coefficients using M x' = D x
        suzerain_bspline_operator_apply(derivative, 1, coefficient, ncoeff, w);
        suzerain_bspline_operator_lu_solve(1, coefficient, ncoeff, mass);

        // Ensure we recover the leading order, scaled monomial coefficients
        for (int i = 0; i < ncoeff; ++i) {
            BOOST_CHECK_CLOSE(1.0 * p->c[1], coefficient[i], 1e-12);
        }

        free(coefficient);
    }

//  {
//      const int derivative = 1;

//      p->c[0] = 1.2; // Constant
//      p->c[1] = 3.4; // Linear
//      p->c[2] = 0.0; // Quadratic

//      // Compute the right hand side coefficients for M x = b
//      double * coefficient = (double *) malloc(ncoeff * sizeof(double));
//      suzerain_bspline_operator_functioncoefficient_rhs(&f, coefficient, w);

//      // Solve for function coefficients using the mass matrix
//      suzerain_bspline_operator_lu_solve(1, coefficient, ncoeff, mass);

//      // Take the n-th derivative of the coefficients using M x' = D x
//      suzerain_bspline_operator_apply(derivative, 1, coefficient, ncoeff, w);
//      suzerain_bspline_operator_lu_solve(1, coefficient, ncoeff, mass);

//      // Ensure we recover the leading order, scaled monomial coefficients
//      for (int i = 0; i < ncoeff; ++i) {
//          BOOST_CHECK_CLOSE(1.0 * p->c[1], coefficient[i], 1e-12);
//      }

//      free(coefficient);
//  }

    suzerain_bspline_operator_lu_free(mass);
    suzerain_bspline_operator_free(w);
    free(p);
}
