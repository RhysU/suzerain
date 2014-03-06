/*--------------------------------------------------------------------------
 *
 * Code copyright and licensing details follow on routine-specific basis.
 *
 *--------------------------------------------------------------------------
 */

/** @file
 * @copydoc pre_gsl.h
 */

#include <suzerain/pre_gsl.h>

#include <gsl/gsl_ieee_utils.h>
#include <stdio.h>
#include <unistd.h>

#include <suzerain/common.h>

/*
 * License identical to gsl_ieee_env_setup() from GSL.
 * Copyright (c) 2012-2014 Rhys Ulerich
 */
void
mpi_gsl_ieee_env_setup(const int rank)
{
    // Used for write(2)-based error messages
    static const char newline[] = "\n";

    // Flush stdout, stderr
    if (fflush(stdout))
        perror("mpi_gsl_ieee_env_setup fflush(stdout) before redirect");
    if (fflush(stderr))
        perror("mpi_gsl_ieee_env_setup fflush(stderr) before redirect");

    // Save stdout, stderr so we may restore them later
    int stdout_copy, stderr_copy;
    if ((stdout_copy = dup(fileno(stdout))) < 0)
        perror("mpi_gsl_ieee_env_setup error duplicating stdout");
    if ((stderr_copy = dup(fileno(stderr))) < 0)
        perror("mpi_gsl_ieee_env_setup error duplicating stderr");

    // On non-root processes redirect stdout, stderr to /dev/null
    if (rank) {
        if (!freopen("/dev/null", "a", stdout))
            perror("mpi_gsl_ieee_env_setup redirecting stdout");
        if (!freopen("/dev/null", "a", stderr))
            perror("mpi_gsl_ieee_env_setup redirecting stderr");
    }

    // Invoke gsl_ieee_env_setup on all ranks.
    gsl_ieee_env_setup();

    // Flush stdout, stderr again.
    // Error messages sent to stderr_copy on a best-effort basis
    if (fflush(stdout)) {
        static const char pre[] = "mpi_gsl_ieee_env_setup fflush(stdout) after redirect: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }
    if (fflush(stderr)) {
        static const char pre[] = "mpi_gsl_ieee_env_setup fflush(stderr) after redirect: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }

    // Restore stdout, stderr
    // Error messages sent to stderr_copy on a best-effort basis
    if (dup2(stdout_copy, fileno(stdout)) < 0) {
        static const char pre[] = "mpi_gsl_ieee_env_setup reopening stdout: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }
    if (dup2(stderr_copy, fileno(stderr)) < 0) {
        static const char pre[] = "mpi_gsl_ieee_env_setup reopening stderr: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }

    // Close saved versions of stdout, stderr
    if (close(stdout_copy))
        perror("mpi_gsl_ieee_env_setup closing stdout_copy");
    if (close(stderr_copy))
        perror("mpi_gsl_ieee_env_setup closing stderr_copy");

    // Clear any errors that may have occurred on stdout, stderr
    clearerr(stdout);
    clearerr(stderr);
}

/*
 * integration/glfixed.c
 *
 * Numerical Integration by Gauss-Legendre Quadrature Formulas of high orders.
 * High-precision abscissas and weights are used.
 *
 * Original project homepage: http://www.holoborodko.com/pavel/?page_id=679
 * Original contact e-mail:   pavel@holoborodko.com
 *
 * Copyright (c)2007-2008 Pavel Holoborodko
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributors:
 * Konstantin Holoborodko - Optimization of Legendre polynomial computing.
 * Rhys Ulerich - Inclusion within GNU Scientific Library.
 *
 */

// gsl_integration_glfixed in GSL 1.14 but
// gsl_integration_glfixed_point is 1.14+ so build it atop 1.14's public API
// Source code lifted verbatim from the GSL
// FIXME: Remove this logic once GSL 1.15 becomes widespread
#if    (!defined GSL_MAJOR_VERSION                     ) \
    || (GSL_MAJOR_VERSION < 2 && GSL_MINOR_VERSION < 15)
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

int
gsl_integration_glfixed_point (
        double a,
        double b,
        size_t i,
        double *xi,
        double *wi,
        const gsl_integration_glfixed_table * t)
{
    const double A = (b - a) / 2;  /* Length of [a,b] */
    const double B = (a + b) / 2;  /* Midpoint of [a,b] */

    if (i >= t->n) {
        GSL_ERROR ("i must be less than t->n", GSL_EINVAL);
    }

    /* See comments above gsl_integration_glfixed for struct's x, w layout. */
    /* Simply unpack that layout into a sorted set of points, weights. */
    if (GSL_IS_ODD(t->n)) {
        const int k = ((int) i) - ((int) t->n) / 2;
        const int s = k < 0 ? -1 : 1;

        *xi = B + s*A*t->x[s*k];
        *wi =       A*t->w[s*k];
    } else if (/* GSL_IS_EVEN(t->n) && */ i < t->n / 2) {
        i = (t->n / 2) - 1 - i;
        *xi = B - A*t->x[i];
        *wi =     A*t->w[i];
    } else /* GSL_IS_EVEN(t->n) && i >= n / 2 */ {
        i  -= t->n / 2;
        *xi = B + A*t->x[i];
        *wi =     A*t->w[i];
    }

    return GSL_SUCCESS;
}


/* bspline/bspline.c
 *
 * Copyright (C) 2006, 2007, 2008, 2009 Patrick Alken
 * Copyright (C) 2008 Rhys Ulerich
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// gsl_bspline_knots_greville is 1.15+ so build it atop 1.15's public API
// FIXME: Remove this logic once GSL 1.16 becomes widespread
#if    (!defined GSL_MAJOR_VERSION                     ) \
    || (GSL_MAJOR_VERSION < 2 && GSL_MINOR_VERSION < 16)
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

int
gsl_bspline_knots_greville (const gsl_vector *abscissae,
                            gsl_bspline_workspace *w,
                            double *abserr)
{
  int s;

  /* Check incoming arguments satisfy mandatory algorithmic assumptions */
  if (w->k < 2)
    GSL_ERROR ("w->k must be at least 2", GSL_EINVAL);
  else if (abscissae->size < 2)
    GSL_ERROR ("abscissae->size must be at least 2", GSL_EINVAL);
  else if (w->nbreak != abscissae->size - w->k + 2)
    GSL_ERROR ("w->nbreak must equal abscissae->size - w->k + 2", GSL_EINVAL);

  if (w->nbreak == 2)
    {
      /* No flexibility in abscissae values possible in this degenerate case */
      s = gsl_bspline_knots_uniform (
              gsl_vector_get (abscissae, 0),
              gsl_vector_get (abscissae, abscissae->size - 1), w);
    }
  else
    {
      double * storage;
      gsl_matrix_view A;
      gsl_vector_view tau, b, x, r;
      size_t i, j;

      /* Constants derived from the B-spline workspace and abscissae details */
      const size_t km2    = w->k - 2;
      const size_t M      = abscissae->size - 2;
      const size_t N      = w->nbreak - 2;
      const double invkm1 = 1.0 / w->km1;

      /* Allocate working storage and prepare multiple, zero-filled views */
      storage = (double *) calloc (M*N + 2*N + 2*M, sizeof (double));
      if (storage == 0)
        GSL_ERROR ("failed to allocate working storage", GSL_ENOMEM);
      A   = gsl_matrix_view_array (storage, M, N);
      tau = gsl_vector_view_array (storage + M*N,             N);
      b   = gsl_vector_view_array (storage + M*N + N,         M);
      x   = gsl_vector_view_array (storage + M*N + N + M,     N);
      r   = gsl_vector_view_array (storage + M*N + N + M + N, M);

      /* Build matrix from interior breakpoints to interior Greville abscissae.
       * For example, when w->k = 4 and w->nbreak = 7 the matrix is
       *   [   1,      0,      0,      0,      0;
       *     2/3,    1/3,      0,      0,      0;
       *     1/3,    1/3,    1/3,      0,      0;
       *       0,    1/3,    1/3,    1/3,      0;
       *       0,      0,    1/3,    1/3,    1/3;
       *       0,      0,      0,    1/3,    2/3;
       *       0,      0,      0,      0,      1  ]
       * but only center formed as first/last breakpoint is known.
       */
      for (j = 0; j < N; ++j)
        for (i = 0; i <= km2; ++i)
          gsl_matrix_set (&A.matrix, i+j, j, invkm1);

      /* Copy interior collocation points from abscissae into b */
      for (i = 0; i < M; ++i)
        gsl_vector_set (&b.vector, i, gsl_vector_get (abscissae, i+1));

      /* Adjust b to account for constraint columns not stored in A */
      for (i = 0; i < km2; ++i)
        {
          double * const v = gsl_vector_ptr (&b.vector, i);
          *v -= (1 - (i+1)*invkm1) * gsl_vector_get (abscissae, 0);
        }
      for (i = 0; i < km2; ++i)
        {
          double * const v = gsl_vector_ptr (&b.vector, M - km2 + i);
          *v -= (i+1)*invkm1 * gsl_vector_get (abscissae, abscissae->size - 1);
        }

      /* Perform linear least squares to determine interior breakpoints */
      s =  gsl_linalg_QR_decomp (&A.matrix, &tau.vector)
        || gsl_linalg_QR_lssolve (&A.matrix, &tau.vector,
                                  &b.vector, &x.vector, &r.vector);
      if (s)
        {
          free (storage);
          return s;
        }

      /* "Expand" solution x by adding known first and last breakpoints. */
      x = gsl_vector_view_array_with_stride (
          gsl_vector_ptr (&x.vector, 0) - x.vector.stride,
          x.vector.stride, x.vector.size + 2);
      gsl_vector_set (&x.vector, 0, gsl_vector_get (abscissae, 0));
      gsl_vector_set (&x.vector, x.vector.size - 1,
                      gsl_vector_get (abscissae, abscissae->size - 1));

      /* Finally, initialize workspace knots using the now-known breakpoints */
      s = gsl_bspline_knots (&x.vector, w);
      free (storage);
    }

  /* Sum absolute errors in the resulting vs requested interior abscissae */
  /* Provided as a fit quality metric which may be monitored by callers */
  if (!s && abserr)
    {
      size_t i;
      *abserr = 0;
      for (i = 1; i < abscissae->size - 1; ++i)
        *abserr += fabs (   gsl_bspline_greville_abscissa (i, w)
                          - gsl_vector_get (abscissae, i) );
    }

  return s;
}
#endif

#endif
