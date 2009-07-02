/*
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Adapted from the GNU Scientific Library
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

#ifndef PECOS_SUZERAIN_ERRNO_H__
#define PECOS_SUZERAIN_ERRNO_H__

#include <stdio.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum {
  SUZERAIN_SUCCESS  = 0,
  SUZERAIN_FAILURE  = -1,
  SUZERAIN_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  SUZERAIN_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  SUZERAIN_EFAULT   = 3,   /* invalid pointer */
  SUZERAIN_EINVAL   = 4,   /* invalid argument supplied by user */
  SUZERAIN_EFAILED  = 5,   /* generic failure */
  SUZERAIN_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  SUZERAIN_ENOMEM   = 8,   /* malloc failed */
  SUZERAIN_EBADFUNC = 9,   /* problem with user-supplied function */
  SUZERAIN_EZERODIV = 12   /* tried to divide by zero */
};

void suzerain_error(const char * reason,
                     const char * file,
                     int line,
                     int suzerain_errno);

void suzerain_stream_printf(const char *label,
                             const char *file,
                             int line,
                             const char *reason);

const char * suzerain_strerror(const int suzerain_errno);

typedef void suzerain_error_handler_t(const char * reason,
                                       const char * file,
                                       int line,
                                       int suzerain_errno);

typedef void suzerain_stream_handler_t(const char * label,
                                        const char * file,
                                        int line,
                                        const char * reason);

suzerain_error_handler_t *
suzerain_set_error_handler(suzerain_error_handler_t * new_handler);

suzerain_error_handler_t *
suzerain_set_error_handler_off(void);

suzerain_stream_handler_t *
suzerain_set_stream_handler(suzerain_stream_handler_t * new_handler);

FILE * suzerain_set_stream(FILE * new_stream);

/* SUZERAIN_ERROR: call the error handler, and return the error code */

#define SUZERAIN_ERROR(reason, suzerain_errno) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return suzerain_errno ; \
       } while (0)

/* SUZERAIN_ERROR_VAL: call the error handler, and return the given value */

#define SUZERAIN_ERROR_VAL(reason, suzerain_errno, value) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return value ; \
       } while (0)

/* SUZERAIN_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define SUZERAIN_ERROR_VOID(reason, suzerain_errno) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return ; \
       } while (0)

/* SUZERAIN_ERROR_NULL suitable for out-of-memory conditions */

#define SUZERAIN_ERROR_NULL(reason, suzerain_errno) \
        SUZERAIN_ERROR_VAL(reason, suzerain_errno, 0)

__END_DECLS

#endif /* PECOS_SUZERAIN_ERRNO_H__ */
