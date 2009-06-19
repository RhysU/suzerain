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

#ifndef PECOS_UNDERLING_ERRNO_H__
#define PECOS_UNDERLING_ERRNO_H__

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
  UNDERLING_SUCCESS  = 0,
  UNDERLING_FAILURE  = -1,
  UNDERLING_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  UNDERLING_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  UNDERLING_EFAULT   = 3,   /* invalid pointer */
  UNDERLING_EINVAL   = 4,   /* invalid argument supplied by user */
  UNDERLING_EFAILED  = 5,   /* generic failure */
  UNDERLING_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  UNDERLING_ENOMEM   = 8,   /* malloc failed */
  UNDERLING_EBADFUNC = 9,   /* problem with user-supplied function */
  UNDERLING_EZERODIV = 12   /* tried to divide by zero */
};

void underling_error(const char * reason,
                     const char * file,
                     int line,
                     int underling_errno);

void underling_stream_printf(const char *label,
                             const char *file,
                             int line,
                             const char *reason);

const char * underling_strerror(const int underling_errno);

typedef void underling_error_handler_t(const char * reason,
                                       const char * file,
                                       int line,
                                       int underling_errno);

typedef void underling_stream_handler_t(const char * label,
                                        const char * file,
                                        int line,
                                        const char * reason);

underling_error_handler_t *
underling_set_error_handler(underling_error_handler_t * new_handler);

underling_error_handler_t *
underling_set_error_handler_off(void);

underling_stream_handler_t *
underling_set_stream_handler(underling_stream_handler_t * new_handler);

FILE * underling_set_stream(FILE * new_stream);

/* UNDERLING_ERROR: call the error handler, and return the error code */

#define UNDERLING_ERROR(reason, underling_errno) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return underling_errno ; \
       } while (0)

/* UNDERLING_ERROR_VAL: call the error handler, and return the given value */

#define UNDERLING_ERROR_VAL(reason, underling_errno, value) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return value ; \
       } while (0)

/* UNDERLING_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define UNDERLING_ERROR_VOID(reason, underling_errno) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return ; \
       } while (0)

/* UNDERLING_ERROR_NULL suitable for out-of-memory conditions */

#define UNDERLING_ERROR_NULL(reason, underling_errno) \
        UNDERLING_ERROR_VAL(reason, underling_errno, 0)

__END_DECLS

#endif /* PECOS_UNDERLING_ERRNO_H__ */
