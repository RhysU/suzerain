/*--------------------------------------------------------------------------
 *
 * Functionality adopted from the GNU Scientific Library.
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * Copyright (C) 2010, 2011 Rhys Ulerich
 * Copyright (C) 2011 The PECOS Development Team
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

/** @file
 * @copydoc error.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/error.h>

#include <suzerain/common.h>

suzerain_error_handler_t * suzerain_error_handler = NULL;

static void
no_error_handler(const char *reason,
                 const char *file,
                 int line,
                 int suzerain_errno);

void
suzerain_error(const char * reason,
               const char * file,
               int line,
               int suzerain_errno)
{
    if (suzerain_error_handler) {
        (*suzerain_error_handler) (reason, file, line, suzerain_errno);
        return ;
    }

    suzerain_stream_printf ("ERROR", file, line, reason);

    fflush (stdout);
    fprintf (stderr, "Default suzerain error handler invoked.\n");
    fflush (stderr);

    abort ();
}

suzerain_error_handler_t *
suzerain_set_error_handler(suzerain_error_handler_t * new_handler)
{
    suzerain_error_handler_t * previous_handler = suzerain_error_handler;
    suzerain_error_handler = new_handler;
    return previous_handler;
}


suzerain_error_handler_t *
suzerain_set_error_handler_off(void)
{
    suzerain_error_handler_t * previous_handler = suzerain_error_handler;
    suzerain_error_handler = no_error_handler;
    return previous_handler;
}

#pragma warning(push,disable:869)
static void
no_error_handler(const char *reason,
                 const char *file,
                 int line,
                 int suzerain_errno)
{
    SUZERAIN_UNUSED(reason);
    SUZERAIN_UNUSED(file);
    SUZERAIN_UNUSED(line);
    SUZERAIN_UNUSED(suzerain_errno);
    return; /* do nothing */
}
#pragma warning(pop)

FILE * suzerain_stream = NULL ;
suzerain_stream_handler_t * suzerain_stream_handler = NULL;

void
suzerain_stream_printf(const char *label,
                       const char *file,
                       int line,
                       const char *reason)
{
    if (suzerain_stream == NULL) {
        suzerain_stream = stderr;
    }
    if (suzerain_stream_handler) {
        (*suzerain_stream_handler) (label, file, line, reason);
        return;
    }
    fprintf(suzerain_stream,
            "suzerain: %s:%d: %s: %s\n", file, line, label, reason);

}

suzerain_stream_handler_t *
suzerain_set_stream_handler(suzerain_stream_handler_t * new_handler)
{
    suzerain_stream_handler_t * previous_handler = suzerain_stream_handler;
    suzerain_stream_handler = new_handler;
    return previous_handler;
}

FILE *
suzerain_set_stream(FILE * new_stream)
{
    FILE * previous_stream;
    if (suzerain_stream == NULL) {
        suzerain_stream = stderr;
    }
    previous_stream = suzerain_stream;
    suzerain_stream = new_stream;
    return previous_stream;
}

const char *
suzerain_strerror(const int suzerain_errno)
{
    switch (suzerain_errno) {
    case SUZERAIN_FAILURE:
        return "failure" ;
    case SUZERAIN_CONTINUE:
        return "iteration has not converged";
    case SUZERAIN_SUCCESS:
        return "success" ;
    case SUZERAIN_EDOM:
        return "input domain error" ;
    case SUZERAIN_ERANGE:
        return "output range error" ;
    case SUZERAIN_EFAULT:
        return "invalid pointer" ;
    case SUZERAIN_EINVAL:
        return "invalid argument supplied by user" ;
    case SUZERAIN_EFAILED:
        return "generic failure" ;
    case SUZERAIN_ESANITY:
        return "sanity check failed - shouldn't happen" ;
    case SUZERAIN_ENOMEM:
        return "malloc failed" ;
    case SUZERAIN_EBADFUNC:
        return "problem with user-supplied function";
    case SUZERAIN_EMAXITER:
        return "exceeded max number of iterations";
    case SUZERAIN_EZERODIV:
        return "tried to divide by zero" ;
    case SUZERAIN_EROUND:
        return "failed because of roundoff error";
    case SUZERAIN_EBADLEN:
        return "matrix or vector lengths are not conformant";
    case SUZERAIN_ESING:
        return "apparent singularity detected";
    case SUZERAIN_EDIVERGE:
        return "integral or series is divergent";
    case SUZERAIN_EUNIMPL:
        return "requested feature not (yet) implemented";
    default:
        return "unknown error code" ;
    }
}
