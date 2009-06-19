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

#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <suzerain/underling_error.h>

underling_error_handler_t * underling_error_handler = NULL;

static void no_error_handler (const char *reason, const char *file, int line, int underling_errno);

void
underling_error (const char * reason, const char * file, int line, int underling_errno)
{
  if (underling_error_handler)
    {
      (*underling_error_handler) (reason, file, line, underling_errno);
      return ;
    }

  underling_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}

underling_error_handler_t *
underling_set_error_handler (underling_error_handler_t * new_handler)
{
  underling_error_handler_t * previous_handler = underling_error_handler;
  underling_error_handler = new_handler;
  return previous_handler;
}


underling_error_handler_t *
underling_set_error_handler_off (void)
{
  underling_error_handler_t * previous_handler = underling_error_handler;
  underling_error_handler = no_error_handler;
  return previous_handler;
}

static void
no_error_handler (const char *reason, const char *file, int line, int underling_errno)
{
  /* do nothing */
  reason = 0;
  file = 0;
  line = 0;
  underling_errno = 0;
  return;
}

FILE * underling_stream = NULL ;
underling_stream_handler_t * underling_stream_handler = NULL;

void
underling_stream_printf(const char *label, const char *file, int line,
                   const char *reason)
{
  if (underling_stream == NULL)
    {
      underling_stream = stderr;
    }
  if (underling_stream_handler)
    {
      (*underling_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (underling_stream, "underling: %s:%d: %s: %s\n", file, line, label, reason);

}

underling_stream_handler_t *
underling_set_stream_handler (underling_stream_handler_t * new_handler)
{
  underling_stream_handler_t * previous_handler = underling_stream_handler;
  underling_stream_handler = new_handler;
  return previous_handler;
}

FILE *
underling_set_stream (FILE * new_stream)
{
  FILE * previous_stream;
  if (underling_stream == NULL) {
    underling_stream = stderr;
  }
  previous_stream = underling_stream;
  underling_stream = new_stream;
  return previous_stream;
}

const char *
underling_strerror(const int underling_errno)
{
  switch (underling_errno)
    {
    case UNDERLING_SUCCESS:
      return "success" ;
    case UNDERLING_FAILURE:
      return "failure" ;
    case UNDERLING_EDOM:
      return "input domain error" ;
    case UNDERLING_ERANGE:
      return "output range error" ;
    case UNDERLING_EFAULT:
      return "invalid pointer" ;
    case UNDERLING_EINVAL:
      return "invalid argument supplied by user" ;
    case UNDERLING_EFAILED:
      return "generic failure" ;
    case UNDERLING_ESANITY:
      return "sanity check failed - shouldn't happen" ;
    case UNDERLING_ENOMEM:
      return "malloc failed" ;
    case UNDERLING_EBADFUNC:
      return "problem with user-supplied function";
    case UNDERLING_EZERODIV:
      return "tried to divide by zero" ;
    default:
      return "unknown error code" ;
    }
}
