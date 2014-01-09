// Adapted from the GNU Scientific Library
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include <underling/config.h>
#endif
#include <underling/error.h>

#define CHECK(x) errors[n].number = x ; errors[n].name = #x ; n++ ;
#define MAX_ERRS 64

int verbose = 0 ;

int
main (void)
{
  int n = 0, r = 0;

#pragma warning(push,disable:2021)
  struct {
    int number;
    const char * name;
  } errors[MAX_ERRS] ;
#pragma warning(pop)

  CHECK(UNDERLING_SUCCESS);
  CHECK(UNDERLING_EFAULT);
  CHECK(UNDERLING_EINVAL);
  CHECK(UNDERLING_EFAILED);
  CHECK(UNDERLING_ESANITY);
  CHECK(UNDERLING_ENOMEM);

  for (int i = 0 ; i < n ; i++)
    {
      if (verbose) {
        printf("%s = %d\n", errors[i].name, errors[i].number);
      }
    }

  for (int i = 0; i < n; i++)
    {
      int status = 0;
      for (int j = 0; j < n; j++)
        {
          if (j != i)
              status |= (errors[i].number == errors[j].number);
        }
      if (status)
        {
          r = 1;
          printf("Found non-distinct error value\n");
        }
    }

  for (int i = 0; i < n; i++)
    {
      int status = 0;
      int e1 = errors[i].number ;
      for (int j = 0; j < n; j++)
        {
          if (j != i)
            {
              int e2 = errors[j].number;
              status |= (underling_strerror(e1) == underling_strerror(e2)) ;
            }
        }
      if (status)
        {
          r = 1;
          printf("Found non-distinct error message\n");
        }
    }

  return r;
}
