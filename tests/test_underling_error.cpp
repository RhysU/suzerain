/*
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
#define BOOST_TEST_MODULE $Id: test_underling_error.cpp$

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <boost/test/included/unit_test.hpp>
#include <suzerain/underling_error.h>

#define CHECK(x) errors[n].number = x ; errors[n].name = #x ; n++ ;
#define MAX_ERRS 64

int verbose = 0 ;

BOOST_AUTO_TEST_CASE( main_test )
{
  int i, j, n = 0 ;

  struct {
    int number;
    const char * name;
  } errors[MAX_ERRS] ;

  CHECK(UNDERLING_SUCCESS);
  CHECK(UNDERLING_FAILURE);
  CHECK(UNDERLING_EDOM);
  CHECK(UNDERLING_ERANGE);
  CHECK(UNDERLING_EFAULT);
  CHECK(UNDERLING_EINVAL);
  CHECK(UNDERLING_EFAILED);
  CHECK(UNDERLING_ESANITY);
  CHECK(UNDERLING_ENOMEM);
  CHECK(UNDERLING_EBADFUNC);
  CHECK(UNDERLING_EZERODIV);

  for (i = 0 ; i < n ; i++)
    {
      if (verbose) printf ("%s = %d\n", errors[i].name, errors[i].number) ;
    }

  for (i = 0; i < n; i++)
    {
      int status = 0;
      for (j = 0; j < n; j++)
        {
          if (j != i)
              status |= (errors[i].number == errors[j].number);
        }

      BOOST_REQUIRE_MESSAGE(!status, "Found non-distinct error value");
    }

  for (i = 0; i < n; i++)
    {
      int status = 0;
      int e1 = errors[i].number ;
      for (j = 0; j < n; j++)
        {
          if (j != i)
            {
              int e2 = errors[j].number;
              status |= (underling_strerror(e1) == underling_strerror(e2)) ;
            }
        }
      BOOST_REQUIRE_MESSAGE(!status, "Found non-distinct error message");
    }
}
