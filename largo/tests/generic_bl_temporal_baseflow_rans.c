//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// largo 0.0.2: largo - slow growth terms for turbulence simulations
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2011-2014 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <largo/largo.h>

// Include FCTX and silence useless warnings
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:981)
#endif
#include "fct.h"
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif


// See FCTX manual at http://fctx.wildbearsoftware.com/
FCT_BGN()
{
    FCT_FIXTURE_SUITE_BGN(esio_file)
    {
        // Code which must run before every test case
        FCT_SETUP_BGN()
        {
            // Nothing
        }
        FCT_SETUP_END();

        // Code which must run after every test case
        FCT_TEARDOWN_BGN()
        {
            // Nothing
        }
        FCT_TEARDOWN_END();

        // A test case
        FCT_TEST_BGN(something)
        {
          const char ransmodel[] = "self-similar";
          const int ntvar        = 2;
          const int ns           = 2;
          const int neq          = 5 + ns + ntvar;
          const int nvar         = neq;
          const int nvar_base    = 5 + ns;

          largo_workspace * work;
          double y       = 1.0/ 10.0;
          double grDelta = 5.0/100.0;

          double mean  [] = \
          {     1.0/ 100.0,  \
               45.0/  10.0,  \
                1.0/1000.0,  \
                5.0/ 100.0,  \
            41200.0       ,  \
                2.0/1000.0,  \
                1.0/1000.0,  \
                5.0/ 100.0,  \
                3.0/ 100.0 
          };

          double dmean [] = \
          {      1.0/ 10.0,  \
                45.0      ,  \
                 1.0/100.0,  \
                 5.0/ 10.0,  \
            412000.0      ,  \
                 2.0/100.0,  \
                 1.0/100.0,  \
                 5.0/ 10.0,  \
                 3.0/ 10.0 
          };

          double grDA   [] = \
          {    2.0/   100.0, \
               1.0/    10.0, \
               1.0/   100.0, \
               3.0/   100.0, \
               4.0/   100.0, \
               1.0/  1000.0, \
               2.0/  1000.0, \
               4.0/   100.0, \
               2.0/   100.0 
          };

          double grDArms[] = \
          {    1.0/  1000.0, \
               5.0/  1000.0, \
               5.0/ 10000.0, \
               2.0/  1000.0, \
               1.0/  1000.0, \
               1.0/100000.0, \
               2.0/100000.0, \
               0.0         , \
               0.0         
          };

          double base   [] = \
          {    5.0/  1000.0, \
              25.0/    10.0, \
               5.0/ 10000.0, \
               2.0/   100.0, \
           20000.0         , \
               1.0/  1000.0, \
               5.0/ 10000.0
           };


          double dybase [] = \
          {    5.0/   100.0, \
              25.0         , \
               5.0/  1000.0, \
               2.0/    10.0, \
          200000.0         , \
               1.0/   100.0, \
               5.0/  1000.0
          };


          double dtbase [] = \
          {    2.0/  1000.0, \
              15.0/    10.0, \
               2.0/ 10000.0, \
               1.0/   100.0, \
           10000.0         , \
               5.0/ 10000.0, \
               2.0/ 10000.0
           };

          double dxbase  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcbase []=\
           {   1.0/  1000.0, \
               5.0/    10.0, \
               1.0/ 10000.0, \
               5.0/  1000.0, \
            5000.0         , \
               2.0/ 10000.0, \
               1.0/ 10000.0, \
               0.0         , \
               0.0         
           };

          double srcmean [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcmean_good [] = \
          { -  17.0/  20000.0,     \
            -  11.0/     10.0,     \
            -   1.0/  12500.0,     \
            -  11.0/   2500.0,     \
            -4788.0          ,     \
            - 251.0/1000000.0,     \
            -  19.0/ 250000.0,     \
                1.0/   2000.0,     \
                9.0/  10000.0  
          };

          const double tolerance = 1.0E-12;

          const char model[] = "bl_temporal";


          // Allocate generic workspace;
          largo_allocate (&work, model, neq, ns, ntvar, ransmodel);

          // Init growth rates
          largo_init_rans     (work, grDelta, grDA, grDArms);

          // Compute prestep values;
          largo_prestep_baseflow  (work,   base,  dybase,
                                          dtbase, dxbase, srcbase);
          largo_prestep_setamean_rans (work, y, mean, dmean);

          // Compute sources using wrapper method;
          largo_seta_rans (work, 0.0, 1.0, &srcmean[1-1]);

          // Check sources
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcmean[1-1]       , srcmean_good[1-1]);
          fct_chk_eq_dbl(srcmean[2-1]/10.0  , srcmean_good[2-1]/10.0);
          fct_chk_eq_dbl(srcmean[3-1]       , srcmean_good[3-1]);
          fct_chk_eq_dbl(srcmean[4-1]       , srcmean_good[4-1]);
          fct_chk_eq_dbl(srcmean[5-1]/10000., srcmean_good[5-1]/10000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcmean[5+is]    , srcmean_good[5+is]);
          }
          for (unsigned int it=0; it < ntvar; ++it)
          {
            fct_chk_eq_dbl(srcmean[5+ns+it] , srcmean_good[5+ns+it]);
          }

          // Deallocate workspace
          largo_deallocate (&work);

        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
