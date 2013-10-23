//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// largo 0.0.1: largo - slow growth terms for turbulence simulations
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2011, 2012, 2013 The PECOS Development Team
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
          const int neq = 7;
          const int ns  = 2;

          void * generic_workspace;
          double y       = 1.0/ 10.0;
          double grDelta = 5.0/100.0;

          double field  [] = \
          {    11.0/ 1000.0,  \
              485.0/  100.0,  \
                2.0/   10.0,  \
                3.0/   10.0,  \
            41500.0        ,  \
               22.0/10000.0,  \
               11.0/10000.0
          };

          double mean  [] = \
          {     1.0/ 100.0,  \
               45.0/  10.0,  \
                1.0/1000.0,  \
                5.0/ 100.0,  \
            41200.0       ,  \
                2.0/1000.0,  \
                1.0/1000.0
          };

          double dmean [] = \
          {      1.0/  5.0,  \
                45.0      ,  \
                 1.0/100.0,  \
                 5.0/ 10.0,  \
            412000.0      ,  \
                 2.0/100.0,  \
                 1.0/100.0
          };

          double mean_rqq [] = \
          {           21.0/  2000.0,  \
                    8505.0/     4.0,  \
                      21.0/200000.0,  \
                      21.0/    80.0,  \
            178231200000.0         ,  \
                      21.0/ 50000.0,  \
                      21.0/200000.0
          };

          double dmean_rqq [] = \
          {            11.0/    50.0, \
                     2025.0         , \
                        1.0/ 10000.0, \
                        1.0/     4.0, \
             169744000000.0         , \
                        1.0/  2500.0, \
                        1.0/ 10000.0
          };

          double rms [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double drms [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double grDA [] = \
          { 2.0/  100.0, \
            1.0/   10.0, \
            1.0/  100.0, \
            3.0/  100.0, \
            4.0/  100.0, \
            1.0/ 1000.0, \
            2.0/ 1000.0
          };

          double grDArms  [] = \
          { 1.0/  1000.0, \
            5.0/  1000.0, \
            5.0/ 10000.0, \
            2.0/  1000.0, \
            1.0/  1000.0, \
            1.0/100000.0, \
            2.0/100000.0
          };

          double base [] = \
          {     5.0/ 1000.0, \
               25.0/   10.0, \
                5.0/10000.0, \
                2.0/  100.0, \
            20000.0        , \
                1.0/ 1000.0, \
                5.0/10000.0
          };

          double dybase  [] = \
          {      5.0/ 100.0, \
                25.0       , \
                 5.0/1000.0, \
                 2.0/  10.0, \
            200000.0       , \
                 1.0/ 100.0, \
                 5.0/1000.0
          };

          double dtbase [] = \
          {     2.0/ 1000.0, \
               15.0/   10.0, \
                2.0/10000.0, \
                1.0/  100.0, \
            10000.0        , \
                5.0/10000.0, \
                2.0/10000.0
          };

          double dxbase [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcbase  [] = \
          {    1.0/ 1000.0, \
               5.0/   10.0, \
               1.0/10000.0, \
               5.0/ 1000.0, \
            5000.0        , \
               2.0/10000.0, \
               1.0/10000.0
          };

          double srcmean [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcfull [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcmean_good [] = \
          {     -27.0/ 20000.0,    \
               -713.0/   400.0,    \
                -37.0/200000.0,    \
               -271.0/ 20000.0,    \
             -11670.0         ,    \
                -57.0/100000.0,    \
                -37.0/200000.0     \
          };

          double srcfull_good  [] = \
          {     -297.0/  200000.0, \
               -7787.0/    4000.0, \
             -543089.0/20000000.0, \
               -4847.0/  100000.0, \
             -307937.0/      25.0, \
                -627.0/ 1000000.0, \
                -407.0/ 2000000.0  \
          };

          const double tolerance = 1.0E-7;

          const char model[] = "bl_temporal_consistent";

          // Allocate generic workspace;
          largo_allocate (&generic_workspace, model, neq, ns, 0, "dns");

          // Init growth rates
          largo_init     (generic_workspace, grDelta, grDA, grDArms);

          // Compute prestep values;
          largo_prestep_baseflow  (generic_workspace,   base,  dybase,
                                             dtbase, dxbase, srcbase);
          largo_prestep_seta_innery  (generic_workspace, y,
                                             mean,  rms,  mean_rqq,
                                            dmean, drms, dmean_rqq);
          largo_prestep_seta_innerxz (generic_workspace, field);

          // Compute mean sources;
          largo_continuity_setamean (generic_workspace, 0.0, 1.0, &srcmean[1-1]);
          largo_xmomentum_setamean  (generic_workspace, 0.0, 1.0, &srcmean[2-1]);
          largo_ymomentum_setamean  (generic_workspace, 0.0, 1.0, &srcmean[3-1]);
          largo_zmomentum_setamean  (generic_workspace, 0.0, 1.0, &srcmean[4-1]);
          largo_energy_setamean     (generic_workspace, 0.0, 1.0, &srcmean[5-1]);
          largo_species_setamean    (generic_workspace, 0.0, 1.0, &srcmean[5-0]);

          // Compute rms sources;
          largo_continuity_seta  (generic_workspace, 0.0, 1.0, &srcfull[1-1]);
          largo_xmomentum_seta   (generic_workspace, 0.0, 1.0, &srcfull[2-1]);
          largo_ymomentum_seta   (generic_workspace, 0.0, 1.0, &srcfull[3-1]);
          largo_zmomentum_seta   (generic_workspace, 0.0, 1.0, &srcfull[4-1]);
          largo_energy_seta      (generic_workspace, 0.0, 1.0, &srcfull[5-1]);
          largo_species_seta     (generic_workspace, 0.0, 1.0, &srcfull[5-0]);

          // Deallocate workspace
          largo_deallocate (&generic_workspace);

          // Check mean part
          // Tolerance for energy source adjusted manually
          // printf("%.16e, %.16e", srcmean[5-1], srcmean_good[5-1]);
          fct_chk_eq_dbl(srcmean[1-1]        , srcmean_good[1-1]        );
          fct_chk_eq_dbl(srcmean[2-1]/10.    , srcmean_good[2-1]/10.    );
          fct_chk_eq_dbl(srcmean[3-1]        , srcmean_good[3-1]        );
          fct_chk_eq_dbl(srcmean[4-1]        , srcmean_good[4-1]        );
          fct_chk_eq_dbl(srcmean[5-1]/100000., srcmean_good[5-1]/100000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcmean[5+is]      , srcmean_good[5+is]);
          }

          // Check full part
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcfull[1-1]        , srcfull_good[1-1]        );
          fct_chk_eq_dbl(srcfull[2-1]/10.    , srcfull_good[2-1]/10.    );
          fct_chk_eq_dbl(srcfull[3-1]        , srcfull_good[3-1]        );
          fct_chk_eq_dbl(srcfull[4-1]        , srcfull_good[4-1]        );
          fct_chk_eq_dbl(srcfull[5-1]/100000., srcfull_good[5-1]/100000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcfull[5+is]      , srcfull_good[5+is]);
          }


          // Recompute sources using wrapper methods
          // Allocate generic workspace;
          largo_allocate (&generic_workspace, model, neq, ns, 0, "dns");

          // Init growth rates
          largo_init     (generic_workspace, grDelta, grDA, grDArms);

          // Compute prestep values;
          largo_prestep_baseflow  (generic_workspace,   base,  dybase,
                                             dtbase, dxbase, srcbase);
          largo_prestep_seta (generic_workspace, y, field,
                                           mean,   rms,  mean_rqq,
                                          dmean,  drms, dmean_rqq);

          // Compute sources using wrapper method;
          largo_seta (generic_workspace, 0.0, 1.0, &srcfull[1-1]);

          // Check sources
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcfull[1-1]        , srcfull_good[1-1]        );
          fct_chk_eq_dbl(srcfull[2-1]/10.    , srcfull_good[2-1]/10.    );
          fct_chk_eq_dbl(srcfull[3-1]        , srcfull_good[3-1]        );
          fct_chk_eq_dbl(srcfull[4-1]        , srcfull_good[4-1]        );
          fct_chk_eq_dbl(srcfull[5-1]/100000., srcfull_good[5-1]/100000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcfull[5+is]      , srcfull_good[5+is]);
          }

          // Deallocate workspace
          largo_deallocate (&generic_workspace);

        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()