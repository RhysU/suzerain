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

          largo_workspace * work;
          double y       = 1.0/ 10.0;
          double grDelta = 5.0/100.0;

          double field  [] = \
          {    11.0/ 1000.0,  \
              485.0/  100.0,  \
                2.0/   10.0,  \
                3.0/   10.0,  \
            41500.0        ,  \
                2.0/10000.0,  \
                1.0/10000.0
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
          {      1.0/ 10.0,  \
                45.0      ,  \
                 1.0/100.0,  \
                 5.0/ 10.0,  \
            412000.0      ,  \
                 2.0/100.0,  \
                 1.0/100.0
          };

          double rms   [] = \
          {   4.0/ 10000.0,  \
             25.0/   100.0,  \
             15.0/   100.0,  \
              1.0/    10.0,  \
            300.0         ,  \
              8.0/100000.0,  \
              4.0/100000.0
          };

          double drms  [] = \
          {    42.0/ 1000.0, \
               24.0/   10.0, \
              153.0/  100.0, \
               12.0/   10.0, \
             3200.0        , \
               82.0/10000.0, \
               41.0/10000.0
          };

          double mean_rqq [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double dmean_rqq  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcmean [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcrms  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcall  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcmean_good [] = \
          {    1.0/ 2000.0,         \
               9.0/   40.0,         \
               1.0/20000.0,         \
               1.0/  400.0,         \
            2060.0        ,         \
               1.0/10000.0,         \
               1.0/20000.0          \
          };

          double srcrms_good  [] = \
          {    21.0/  40000.0,      \
               21.0/   1250.0,      \
            10149.0/1000000.0,      \
                3.0/    200.0,      \
               16.0          ,      \
             -369.0/ 400000.0,      \
             -369.0/ 800000.0       \
          };

          double srcall_good  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          for (unsigned int i=0; i < neq; ++i)
          {
            srcall_good[i] = srcmean_good[i] + srcrms_good[i];
          }

          const double tolerance = 1.0E-12;

          const char model[] = "bl_temporal";
          largo_allocate (&work, model, neq, ns, 0, "dns");

          // Compute prestep values;
          largo_prestep_seta_innery  (work, y*grDelta,
                                             mean,  rms,  mean_rqq,
                                            dmean, drms, dmean_rqq);
          largo_prestep_seta_innerxz (work, field);

          // Compute mean sources;
          largo_continuity_setamean (work, 0.0, 1.0, &srcmean[1-1]);
          largo_xmomentum_setamean  (work, 0.0, 1.0, &srcmean[2-1]);
          largo_ymomentum_setamean  (work, 0.0, 1.0, &srcmean[3-1]);
          largo_zmomentum_setamean  (work, 0.0, 1.0, &srcmean[4-1]);
          largo_energy_setamean     (work, 0.0, 1.0, &srcmean[5-1]);
          largo_species_setamean    (work, 0.0, 1.0, &srcmean[5-0]);
//           for (int is=0; is < ns; ++is)
//           {
//             largo_bl_temporal_ispecies_setamean (work, 0.0, 1.0, &srcmean[5+is], is+1);
//           }

          // Compute rms sources;
          largo_continuity_setarms  (work, 0.0, 1.0, &srcrms [1-1]);
          largo_xmomentum_setarms   (work, 0.0, 1.0, &srcrms [2-1]);
          largo_ymomentum_setarms   (work, 0.0, 1.0, &srcrms [3-1]);
          largo_zmomentum_setarms   (work, 0.0, 1.0, &srcrms [4-1]);
          largo_energy_setarms      (work, 0.0, 1.0, &srcrms [5-1]);
          largo_species_setarms     (work, 0.0, 1.0, &srcrms [5-0]);
//           for (int is=0; is < ns; ++is)
//           {
//             largo_bl_temporal_ispecies_setarms (work, 0.0, 1.0, &srcrms[5+is], is+1);
//           }

          // Deallocate workspace
          largo_deallocate (&work);

          // Check mean part
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcmean[1-1]       , srcmean_good[1-1]);
          fct_chk_eq_dbl(srcmean[2-1]       , srcmean_good[2-1]);
          fct_chk_eq_dbl(srcmean[3-1]       , srcmean_good[3-1]);
          fct_chk_eq_dbl(srcmean[4-1]       , srcmean_good[4-1]);
          fct_chk_eq_dbl(srcmean[5-1]/10000., srcmean_good[5-1]/10000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcmean[5+is]      , srcmean_good[5+is]);
          }

          // Check rms part
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcrms[1-1]     , srcrms_good[1-1]);
          fct_chk_eq_dbl(srcrms[2-1]     , srcrms_good[2-1]);
          fct_chk_eq_dbl(srcrms[3-1]     , srcrms_good[3-1]);
          fct_chk_eq_dbl(srcrms[4-1]     , srcrms_good[4-1]);
          fct_chk_eq_dbl(srcrms[5-1]/100., srcrms_good[5-1]/100.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcrms[5+is]      , srcrms_good[5+is]);
          }


          // Recompute sources using wrapper methods
          // Allocate generic workspace;
          largo_allocate (&work, model, neq, ns, 0, "dns");

          // Compute prestep values;
          largo_prestep_seta (work, y*grDelta, field,
                                           mean,   rms,  mean_rqq,
                                          dmean,  drms, dmean_rqq);

          // Compute sources using wrapper method;
          largo_seta (work, 0.0, 1.0, &srcall[1-1]);

          // Check sources
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcall[1-1]       , srcall_good[1-1]);
          fct_chk_eq_dbl(srcall[2-1]       , srcall_good[2-1]);
          fct_chk_eq_dbl(srcall[3-1]       , srcall_good[3-1]);
          fct_chk_eq_dbl(srcall[4-1]       , srcall_good[4-1]);
          fct_chk_eq_dbl(srcall[5-1]/10000., srcall_good[5-1]/10000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcall[5+is]      , srcall_good[5+is]);
          }

          // Deallocate workspace
          largo_deallocate (&work);

        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
