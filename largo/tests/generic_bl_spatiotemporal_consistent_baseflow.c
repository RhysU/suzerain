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
          const int neq = 7;
          const int ns  = 2;

          largo_workspace * work;
          double y        =   1.0/ 10.0;
          double grtDelta =   5.0/100.0;
          double uIw      = 460.0      ;
          double grxDelta = grtDelta / uIw;

          double field  [] = \
          {    11.0/  1000.0,  \
              485.0/   100.0,  \
                2.0/    10.0,  \
                3.0/    10.0,  \
            41500.0         ,  \
               22.0/ 10000.0,  \
               11.0/ 10000.0,  \
             4100.0
          };

          double mean  [] = \
          {     1.0/ 100.0,  \
               45.0/  10.0,  \
                1.0/1000.0,  \
                5.0/ 100.0,  \
            41200.0       ,  \
                3.0/1000.0,  \
                2.0/1000.0,  \
             4000.0
          };

          double dmean [] = \
          {      1.0/  5.0,  \
                45.0      ,  \
                 1.0/100.0,  \
                 5.0/ 10.0,  \
            412000.0      ,  \
                 2.0/100.0,  \
                 1.0/100.0,  \
             42000.0
          };

          double rms   [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double drms  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double mean_rqq [] =  \
          {           21.0/  2000.0,  \
                    8505.0/     4.0,  \
                      21.0/200000.0,  \
                      21.0/    80.0,  \
            178231200000.0         ,  \
                     189.0/200000.0,  \
                      21.0/ 50000.0,  \
                       0.0         
          };

          double dmean_rqq [] =  \
          {           11.0/    50.0, \
                    2025.0         , \
                       1.0/ 10000.0, \
                       1.0/     4.0, \
            169744000000.0         , \
                     -27.0/  5000.0, \
                     -19.0/  5000.0, \
                       0.0         
          };

          double grxDA   [] = \
          {   1.0/ 23000.0, \
              1.0/  4600.0, \
              1.0/ 46000.0, \
              3.0/ 46000.0, \
              1.0/ 11500.0, \
              1.0/460000.0, \
              1.0/230000.0, \
              1.0/ 46000.0  \
          };

          double grtDArms[] = \
          {    1.0/  1000.0, \
               5.0/  1000.0, \
               5.0/ 10000.0, \
               2.0/  1000.0, \
               1.0/  1000.0, \
               1.0/100000.0, \
               2.0/100000.0, \
               0.0           \
          };

          double grxDArms[] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          for (unsigned int i=0; i < neq+1; ++i)
          {
            grxDArms[i] = grtDArms[i] / uIw;
          }

          double base   [] = \
          {    5.0/  1000.0, \
              25.0/    10.0, \
               5.0/ 10000.0, \
               2.0/   100.0, \
           20000.0         , \
               1.0/  1000.0, \
               5.0/ 10000.0, \
            2000.0           \   
           };


          double dybase [] = \
          {    5.0/   100.0, \
              25.0         , \
               5.0/  1000.0, \
               2.0/    10.0, \
          200000.0         , \
               1.0/   100.0, \
               5.0/  1000.0, \
           20000.0           \    
          };


          double dtbase [] = \
          {    2.0/  1000.0, \
              15.0/    10.0, \
               2.0/ 10000.0, \
               1.0/   100.0, \
           10000.0         , \
               5.0/ 10000.0, \
               2.0/ 10000.0, \
            1000.0           \
           };

          double dxbase  [] = \
          {   1.0/ 230000.0, \
              3.0/    920.0, \
              1.0/2300000.0, \
              1.0/  46000.0, \
            500.0/     23.0, \
              1.0/ 920000.0, \
              1.0/2300000.0, \
             50.0/     23.0  \
          };

          double srcbase []=\
           {   1.0/  100000.0, \
               5.0/    1000.0, \
               1.0/ 1000000.0, \
               5.0/  100000.0, \
              50.0           , \
               2.0/ 1000000.0, \
               1.0/ 1000000.0
           };

          double wall_base [] = \
          { 1.0, uIw, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
                                                   
          double wall_ddy_base [] = \              
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
                                                   
          double wall_ddt_base [] = \              
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
                                                   
          double wall_ddx_base [] = \              
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
                                                   
          double wall_src_base [] = \              
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };

          double srcmean [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcfull  [] = \
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

          double srcmean_good [] = \
          { -   773.0/  200000.0,     \
            - 89543.0/   18400.0,     \
            - 20029.0/46000000.0,     \
            - 23899.0/  920000.0,     \
            -547450.0/      23.0,     \
            - 17857.0/11500000.0,     \
            - 10611.0/11500000.0
          };

          double srcfull_good [] = \
          { -     1701.0/      400000.0, \
            -   977843.0/      184000.0, \
            - 35720037.0/   460000000.0, \
            -   285217.0/     2300000.0, \
            - 56817639.0/        2300.0, \
            - 32015997.0/ 23000000000.0, \
            -  4018061.0/  5750000000.0  \
          };

          const double tolerance = 1.0E-12;

          const char model[] = "bl_spatiotemporal_consistent";

          largo_allocate (&work, model, neq, ns, 0, "dns");

          // Init growth rates
          largo_init     (work, grxDelta, grxDA, grxDArms);

          // Init wall baseflow
          largo_init_wall_baseflow(work   \ 
            ,wall_base ,wall_ddy_base ,wall_ddt_base \
                       ,wall_ddx_base ,wall_src_base \
            );

          // Compute prestep values;
          largo_prestep_baseflow  (work,   base,  dybase,
                                              dtbase, dxbase, srcbase);
          largo_prestep_seta_innery  (work, y,
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

          // Compute full sources;
          largo_continuity_seta     (work, 0.0, 1.0, &srcfull [1-1]);
          largo_xmomentum_seta      (work, 0.0, 1.0, &srcfull [2-1]);
          largo_ymomentum_seta      (work, 0.0, 1.0, &srcfull [3-1]);
          largo_zmomentum_seta      (work, 0.0, 1.0, &srcfull [4-1]);
          largo_energy_seta         (work, 0.0, 1.0, &srcfull [5-1]);
          largo_species_seta        (work, 0.0, 1.0, &srcfull [5-0]);

          // Deallocate workspace
          largo_deallocate (&work);

          // Check mean part
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcmean[1-1]        , srcmean_good[1-1]);
          fct_chk_eq_dbl(srcmean[2-1]/10.0   , srcmean_good[2-1]/10.0);
          fct_chk_eq_dbl(srcmean[3-1]        , srcmean_good[3-1]);
          fct_chk_eq_dbl(srcmean[4-1]        , srcmean_good[4-1]);
          fct_chk_eq_dbl(srcmean[5-1]/100000., srcmean_good[5-1]/100000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcmean[5+is]      , srcmean_good[5+is]);
          }

          // Check full source
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcfull[1-1]     , srcfull_good[1-1]);
          fct_chk_eq_dbl(srcfull[2-1]/100., srcfull_good[2-1]/100.);
          fct_chk_eq_dbl(srcfull[3-1]     , srcfull_good[3-1]);
          fct_chk_eq_dbl(srcfull[4-1]     , srcfull_good[4-1]);
          fct_chk_eq_dbl(srcfull[5-1]/100., srcfull_good[5-1]/100.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcfull[5+is]  , srcfull_good[5+is]);
          }


          // Recompute sources using wrapper methods
          // Allocate generic workspace;
          largo_allocate (&work, model, neq, ns, 0, "dns");

          // Init wall baseflow
          largo_init_wall_baseflow(work \ 
            ,wall_base ,wall_ddy_base ,wall_ddt_base \
                       ,wall_ddx_base ,wall_src_base \
            );

          // Init growth rates
          largo_init     (work, grxDelta, grxDA, grxDArms);

          // Compute prestep values;
          largo_prestep_baseflow  (work,   base,  dybase,
                                              dtbase, dxbase, srcbase);
          largo_prestep_seta (work, y, field,
                                           mean,   rms,  mean_rqq,
                                          dmean,  drms, dmean_rqq);

          // Compute sources using wrapper method;
          largo_seta (work, 0.0, 1.0, &srcfull[1-1]);

          // Check sources
          // Tolerance for energy source adjusted manually
          fct_chk_eq_dbl(srcfull[1-1]       , srcfull_good[1-1]);
          fct_chk_eq_dbl(srcfull[2-1]/100.0 , srcfull_good[2-1]/100.0);
          fct_chk_eq_dbl(srcfull[3-1]       , srcfull_good[3-1]);
          fct_chk_eq_dbl(srcfull[4-1]       , srcfull_good[4-1]);
          fct_chk_eq_dbl(srcfull[5-1]/10000., srcfull_good[5-1]/10000.);
          for (unsigned int is=0; is < ns; ++is)
          {
            fct_chk_eq_dbl(srcfull[5+is]    , srcfull_good[5+is]);
          }

          // Deallocate workspace
          largo_deallocate (&work);

        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
