//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// largo 0.0.1: largo - slow growth terms for turbulence simulations
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

#ifndef LARGO_LARGO_H
#define LARGO_LARGO_H

#include <largo/version.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declare opaque generic Largo model workspace */
struct largo_workspace;
typedef struct largo_workspace largo_workspace;

void largo_allocate ( largo_workspace** work,
                      const char * model,
                      int    neq,
                      int    ns,
                      int    ntvar,
                      const char * ransmodel );

void largo_init ( largo_workspace* work,
                  double  grDelta,
                  const double* grDA,
                  const double* grDArms );

void largo_init_rans ( largo_workspace* work,
                       double  grDelta,
                       double* grDA,
                       double* grDArms );

void largo_init_wall_baseflow ( largo_workspace* work,
                                double*     wall,
                                double* ddy_wall,
                                double* ddt_wall,
                                double* ddx_wall,
                                double* src_wall );

void largo_prestep_baseflow ( largo_workspace* work,
                              double*     base,
                              double* ddy_base,
                              double* ddt_base,
                              double* ddx_base,
                              double* src_base );

void largo_prestep_setamean ( largo_workspace* work,
                              double  y,
                              double* mean,
                              double* ddy_mean );

void largo_prestep_setamean_rans ( largo_workspace* work,
                                   double  y,
                                   double* mean,
                                   double* ddy_mean );

void largo_prestep_setarms ( largo_workspace* work,
                             double  y,
                             double* rms,
                             double* ddy_rms );

void largo_prestep_seta_innerxz ( largo_workspace* work,
                                  double* qflow );

void largo_prestep_seta_innery ( largo_workspace* work,
                                 double  y,
                                 double* mean,
                                 double* rms,
                                 double* mean_rqq,
                                 double* ddy_mean,
                                 double* ddy_rms,
                                 double* ddy_mean_rqq );

void largo_prestep_seta ( largo_workspace* work,
                          double  y,
                          double* qflow,
                          double* mean,
                          double* rms,
                          double* mean_rqq,
                          double* ddy_mean,
                          double* ddy_rms,
                          double* ddy_mean_rqq );


void largo_continuity_setamean ( largo_workspace* work, double A, double B, double* src );
void largo_xmomentum_setamean  ( largo_workspace* work, double A, double B, double* src );
void largo_ymomentum_setamean  ( largo_workspace* work, double A, double B, double* src );
void largo_zmomentum_setamean  ( largo_workspace* work, double A, double B, double* src );
void largo_energy_setamean     ( largo_workspace* work, double A, double B, double* src );
void largo_species_setamean    ( largo_workspace* work, double A, double B, double* src );
void largo_ispecies_setamean   ( largo_workspace* work, double A, double B, double* src, int is );
void largo_continuity_setarms  ( largo_workspace* work, double A, double B, double* src );
void largo_xmomentum_setarms   ( largo_workspace* work, double A, double B, double* src );
void largo_ymomentum_setarms   ( largo_workspace* work, double A, double B, double* src );
void largo_zmomentum_setarms   ( largo_workspace* work, double A, double B, double* src );
void largo_energy_setarms      ( largo_workspace* work, double A, double B, double* src );
void largo_species_setarms     ( largo_workspace* work, double A, double B, double* src );
void largo_ispecies_setarms    ( largo_workspace* work, double A, double B, double* src, int is );
void largo_continuity_seta     ( largo_workspace* work, double A, double B, double* src );
void largo_xmomentum_seta      ( largo_workspace* work, double A, double B, double* src );
void largo_ymomentum_seta      ( largo_workspace* work, double A, double B, double* src );
void largo_zmomentum_seta      ( largo_workspace* work, double A, double B, double* src );
void largo_energy_seta         ( largo_workspace* work, double A, double B, double* src );
void largo_species_seta        ( largo_workspace* work, double A, double B, double* src );

void largo_setamean            ( largo_workspace* work, double A, double B, double* src );
void largo_seta                ( largo_workspace* work, double A, double B, double* src );
void largo_seta_rans           ( largo_workspace* work, double A, double B, double* src );

void largo_deallocate          ( largo_workspace** work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LARGO_LARGO_H */
