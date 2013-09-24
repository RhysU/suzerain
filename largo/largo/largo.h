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

#ifndef LARGO_LARGO_H
#define LARGO_LARGO_H

#include <largo/version.h>

#ifdef __cplusplus
extern "C" {
#endif

void  largo_allocate (void * generic_workspace,
                      int    neq,
                      int    ns,
                      const char * model);

void  largo_init ( void * generic_workspace,
                   double grDelta,
                   double * grDA,
                   double * grDArms);

void  largo_init_wall_baseflow ( void * generic_workspace,
                                 double *     wall,
                                 double * ddy_wall,
                                 double * ddt_wall,
                                 double * ddx_wall,
                                 double * src_wall);

void  largo_prestep_baseflow ( void * generic_workspace,
                               double *     base,
                               double * ddy_base,
                               double * ddt_base,
                               double * ddx_base,
                               double * src_base);

void  largo_prestep_setamean ( void * generic_workspace,
                               double y,
                               double * mean,
                               double * ddy_mean);

void  largo_prestep_setarms ( void * generic_workspace,
                              double y,
                              double * rms,
                              double * ddy_rms);

void  largo_prestep_seta_innerxz( void * generic_workspace,
                                  double * qflow);

void  largo_prestep_seta_innery ( void * generic_workspace,
                                  double y,
                                  double * mean,
                                  double * rms,
                                  double * mean_rqq,
                                  double * ddy_mean,
                                  double * ddy_rms,
                                  double * ddy_mean_rqq);

void  largo_prestep_seta ( void * generic_workspace,
                           double y,
                           double * qflow,
                           double * mean,
                           double * rms,
                           double * mean_rqq,
                           double * ddy_mean,
                           double * ddy_rms,
                           double * ddy_mean_rqq);


void  largo_continuity_setamean (void * generic_workspace, double A, double B, double * src);
void  largo_xmomentum_setamean  (void * generic_workspace, double A, double B, double * src);
void  largo_ymomentum_setamean  (void * generic_workspace, double A, double B, double * src);
void  largo_zmomentum_setamean  (void * generic_workspace, double A, double B, double * src);
void  largo_energy_setamean     (void * generic_workspace, double A, double B, double * src);
void  largo_species_setamean    (void * generic_workspace, double A, double B, double * src);
void  largo_ispecies_setamean   (void * generic_workspace, double A, double B, double * src, int is);
void  largo_continuity_setarms  (void * generic_workspace, double A, double B, double * src);
void  largo_xmomentum_setarms   (void * generic_workspace, double A, double B, double * src);
void  largo_ymomentum_setarms   (void * generic_workspace, double A, double B, double * src);
void  largo_zmomentum_setarms   (void * generic_workspace, double A, double B, double * src);
void  largo_energy_setarms      (void * generic_workspace, double A, double B, double * src);
void  largo_species_setarms     (void * generic_workspace, double A, double B, double * src);
void  largo_ispecies_setarms    (void * generic_workspace, double A, double B, double * src, int is);
void  largo_continuity_seta     (void * generic_workspace, double A, double B, double * src);
void  largo_xmomentum_seta      (void * generic_workspace, double A, double B, double * src);
void  largo_ymomentum_seta      (void * generic_workspace, double A, double B, double * src);
void  largo_zmomentum_seta      (void * generic_workspace, double A, double B, double * src);
void  largo_energy_seta         (void * generic_workspace, double A, double B, double * src);
void  largo_species_seta        (void * generic_workspace, double A, double B, double * src);

void  largo_setamean            (void * generic_workspace, double A, double B, double * src);
void  largo_seta                (void * generic_workspace, double A, double B, double * src);

void  largo_deallocate          (void * generic_workspace);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LARGO_LARGO_H */
