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

#ifndef LARGO_LARGO_BL_temporal_tconsistent_H
#define LARGO_LARGO_BL_temporal_tconsistent_H

#ifdef __cplusplus
extern "C" {
#endif

/* Add declarations here for the public C API */
void  largo_bl_temporal_tconsistent_allocate (void * workspace, 
                                  int    neq, 
                                  int    ns);

void  largo_bl_temporal_tconsistent_prestep_setamean ( double y, 
                                           double * mean, 
                                           double * ddy_mean, 
                                           void * workspace);

void  largo_bl_temporal_tconsistent_prestep_seta_ ( double y, 
                                          double * rms,
                                          double * meanrqq, 
                                          double * ddy_rms,
                                          double * dmeanrqq, 
                                          void * workspace);

void  largo_BL_temporal_tconsistent_preStep_seta_innerxz( double * qflow, 
                                              double * mean,     
                                              void * workspace); 
                                              
void  largo_BL_temporal_tconsistent_preStep_seta_innery ( double y, 
                                              double * mean, 
                                              double * rms,
                                              double * meanrqq, 
                                              double * ddy_mean, 
                                              double * ddy_rms,
                                              double * dmeanrqq, 
                                              void * workspace);

void  largo_bl_temporal_tconsistent_prestep_seta ( double y, 
                                       double * qflow, 
                                       double * mean, 
                                       double * rms, 
                                       double * meanrqq, 
                                       double * ddy_mean, 
                                       double * ddy_rms,
                                       double * dmeanrqq, 
                                       void * workspace);

void  largo_bl_temporal_tconsistent_continuity_setamean (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_xmomentum_setamean  (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_ymomentum_setamean  (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_zmomentum_setamean  (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_energy_setamean     (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_species_setamean    (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_ispecies_setamean   (void * workspace, double A, double B, double * src, int is);
void  largo_bl_temporal_tconsistent_continuity_setarms  (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_xmomentum_setarms   (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_ymomentum_setarms   (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_zmomentum_setarms   (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_energy_setarms      (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_species_setarms     (void * workspace, double A, double B, double * src);
void  largo_bl_temporal_tconsistent_ispecies_setarms    (void * workspace, double A, double B, double * src, int is);
void  largo_bl_temporal_tconsistent_deallocate (void * workspace);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LARGO_LARGO_BL_temporal_tconsistent_H */
