!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! largo 0.0.1: largo - slow growth terms for turbulence simulations
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2011, 2012, 2013 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!! $Id: generic_bl_temporal_f.f90 40639 2013-08-07 19:07:34Z topalian $

#include "testframework_assert.h"

program generic_bl_temporal_baseflow_f

    use largo
    use testframework

    use, intrinsic :: iso_c_binding, only: c_associated,   &
                                           WP => c_double, &
                                           c_int
    implicit none

    integer(c_int), parameter :: neq = 7
    integer(c_int), parameter :: ns  = 2

    type(largo_ptr)                     :: generic_workspace
    real(WP), parameter                 :: y       = 1.0_WP/ 10.0_WP
    real(WP), parameter                 :: grDelta = 5.0_WP/100.0_WP

    character(len=255), parameter       :: turbmodel = "k_epsilon"
    integer(c_int), parameter           :: ntvar  = 2

    real(WP), dimension(neq), parameter :: &
      field   = (/                  &
      &        11.0_WP/ 1000.0_WP,  &
      &       485.0_WP/  100.0_WP,  &
      &         2.0_WP/   10.0_WP,  &
      &         3.0_WP/   10.0_WP,  &
      &     41500.0_WP           ,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/10000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      mean    = (/                &
      &        1.0_WP/ 100.0_WP,  &
      &       45.0_WP/  10.0_WP,  &
      &        1.0_WP/1000.0_WP,  &
      &        5.0_WP/ 100.0_WP,  &
      &    41200.0_WP          ,  &
      &        2.0_WP/1000.0_WP,  &
      &        1.0_WP/1000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      dmean   = (/                &
      &         1.0_WP/ 10.0_WP,  &
      &        45.0_WP         ,  &
      &         1.0_WP/100.0_WP,  &
      &         5.0_WP/ 10.0_WP,  &
      &    412000.0_WP         ,  &
      &         2.0_WP/100.0_WP,  &
      &         1.0_WP/100.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: & 
      rms     = (/                &
      &      4.0_WP/ 10000.0_WP,  &
      &     25.0_WP/   100.0_WP,  &
      &     15.0_WP/   100.0_WP,  &
      &      1.0_WP/    10.0_WP,  &
      &    300.0_WP            ,  &
      &      8.0_WP/100000.0_WP,  &
      &      4.0_WP/100000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      drms    = (/                 &
      &      42.0_WP/ 1000.0_WP,   & 
      &      24.0_WP/   10.0_WP,   &
      &     153.0_WP/  100.0_WP,   &
      &      12.0_WP/   10.0_WP,   &
      &    3200.0_WP           ,   &
      &      82.0_WP/10000.0_WP,   & 
      &      41.0_WP/10000.0_WP    & 
      /)

    real(WP), dimension(neq), parameter :: &
    mean_rqq  = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP /)

    real(WP), dimension(neq)            :: &
    dmean_rqq = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP /)

    real(WP), dimension(neq), parameter :: &
      grDA    = (/                 &
      &       2.0_WP/  100.0_WP,   & 
      &       1.0_WP/   10.0_WP,   &
      &       1.0_WP/  100.0_WP,   &
      &       3.0_WP/  100.0_WP,   &
      &       4.0_WP/  100.0_WP,   &
      &       1.0_WP/ 1000.0_WP,   & 
      &       2.0_WP/ 1000.0_WP    & 
      /)

    real(WP), dimension(neq), parameter :: &
      grDArms = (/                 &
      &      1.0_WP/  1000.0_WP,   & 
      &      5.0_WP/  1000.0_WP,   &
      &      5.0_WP/ 10000.0_WP,   &
      &      2.0_WP/  1000.0_WP,   &
      &      1.0_WP/  1000.0_WP,   &
      &      1.0_WP/100000.0_WP,   & 
      &      2.0_WP/100000.0_WP    & 
      /)

    real(WP), dimension(neq), parameter :: &
      base    = (/                  &
      &         5.0_WP/ 1000.0_WP,  &
      &        25.0_WP/   10.0_WP,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/  100.0_WP,  &
      &     20000.0_WP           ,  &
      &         1.0_WP/ 1000.0_WP,  &
      &         5.0_WP/10000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      dybase  = (/                &
      &         5.0_WP/ 100.0_WP,  &
      &        25.0_WP          ,  &
      &         5.0_WP/1000.0_WP,  &
      &         2.0_WP/  10.0_WP,  &
      &    200000.0_WP          ,  &
      &         1.0_WP/ 100.0_WP,  &
      &         5.0_WP/1000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      dtbase  = (/                &
      &         2.0_WP/ 1000.0_WP,  &
      &        15.0_WP/   10.0_WP,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/  100.0_WP,  &
      &     10000.0_WP           ,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/10000.0_WP   &
      /)

    real(WP), dimension(neq),parameter  :: &
      dxbase = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP /)

    real(WP), dimension(ntvar), parameter :: &
      meanTurb = (/               &
      &        5.0_WP/ 100.0_WP,  &
      &        3.0_WP/ 100.0_WP   &
      /)

    real(WP), dimension(ntvar), parameter :: &
      dymeanTurb = (/             &
      &        5.0_WP/  10.0_WP,  &
      &        3.0_WP/  10.0_WP   &
      /)

    real(WP), dimension(ntvar), parameter :: &
      grDAturb = (/               &
      &        4.0_WP/ 100.0_WP,  &
      &        2.0_WP/ 100.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: & 
      srcbase = (/                &
      &         1.0_WP/ 1000.0_WP,  &
      &         5.0_WP/   10.0_WP,  &
      &         1.0_WP/10000.0_WP,  &
      &         5.0_WP/ 1000.0_WP,  &
      &      5000.0_WP           ,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/10000.0_WP   &
      /)

    real(WP), dimension(neq)            :: srcmean
    real(WP), dimension(neq)            :: srcrms
    real(WP), dimension(neq)            :: srcall
    real(WP), dimension(ntvar)          :: srcturb

    real(WP), dimension(neq), parameter :: &
      srcmean_good = (/              &
      & -    17.0_WP/  20000.0_WP,   &
      & -    11.0_WP/     10.0_WP,   &
      & -     1.0_WP/  12500.0_WP,   &
      & -    11.0_WP/   2500.0_WP,   &
      & -  4788.0_WP             ,   &
      & -   251.0_WP/1000000.0_WP,   &
      & -    19.0_WP/ 250000.0_WP    &
      /)

    real(WP), dimension(neq), parameter :: &
      srcrms_good  = (/                 &
      &      131.0_WP/   250000.0_WP,   &
      &      301.0_WP/    20000.0_WP,   &
      &    20099.0_WP/  2000000.0_WP,   &
      &       29.0_WP/     2000.0_WP,   &
      &      157.0_WP/       10.0_WP,   &
      &  -461241.0_WP/500000000.0_WP,   &
      &   -28827.0_WP/ 62500000.0_WP    &
      /)

    real(WP), dimension(ntvar), parameter :: &
      srcturb_good  = (/                &
      &        1.0_WP/     2000.0_WP,   &
      &        9.0_WP/    10000.0_WP    &
      /)

    real(WP), dimension(neq)            :: &
      srcall_good = srcmean_good + srcrms_good

    real(WP), parameter :: tolerance = 1.0E-14

    integer(c_int) :: is, it

    call testframework_setup(__FILE__)

    ! Initialize srcmean and srcrms
    srcmean = 0.0_WP
    srcrms  = 0.0_WP
    srcall  = 0.0_WP
    srcturb = 0.0_WP

    ! Allocate workspace through generic interface
    call largo_allocate (generic_workspace, "bl_temporal", neq, ns, 0, "dns")

    ! Init growth rates
    call largo_init     (generic_workspace, grDelta, grDA, grDArms)

    ! Compute prestep values
    call largo_preStep_baseflow  (generic_workspace,   base,  dybase,  &
                                             dtbase, dxbase, srcbase)
    call largo_preStep_sEta_innery  (generic_workspace, y,         & 
                                              mean, rms, mean_rqq, &
                                           dmean, drms, dmean_rqq)
    call largo_preStep_sEta_innerxz (generic_workspace, field)


    ! Compute mean sources
    call largo_continuity_sEtaMean (generic_workspace, 0.0_WP, 1.0_WP, srcmean(1))
    call largo_xMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(2))
    call largo_yMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(3))
    call largo_zMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(4))
    call largo_energy_sEtaMean     (generic_workspace, 0.0_WP, 1.0_WP, srcmean(5))
    call largo_species_sEtaMean    (generic_workspace, 0.0_WP, 1.0_WP, srcmean(6))

!    do is=1, ns
!      call largo_BL_temporal_species_sEtaMean (generic_workspace, 0.0_WP, 1.0_WP, srcrms(5+is), is)
!    end do

    ! Compute rms sources
    call largo_continuity_sEtaRms  (generic_workspace, 0.0_WP, 1.0_WP, srcrms (1))
    call largo_xMomentum_sEtaRms   (generic_workspace, 0.0_WP, 1.0_WP, srcrms (2))
    call largo_yMomentum_sEtaRms   (generic_workspace, 0.0_WP, 1.0_WP, srcrms (3))
    call largo_zMomentum_sEtaRms   (generic_workspace, 0.0_WP, 1.0_WP, srcrms (4))
    call largo_energy_sEtaRms      (generic_workspace, 0.0_WP, 1.0_WP, srcrms (5))
    call largo_species_sEtaRms     (generic_workspace, 0.0_WP, 1.0_WP, srcrms (6))

!    do is=1, ns
!      call largo_BL_temporal_species_sEtaRms (generic_workspace, 0.0_WP, 1.0_WP, srcrms(5+is), is)
!    end do

    ! Check mean part
    ASSERT(abs((srcmean(1)/srcmean_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(2)/srcmean_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(3)/srcmean_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(4)/srcmean_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(5)/srcmean_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcmean(5+is)/srcmean_good(5+is))-1.0_WP) < tolerance )
    end do

    ! Check rms part
    ASSERT(abs((srcrms(1) /srcrms_good(1))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(2) /srcrms_good(2))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(3) /srcrms_good(3))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(4) /srcrms_good(4))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(5) /srcrms_good(5))-1.0_WP)  < tolerance )
    do is=1, ns
      ASSERT(abs((srcrms(5+is)/srcrms_good(5+is))-1.0_WP) < tolerance )
    end do

    ! Deallocate workspace
    call largo_deallocate (generic_workspace)

    ! Recompute using wrapper routines
    ! Allocate workspace through generic interface
    call largo_allocate (generic_workspace, "bl_temporal", neq, ns, &
      &    ntvar, turbmodel)

    ! Init growth rates
    call largo_init     (generic_workspace, grDelta, grDA, grDArms)
    call largo_init_rans(generic_workspace, grDAturb)

    ! Compute prestep values
    call largo_preStep_baseflow  (generic_workspace,   base,  dybase,  &
                                             dtbase, dxbase, srcbase)
    call largo_preStep_sEta (generic_workspace, y, field,         & 
                                           mean,  rms,  mean_rqq, &
                                          dmean, drms, dmean_rqq)
    call largo_preStep_sEtaMean_rans (generic_workspace, y,         &
                                                  meanTurb, dymeanTurb)

    ! Compute sources
    call largo_sEta      (generic_workspace, 0.0_WP, 1.0_WP, srcall(1))
    call largo_sEta_rans (generic_workspace, 0.0_WP, 1.0_WP, srcturb(1))

    ! Check all
    ASSERT(abs((srcall(1)/srcall_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(2)/srcall_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(3)/srcall_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(4)/srcall_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(5)/srcall_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcall(5+is)/srcall_good(5+is))-1.0_WP) < tolerance )
    end do
    do it=1, ntvar
      ASSERT(abs((srcturb(it)/srcturb_good(it))-1.0_WP) < tolerance )
    end do

    ! Deallocate workspace
    call largo_deallocate (generic_workspace)


    call testframework_teardown()

end program generic_bl_temporal_baseflow_f
