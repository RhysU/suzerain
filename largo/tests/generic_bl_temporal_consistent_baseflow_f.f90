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
!! $Id: bl_temporal_consistent_f.f90 38396 2013-03-29 16:41:25Z topalian $

#include "testframework_assert.h"

program generic_bl_temporal_consistent_f

    use largo
    use testframework

    use, intrinsic :: iso_c_binding, only: c_associated,   &
                                           WP => c_double, &
                                           c_int,          &
                                           largo_ptr => c_ptr, &
                                           largo_workspace_ptr => c_ptr
    implicit none

    integer(c_int), parameter :: neq = 7
    integer(c_int), parameter :: ns  = 2

    type(largo_ptr)           :: generic_workspace
    real(WP), parameter       :: y       = 1.0_WP/ 10.0_WP
    real(WP), parameter       :: grDelta = 5.0_WP/100.0_WP

    real(WP), dimension(neq), parameter :: &
      field   = (/                  &
      &        11.0_WP/ 1000.0_WP,  &
      &       485.0_WP/  100.0_WP,  &
      &         2.0_WP/   10.0_WP,  &
      &         3.0_WP/   10.0_WP,  &
      &     41500.0_WP           ,  &
      &        22.0_WP/10000.0_WP,  &
      &        11.0_WP/10000.0_WP   &
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
      &         1.0_WP/  5.0_WP,  &
      &        45.0_WP         ,  &
      &         1.0_WP/100.0_WP,  &
      &         5.0_WP/ 10.0_WP,  &
      &    412000.0_WP         ,  &
      &         2.0_WP/100.0_WP,  &
      &         1.0_WP/100.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      meanrqq = (/                &
      &            21.0_WP/  2000.0_WP, &
      &          8505.0_WP/     4.0_WP, &
      &            21.0_WP/200000.0_WP, &
      &            21.0_WP/    80.0_WP, &
      &  178231200000.0_WP            , &
      &            21.0_WP/ 50000.0_WP, &
      &            21.0_WP/200000.0_WP  &
      /)

    real(WP), dimension(neq), parameter :: &
      dmeanrqq = (/               &
      &            11.0_WP/    50.0_WP, &
      &          2025.0_WP            , &
      &             1.0_WP/ 10000.0_WP, &
      &             1.0_WP/     4.0_WP, &
      &  169744000000.0_WP            , &
      &             1.0_WP/  2500.0_WP, &
      &             1.0_WP/ 10000.0_WP  &
      /)


    real(WP), dimension(neq), parameter :: &
    rms  = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /)

    real(WP), dimension(neq)            :: &
    drms = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /)

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
      grDArms = (/                  &
      &       1.0_WP/  1000.0_WP,   & 
      &       5.0_WP/  1000.0_WP,   &
      &       5.0_WP/ 10000.0_WP,   &
      &       2.0_WP/  1000.0_WP,   &
      &       1.0_WP/  1000.0_WP,   &
      &       1.0_WP/100000.0_WP,   & 
      &       2.0_WP/100000.0_WP    & 
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

    real(WP), dimension(neq)            :: &
      dxbase  = 0d0


    real(WP), dimension(neq)            :: srcmean
    real(WP), dimension(neq)            :: srcfull

    real(WP), dimension(neq), parameter :: &
      srcmean_good = (/             &
      &     -27.0_WP/ 20000.0_WP,   &
      &    -713.0_WP/   400.0_WP,   &
      &     -37.0_WP/200000.0_WP,   &
      &    -271.0_WP/ 20000.0_WP,   &
      &  -11670.0_WP            ,   &
      &     -57.0_WP/100000.0_WP,   &
      &     -37.0_WP/200000.0_WP    &
      /)

    real(WP), dimension(neq), parameter :: &
      srcfull_good = (/             &
      &    -297.0_WP/  200000.0_WP,  &
      &   -7787.0_WP/    4000.0_WP,  &
      & -543089.0_WP/20000000.0_WP,  &
      &   -4847.0_WP/  100000.0_WP,  &
      & -307937.0_WP/      25.0_WP,  &
      &    -627.0_WP/ 1000000.0_WP,  &
      &    -407.0_WP/ 2000000.0_WP   &
      /)

    real(WP), parameter :: tolerance = 1.0E-14

    integer(c_int) :: is

    call testframework_setup(__FILE__)

    ! Initialize srcmean and srcrms
    srcmean = 0.0_WP
    srcfull = 0.0_WP

    ! Allocate workspace
    call largo_allocate (generic_workspace, "bl_temporal_consistent" , &
      &  neq, ns, 0, "dns")

    ! Init growth rates
    call largo_init  (generic_workspace, grDelta, grDA, grDArms)

    ! Compute prestep values
    call largo_preStep_baseflow  (generic_workspace,   base,  dybase,  &
                                             dtbase, dxbase, srcbase)
    call largo_preStep_sEta_innery  (generic_workspace, y, & 
                                             mean,  rms,  meanrqq, &
                                            dmean, drms, dmeanrqq)
    call largo_preStep_sEta_innerxz (generic_workspace, field)

    ! Compute mean sources
    call largo_continuity_sEtaMean (generic_workspace, 0.0_WP, 1.0_WP, srcmean(1))
    call largo_xMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(2))
    call largo_yMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(3))
    call largo_zMomentum_sEtaMean  (generic_workspace, 0.0_WP, 1.0_WP, srcmean(4))
    call largo_energy_sEtaMean     (generic_workspace, 0.0_WP, 1.0_WP, srcmean(5))
    call largo_species_sEtaMean    (generic_workspace, 0.0_WP, 1.0_WP, srcmean(6))


    ! Compute full sources
    call largo_continuity_sEta  (generic_workspace, 0.0_WP, 1.0_WP, srcfull (1))
    call largo_xMomentum_sEta   (generic_workspace, 0.0_WP, 1.0_WP, srcfull (2))
    call largo_yMomentum_sEta   (generic_workspace, 0.0_WP, 1.0_WP, srcfull (3))
    call largo_zMomentum_sEta   (generic_workspace, 0.0_WP, 1.0_WP, srcfull (4))
    call largo_energy_sEta      (generic_workspace, 0.0_WP, 1.0_WP, srcfull (5))
    call largo_species_sEta     (generic_workspace, 0.0_WP, 1.0_WP, srcfull (6))


    ! Check mean part
    ASSERT(abs((srcmean(1)/srcmean_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(2)/srcmean_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(3)/srcmean_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(4)/srcmean_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(5)/srcmean_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcmean(is)/srcmean_good(is))-1.0_WP) < tolerance )
    end do

    ! Check full part
    ASSERT(abs((srcfull(1) /srcfull_good(1))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(2) /srcfull_good(2))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(3) /srcfull_good(3))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(4) /srcfull_good(4))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(5) /srcfull_good(5))-1.0_WP)  < tolerance * 100_WP)
    do is=1, ns
      ASSERT(abs((srcfull(is)/srcfull_good(is))-1.0_WP) < tolerance )
    end do

    ! Deallocate workspace
    call largo_deallocate (generic_workspace)


    ! Recompute using wrapper methods
    srcfull = 0.0_WP

    ! Allocate workspace
    call largo_allocate (generic_workspace, "bl_temporal_consistent" , &
      &  neq, ns, 0, "dns")

    ! Init growth rates
    call largo_init  (generic_workspace, grDelta, grDA, grDArms)

    ! Compute prestep values
    call largo_preStep_baseflow  (generic_workspace,   base,  dybase,  &
                                             dtbase, dxbase, srcbase)
    call largo_preStep_sEta (generic_workspace, y, field, mean, rms, meanrqq, dmean, drms, dmeanrqq)

    ! Compute sources
    call largo_sEta (generic_workspace, 0.0_WP, 1.0_WP, srcfull(1))

    ! Check full part
    ASSERT(abs((srcfull(1) /srcfull_good(1))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(2) /srcfull_good(2))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(3) /srcfull_good(3))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(4) /srcfull_good(4))-1.0_WP)  < tolerance * 100_WP)
    ASSERT(abs((srcfull(5) /srcfull_good(5))-1.0_WP)  < tolerance * 100_WP)
    do is=1, ns
      ASSERT(abs((srcfull(is)/srcfull_good(is))-1.0_WP) < tolerance )
    end do

    ! Deallocate workspace
    call largo_deallocate (generic_workspace)

    call testframework_teardown()

end program generic_bl_temporal_consistent_f
