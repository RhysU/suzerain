!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! largo 0.0.1: largo - slow growth terms for turbulence simulations
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2011-2014 The PECOS Development Team
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
!! $Id: bl_temporal_tconsistent_f.f90 38396 2013-03-29 16:41:25Z topalian $

#include "testframework_assert.h"

program bl_temporal_tconsistent_f

    use largo
    use testframework

    use, intrinsic :: iso_c_binding, only: c_associated,   &
                                           WP => c_double, &
                                           c_int,          &
                                           largo_workspace_ptr => c_ptr
    implicit none

    integer(c_int), parameter :: neq = 7
    integer(c_int), parameter :: ns  = 2

    type(largo_workspace_ptr) :: workspace
    real(WP), parameter               :: y       = 1.0_WP/ 10.0_WP
    real(WP), parameter               :: grDelta = 5.0_WP/100.0_WP

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
      &         1.0_WP/ 10.0_WP,  &
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
      &            11.0_WP/   100.0_WP, &
      &         22275.0_WP            , &
      &            11.0_WP/ 10000.0_WP, &
      &            11.0_WP/     4.0_WP, &
      & 1867184000000.0_WP            , &
      &             1.0_WP/ 50000.0_WP, &
      &             1.0_WP/200000.0_WP  &
      /)


    real(WP), dimension(neq), parameter :: &
      rms     = (/                &
      &      4.0_WP/ 10000.0_WP,  &
      &     25.0_WP/   100.0_WP,  &
      &     15.0_WP/   100.0_WP,  &
      &      1.0_WP/    10.0_WP,  &
      &    300.0_WP            ,  &
      &      0.0_WP/100000.0_WP,  &
      &      0.0_WP/100000.0_WP   &
      /)

    real(WP), dimension(neq), parameter :: &
      drms    = (/                 &
      &      42.0_WP/ 1000.0_WP,   &
      &      24.0_WP/   10.0_WP,   &
      &     153.0_WP/  100.0_WP,   &
      &      12.0_WP/   10.0_WP,   &
      &    3200.0_WP           ,   &
      &       0.0_WP/10000.0_WP,   &
      &       0.0_WP/10000.0_WP    &
      /)

    real(WP), dimension(neq)            :: srcmean
    real(WP), dimension(neq)            :: srcfull

    real(WP), dimension(neq), parameter :: &
      srcmean_good = (/            &
      &       1.0_WP/ 2000.0_WP,   &
      &       9.0_WP/   40.0_WP,   &
      &       1.0_WP/20000.0_WP,   &
      &       1.0_WP/  400.0_WP,   &
      &    2060.0_WP           ,   &
      &       1.0_WP/10000.0_WP,   &
      &       1.0_WP/20000.0_WP    &
      /)

    real(WP), dimension(neq), parameter :: &
      srcfull_good = (/                &
      &       11.0_WP/     20000.0_WP, &
      &       19.0_WP/        80.0_WP, &
      &     3989.0_WP/    200000.0_WP, &
      &      109.0_WP/      4000.0_WP, &
      &     1884.0_WP                , &
      &       11.0_WP/    100000.0_WP, &
      &       11.0_WP/    200000.0_WP  &
      /)

    real(WP), parameter :: tolerance = 1.0E-14

    integer(c_int) :: is

    call testframework_setup(__FILE__)

    ! Initialize srcmean and srcrms
    srcmean = 0.0_WP
    srcfull = 0.0_WP

    ! Allocate workspace
    call largo_BL_temporal_tconsistent_allocate (workspace, neq, ns)

    ! Compute prestep values
    call largo_BL_temporal_tconsistent_preStep_sEta ( workspace, y*grDelta, field, mean, rms, meanrqq, dmean, drms, dmeanrqq)

    ! Compute mean sources
    call largo_BL_temporal_tconsistent_continuity_sEtaMean (workspace, 0.0_WP, 1.0_WP, srcmean(1))
    call largo_BL_temporal_tconsistent_xMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(2))
    call largo_BL_temporal_tconsistent_yMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(3))
    call largo_BL_temporal_tconsistent_zMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(4))
    call largo_BL_temporal_tconsistent_energy_sEtaMean     (workspace, 0.0_WP, 1.0_WP, srcmean(5))
    call largo_BL_temporal_tconsistent_species_sEtaMean    (workspace, 0.0_WP, 1.0_WP, srcmean(6))


    ! Compute full sources
    call largo_BL_temporal_tconsistent_continuity_sEta_  (workspace, 0.0_WP, 1.0_WP, srcfull (1))
    call largo_BL_temporal_tconsistent_xMomentum_sEta_   (workspace, 0.0_WP, 1.0_WP, srcfull (2))
    call largo_BL_temporal_tconsistent_yMomentum_sEta_   (workspace, 0.0_WP, 1.0_WP, srcfull (3))
    call largo_BL_temporal_tconsistent_zMomentum_sEta_   (workspace, 0.0_WP, 1.0_WP, srcfull (4))
    call largo_BL_temporal_tconsistent_energy_sEta_      (workspace, 0.0_WP, 1.0_WP, srcfull (5))
    call largo_BL_temporal_tconsistent_species_sEta_     (workspace, 0.0_WP, 1.0_WP, srcfull (6))


    ! Check mean part
    ASSERT(abs((srcmean(1)/srcmean_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(2)/srcmean_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(3)/srcmean_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(4)/srcmean_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(5)/srcmean_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcmean(5+is)/srcmean_good(5+is))-1.0_WP) < tolerance )
    end do

    ! Check full part
    ASSERT(abs((srcfull(1) /srcfull_good(1))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(2) /srcfull_good(2))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(3) /srcfull_good(3))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(4) /srcfull_good(4))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(5) /srcfull_good(5))-1.0_WP)  < tolerance * 100.0_WP)
    do is=1, ns
      ASSERT(abs((srcfull(5+is)/srcfull_good(5+is))-1.0_WP) < tolerance * 10.0_WP )
    end do


    ! Deallocate workspace
    call largo_BL_temporal_tconsistent_deallocate (workspace)


    ! Recompute using wrapper routines
    srcfull = 0.0_WP

    ! Allocate workspace
    call largo_BL_temporal_tconsistent_allocate (workspace, neq, ns)

    ! Compute prestep values
    call largo_BL_temporal_tconsistent_preStep_sEta (workspace, y*grDelta, field, mean, rms, meanrqq, dmean, drms, dmeanrqq)

    ! Compute sources
    call largo_BL_temporal_tconsistent_sEta (workspace, 0.0_WP, 1.0_WP, srcfull(1))

    ! Check full part
    ASSERT(abs((srcfull(1) /srcfull_good(1))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(2) /srcfull_good(2))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(3) /srcfull_good(3))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(4) /srcfull_good(4))-1.0_WP)  < tolerance * 100.0_WP)
    ASSERT(abs((srcfull(5) /srcfull_good(5))-1.0_WP)  < tolerance * 100.0_WP)
    do is=1, ns
      ASSERT(abs((srcfull(5+is)/srcfull_good(5+is))-1.0_WP) < tolerance )
    end do

    ! Deallocate workspace
    call largo_BL_temporal_tconsistent_deallocate (workspace)

    call testframework_teardown()

end program bl_temporal_tconsistent_f
