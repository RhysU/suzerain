!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! largo 0.0.2: largo - slow growth terms for turbulence simulations
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
!! $Id: bl_spatiotemporal_f.f90 40608 2013-08-06 23:12:28Z topalian $

#include "testframework_assert.h"

program bl_spatiotemporal_baseflow_f

    use largo
    use testframework

    use, intrinsic :: iso_c_binding, only: c_associated,   &
                                           WP => c_double, &
                                           c_int
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

    implicit none

    integer(c_int), parameter :: neq = 7
    integer(c_int), parameter :: ns  = 2

    type(largo_workspace_ptr)         :: workspace
    real(WP), parameter               :: y        = 1.0_WP/ 10.0_WP
    real(WP), parameter               :: grtDelta = 5.0_WP/100.0_WP
    real(WP), parameter               :: uIw      = 460.0_WP
    real(WP), parameter               :: grxDelta = grtDelta / uIw

    real(WP), dimension(neq+1), parameter :: &
      field   = (/                  &
      &        11.0_WP/ 1000.0_WP,  &
      &       485.0_WP/  100.0_WP,  &
      &         2.0_WP/   10.0_WP,  &
      &         3.0_WP/   10.0_WP,  &
      &     41500.0_WP           ,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/10000.0_WP,  &
      &      4100.0_WP              &
      /)

    real(WP), dimension(neq+1), parameter :: &
      mean    = (/                &
      &        1.0_WP/ 100.0_WP,  &
      &       45.0_WP/  10.0_WP,  &
      &        1.0_WP/1000.0_WP,  &
      &        5.0_WP/ 100.0_WP,  &
      &    41200.0_WP          ,  &
      &        2.0_WP/1000.0_WP,  &
      &        1.0_WP/1000.0_WP,  &
      &     4000.0_WP             &
      /)

    real(WP), dimension(neq+1), parameter :: &
      dmean   = (/                &
      &         1.0_WP/ 10.0_WP,  &
      &        45.0_WP         ,  &
      &         1.0_WP/100.0_WP,  &
      &         5.0_WP/ 10.0_WP,  &
      &    412000.0_WP         ,  &
      &         2.0_WP/100.0_WP,  &
      &         1.0_WP/100.0_WP,  &
      &     42000.0_WP            &
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

    real(WP), dimension(neq+1), parameter :: &
    mean_rqq  = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /)

    real(WP), dimension(neq+1)            :: &
    dmean_rqq = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /)

    real(WP), dimension(neq+1), parameter :: &
      grxDA    = (/             &
      &       1.0_WP/ 23000.0_WP,   &
      &       1.0_WP/  4600.0_WP,   &
      &       1.0_WP/ 46000.0_WP,   &
      &       3.0_WP/ 46000.0_WP,   &
      &       1.0_WP/ 11500.0_WP,   &
      &       1.0_WP/460000.0_WP,   &
      &       1.0_WP/230000.0_WP,   &
      &       1.0_WP/ 46000.0_WP    &
      /)

    real(WP), dimension(neq), parameter :: &
      grtDArms = (/             &
      &       1.0_WP/  1000.0_WP,   &
      &       5.0_WP/  1000.0_WP,   &
      &       5.0_WP/ 10000.0_WP,   &
      &       2.0_WP/  1000.0_WP,   &
      &       1.0_WP/  1000.0_WP,   &
      &       1.0_WP/100000.0_WP,   &
      &       2.0_WP/100000.0_WP    &
      /)

    real(WP), dimension(neq), parameter :: &
      grxDArms = grtDArms / uIw

    real(WP), dimension(neq+1), parameter :: &
      base    = (/                 &
      &         5.0_WP/ 1000.0_WP,  &
      &        25.0_WP/   10.0_WP,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/  100.0_WP,  &
      &     20000.0_WP           ,  &
      &         1.0_WP/ 1000.0_WP,  &
      &         5.0_WP/10000.0_WP,  &
      &      2000.0_WP              &
      /)

#ifndef BASEFLOW_UNIFORM
    real(WP), dimension(neq+1), parameter :: &
      dybase  = (/                &
      &         5.0_WP/ 100.0_WP,  &
      &        25.0_WP          ,  &
      &         5.0_WP/1000.0_WP,  &
      &         2.0_WP/  10.0_WP,  &
      &    200000.0_WP          ,  &
      &         1.0_WP/ 100.0_WP,  &
      &         5.0_WP/1000.0_WP,  &
      &     20000.0_WP             &
      /)

    real(WP), dimension(neq+1), parameter :: &
      dtbase  = (/                &
      &         2.0_WP/ 1000.0_WP,  &
      &        15.0_WP/   10.0_WP,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/  100.0_WP,  &
      &     10000.0_WP           ,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/10000.0_WP,  &
      &      1000.0_WP              &
      /)

    real(WP), dimension(neq+1)            :: &
      dxbase  = (/                &
      &       1.0_WP/ 230000.0_WP,  &
      &       3.0_WP/    920.0_WP,  &
      &       1.0_WP/2300000.0_WP,  &
      &       1.0_WP/  46000.0_WP,  &
      &     500.0_WP/     23.0_WP,  &
      &       1.0_WP/ 920000.0_WP,  &
      &       1.0_WP/2300000.0_WP,  &
      &      50.0_WP/     23.0_WP   &
      /)
#else
    real(WP), dimension(neq), parameter ::  dybase = 0.0_WP &
                                         ,  dtbase = 0.0_WP &
                                         ,  dxbase = 0.0_WP
#endif

    real(WP), dimension(neq), parameter :: &
      srcbase = (/                &
      &         1.0_WP/ 100000.0_WP,  &
      &         5.0_WP/   1000.0_WP,  &
      &         1.0_WP/1000000.0_WP,  &
      &         5.0_WP/ 100000.0_WP,  &
      &        50.0_WP           ,  &
      &         2.0_WP/1000000.0_WP,  &
      &         1.0_WP/1000000.0_WP   &
      /)

    real(WP), dimension(neq)            :: srcmean
    real(WP), dimension(neq)            :: srcrms
    real(WP), dimension(neq)            :: srcall

    real(WP), dimension(neq), parameter :: &
      srcmean_good = (/              &
      & -  7977.0_WP/ 2300000.0_WP,   &
      & - 78923.0_WP/   18400.0_WP,   &
      & - 15729.0_WP/46000000.0_WP,   &
      & - 16089.0_WP/  920000.0_WP,   &
      & -413010.0_WP/      23.0_WP,   &
      & - 35553.0_WP/46000000.0_WP,   &
      & - 15549.0_WP/46000000.0_WP    &
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

    real(WP), dimension(neq)            :: &
      srcall_good = srcmean_good + srcrms_good

    real(WP), parameter :: tolerance = 1.0E-14

    integer(c_int) :: is

    call testframework_setup(__FILE__)

    ! Initialize srcmean and srcrms
    srcmean = 0.0_WP
    srcrms  = 0.0_WP
    srcall  = 0.0_WP

    ! Allocate workspace
    call largo_BL_spatiotemporal_allocate (workspace, neq, ns)

    ! Init growth rates
    call largo_BL_spatiotemporal_init  (workspace, grxDelta, grxDA, grxDArms)

    ! Init wall baseflow
    call largo_BL_spatiotemporal_init_wall_baseflow(workspace   &
      ,(/ 1.0_WP,    uIw, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      )

    ! Compute prestep values
    call largo_BL_spatiotemporal_preStep_baseflow  (workspace,   base,  dybase,  &
                                                 dtbase, dxbase, srcbase)
    call largo_BL_spatiotemporal_preStep_sEta_innery  (workspace, y,           &
                                                  mean,  rms,  mean_rqq, &
                                                 dmean, drms, dmean_rqq)
    call largo_BL_spatiotemporal_preStep_sEta_innerxz (workspace, field)

    ! Compute mean sources
    call largo_BL_spatiotemporal_continuity_sEtaMean (workspace, 0.0_WP, 1.0_WP, srcmean(1))
    call largo_BL_spatiotemporal_xMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(2))
    call largo_BL_spatiotemporal_yMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(3))
    call largo_BL_spatiotemporal_zMomentum_sEtaMean  (workspace, 0.0_WP, 1.0_WP, srcmean(4))
    call largo_BL_spatiotemporal_energy_sEtaMean     (workspace, 0.0_WP, 1.0_WP, srcmean(5))
    call largo_BL_spatiotemporal_species_sEtaMean    (workspace, 0.0_WP, 1.0_WP, srcmean(6))

!    do is=1, ns
!      call largo_BL_spatiotemporal_species_sEtaMean (workspace, 0.0_WP, 1.0_WP, srcrms(5+is), is)
!    end do

    ! Compute rms sources
    call largo_BL_spatiotemporal_continuity_sEtaRms  (workspace, 0.0_WP, 1.0_WP, srcrms (1))
    call largo_BL_spatiotemporal_xMomentum_sEtaRms   (workspace, 0.0_WP, 1.0_WP, srcrms (2))
    call largo_BL_spatiotemporal_yMomentum_sEtaRms   (workspace, 0.0_WP, 1.0_WP, srcrms (3))
    call largo_BL_spatiotemporal_zMomentum_sEtaRms   (workspace, 0.0_WP, 1.0_WP, srcrms (4))
    call largo_BL_spatiotemporal_energy_sEtaRms      (workspace, 0.0_WP, 1.0_WP, srcrms (5))
    call largo_BL_spatiotemporal_species_sEtaRms     (workspace, 0.0_WP, 1.0_WP, srcrms (6))

!    do is=1, ns
!      call largo_BL_spatiotemporal_species_sEtaRms (workspace, 0.0_WP, 1.0_WP, srcrms(5+is), is)
!    end do

    ! Check mean part
    if (any(isnan(srcmean))) write (error_unit, *) "srcmean: ", srcmean
    ASSERT(.not.any(isnan(srcmean)))
#ifndef BASEFLOW_UNIFORM
    ASSERT(abs((srcmean(1)/srcmean_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(2)/srcmean_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(3)/srcmean_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(4)/srcmean_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcmean(5)/srcmean_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcmean(5+is)/srcmean_good(5+is))-1.0_WP) < tolerance )
    end do
#endif

    ! Check rms part
    if (any(isnan(srcrms))) write (error_unit, *) "srcrms: ", srcrms
    ASSERT(.not.any(isnan(srcrms)))
#ifndef BASEFLOW_UNIFORM
    ASSERT(abs((srcrms(1) /srcrms_good(1))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(2) /srcrms_good(2))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(3) /srcrms_good(3))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(4) /srcrms_good(4))-1.0_WP)  < tolerance )
    ASSERT(abs((srcrms(5) /srcrms_good(5))-1.0_WP)  < tolerance )
    do is=1, ns
      ASSERT(abs((srcrms(5+is)/srcrms_good(5+is))-1.0_WP) < tolerance )
    end do
#endif

    ! Deallocate workspace
    call largo_BL_spatiotemporal_deallocate (workspace)


    ! Recompute using wrapper routines
    ! Allocate workspace (same pointer as before)
    call largo_BL_spatiotemporal_allocate (workspace, neq, ns)

    ! Init growth rates
    call largo_BL_spatiotemporal_init  (workspace, grxDelta, grxDA, grxDArms)

    ! Init wall baseflow
    call largo_BL_spatiotemporal_init_wall_baseflow(workspace   &
      ,(/ 1.0_WP,    uIw, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      ,(/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP  /) &
      )

    ! Compute prestep values
    call largo_BL_spatiotemporal_preStep_baseflow  (workspace,   base,  dybase,  &
                                                 dtbase, dxbase, srcbase)
    call largo_BL_spatiotemporal_preStep_sEta (workspace, y, field,         &
                                               mean,  rms,  mean_rqq, &
                                              dmean, drms, dmean_rqq)

    ! Compute sources
    call largo_BL_spatiotemporal_sEta (workspace, 0.0_WP, 1.0_WP, srcall(1))

    ! Check all
    if (any(isnan(srcall))) write (error_unit, *) "srcall: ", srcall
    ASSERT(.not.any(isnan(srcall)))
#ifndef BASEFLOW_UNIFORM
    ASSERT(abs((srcall(1)/srcall_good(1))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(2)/srcall_good(2))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(3)/srcall_good(3))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(4)/srcall_good(4))-1.0_WP) < tolerance )
    ASSERT(abs((srcall(5)/srcall_good(5))-1.0_WP) < tolerance )
    do is=1, ns
      ASSERT(abs((srcall(5+is)/srcall_good(5+is))-1.0_WP) < tolerance )
    end do
#endif

    ! Deallocate workspace
    call largo_BL_spatiotemporal_deallocate (workspace)

    call testframework_teardown()

end program bl_spatiotemporal_baseflow_f
