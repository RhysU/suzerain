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
!! $Id: bl_temporal_f.f90 40608 2013-08-06 23:12:28Z topalian $

#include "testframework_assert.h"

program bl_temporal_baseflow_rans_f

    use largo
    use testframework

    use, intrinsic :: iso_c_binding, only: c_associated,   &
                                           WP => c_double, &
                                           c_int
    implicit none

    type(largo_workspace_ptr)         :: workspace
    real(WP), parameter               :: y       = 1.0_WP/ 10.0_WP
    real(WP), parameter               :: grDelta = 5.0_WP/100.0_WP

    character(len=255), parameter     :: turbmodel = "k_epsilon"
    integer(c_int), parameter         :: ntvar     = 2
    integer(c_int), parameter         :: ns        = 2
    integer(c_int), parameter         :: neq       = 5 + ns + ntvar
    integer(c_int), parameter         :: nvar      = neq
    integer(c_int), parameter         :: nvar_base = 5 + ns


    real(WP), dimension(nvar), parameter :: &
      mean    = (/                &
      &        1.0_WP/ 100.0_WP,  &
      &       45.0_WP/  10.0_WP,  &
      &        1.0_WP/1000.0_WP,  &
      &        5.0_WP/ 100.0_WP,  &
      &    41200.0_WP          ,  &
      &        2.0_WP/1000.0_WP,  &
      &        1.0_WP/1000.0_WP,  &
      &        5.0_WP/ 100.0_WP,  &
      &        3.0_WP/ 100.0_WP   &
      /)

    real(WP), dimension(nvar), parameter :: &
      dmean   = (/                &
      &         1.0_WP/ 10.0_WP,  &
      &        45.0_WP         ,  &
      &         1.0_WP/100.0_WP,  &
      &         5.0_WP/ 10.0_WP,  &
      &    412000.0_WP         ,  &
      &         2.0_WP/100.0_WP,  &
      &         1.0_WP/100.0_WP,  &
      &         5.0_WP/ 10.0_WP,  &
      &         3.0_WP/ 10.0_WP   &
      /)

    real(WP), dimension(nvar), parameter :: &
      grDA    = (/                 &
      &       2.0_WP/  100.0_WP,   &
      &       1.0_WP/   10.0_WP,   &
      &       1.0_WP/  100.0_WP,   &
      &       3.0_WP/  100.0_WP,   &
      &       4.0_WP/  100.0_WP,   &
      &       1.0_WP/ 1000.0_WP,   &
      &       2.0_WP/ 1000.0_WP,   &
      &       4.0_WP/  100.0_WP,   &
      &       2.0_WP/  100.0_WP    &

      /)

    real(WP), dimension(nvar), parameter :: &
      grDArms = (/                  &
      &       1.0_WP/  1000.0_WP,   &
      &       5.0_WP/  1000.0_WP,   &
      &       5.0_WP/ 10000.0_WP,   &
      &       2.0_WP/  1000.0_WP,   &
      &       1.0_WP/  1000.0_WP,   &
      &       1.0_WP/100000.0_WP,   &
      &       2.0_WP/100000.0_WP,   &
      &       0.0_WP            ,   &
      &       0.0_WP                &
      /)

    real(WP), dimension(nvar_base), parameter :: &
      base    = (/                  &
      &         5.0_WP/ 1000.0_WP,  &
      &        25.0_WP/   10.0_WP,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/  100.0_WP,  &
      &     20000.0_WP           ,  &
      &         1.0_WP/ 1000.0_WP,  &
      &         5.0_WP/10000.0_WP   &
      /)

#ifndef BASEFLOW_UNIFORM
    real(WP), dimension(nvar_base), parameter :: &
      dybase  = (/                &
      &         5.0_WP/ 100.0_WP,  &
      &        25.0_WP          ,  &
      &         5.0_WP/1000.0_WP,  &
      &         2.0_WP/  10.0_WP,  &
      &    200000.0_WP          ,  &
      &         1.0_WP/ 100.0_WP,  &
      &         5.0_WP/1000.0_WP   &
      /)

    real(WP), dimension(nvar_base), parameter :: &
      dtbase  = (/                &
      &         2.0_WP/ 1000.0_WP,  &
      &        15.0_WP/   10.0_WP,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/  100.0_WP,  &
      &     10000.0_WP           ,  &
      &         5.0_WP/10000.0_WP,  &
      &         2.0_WP/10000.0_WP   &
      /)

    real(WP), dimension(nvar_base)            :: &
      dxbase  = 0.0_WP
#else
    real(WP), dimension(nvar_base), parameter :: &
        dybase = 0.0_WP &
      , dtbase = 0.0_WP &
      , dxbase = 0.0_WP
#endif

    real(WP), dimension(neq), parameter :: &
      srcbase = (/                &
      &         1.0_WP/ 1000.0_WP,  &
      &         5.0_WP/   10.0_WP,  &
      &         1.0_WP/10000.0_WP,  &
      &         5.0_WP/ 1000.0_WP,  &
      &      5000.0_WP           ,  &
      &         2.0_WP/10000.0_WP,  &
      &         1.0_WP/10000.0_WP,  &
      &         0.0_WP           ,  &
      &         0.0_WP              &
      /)

    real(WP), dimension(neq)            :: srcmean

    real(WP), dimension(neq), parameter :: &
      srcmean_good = (/              &
      & -    17.0_WP/  20000.0_WP,   &
      & -    11.0_WP/     10.0_WP,   &
      & -     1.0_WP/  12500.0_WP,   &
      & -    11.0_WP/   2500.0_WP,   &
      & -  4788.0_WP             ,   &
      & -   251.0_WP/1000000.0_WP,   &
      & -    19.0_WP/ 250000.0_WP,   &
      &       1.0_WP/   2000.0_WP,   &
      &       9.0_WP/  10000.0_WP    &
      /)

    real(WP), parameter :: tolerance = 1.0E-14

    integer(c_int) :: is, it

    call testframework_setup(__FILE__)

    ! Initialize srcmean and srcrms
    srcmean = 0.0_WP

    ! Allocate workspace
    call largo_BL_temporal_allocate (workspace, neq, ns, ntvar, turbmodel)

    ! Init growth rates
    call largo_BL_temporal_init_rans (workspace, grDelta, grDA, grDArms)

    ! Compute prestep values
    call largo_BL_temporal_preStep_baseflow  (workspace,   base,  dybase,  &
                                                 dtbase, dxbase, srcbase)
    call largo_BL_temporal_preStep_sEtaMean_rans (workspace, y, mean, dmean)

    ! Compute sources
    call largo_BL_temporal_sEta_rans (workspace, 0.0_WP, 1.0_WP, srcmean(1))

    ! Check mean
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
    do it=1, ntvar
      ASSERT(abs((srcmean(5+ns+it)/srcmean_good(5+ns+it))-1.0_WP) < tolerance )
    end do
#endif

    ! Deallocate workspace
    call largo_BL_temporal_deallocate (workspace)

    call testframework_teardown()

end program bl_temporal_baseflow_rans_f
