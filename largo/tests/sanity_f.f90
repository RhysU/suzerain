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
!! $Id: sanity_f.f90 41305 2013-09-11 18:37:35Z rhys $

#include "testframework_assert.h"

program sanity_f

    use testframework

    implicit none

    call testframework_setup(__FILE__)

    ! Automake's XFAIL_TESTS in Makefile.am is used to ensure
    ! assertion failures return non-zero exit codes to the OS.
    ! Otherwise, our Fortran-based test suite is wholly useless.
    ASSERT(1 == 0)

    call testframework_teardown()

end program sanity_f
