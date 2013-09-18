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
!! $Id$

module testframework

#ifdef __INTEL_COMPILER
  use ifport, only: abort
#endif /* __INTEL_COMPILER */

  implicit none  ! Nothing implicit
  private        ! Everything default private

  logical, public :: verbose = .false.

  public :: testframework_setup
  public :: testframework_teardown
  public :: testframework_assert

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_setup (testsource)

    character(len=*), intent(in) :: testsource
    character(LEN=255)           :: buf

!   Check if assertions should be verbose in nature
    call get_environment_variable("LARGO_TEST_VERBOSE", buf)
    verbose = len_trim(buf) > 0

  end subroutine testframework_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_teardown ()

!   Nothing

  end subroutine testframework_teardown

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_assert (condition, sourcefile, line)

    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

    logical,          intent(in) :: condition
    character(len=*), intent(in) :: sourcefile
    integer,          intent(in) :: line

! Really?!  http://software.intel.com/en-us/forums/showthread.php?t=66887
#ifndef __INTEL_COMPILER
    flush (output_unit)
    flush (error_unit)
#endif
    if (verbose .and. condition) then
      write (output_unit, *) "Assertion passed at ", sourcefile, ":", line
    else if (.not. condition) then
      write (error_unit, *) "Assertion failed at ", sourcefile, ":", line
      stop
    end if
#ifndef __INTEL_COMPILER
    flush (error_unit)
    flush (output_unit)
#endif

  end subroutine testframework_assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module testframework
