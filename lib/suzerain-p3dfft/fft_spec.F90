! This file is part of P3DFFT library
!
! Version 2.3
!
! Copyright (C) 2006-2008 Dmitry Pekurovsky
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!----------------------------------------------------------------------------

      module fft_spec

#ifdef FFTW
      include "fftw3.f"
      integer*8 :: plan1_frc = 0
      integer*8 :: plan1_bcr = 0
      integer*8 :: plan1_fc  = 0
      integer*8 :: plan2_fc  = 0
      integer*8 :: plan1_bc  = 0
      integer*8 :: plan2_bc  = 0

!     Default rigor set at compile time but may be modified at runtime
#if   defined ESTIMATE
      integer, public :: fftw_rigor = FFTW_ESTIMATE
#elif defined PATIENT
      integer, public :: fftw_rigor = FFTW_PATIENT
#else
      integer, public :: fftw_rigor = FFTW_MEASURE
#endif
      integer, parameter :: NULL = 0

#endif

#ifdef ESSL
      integer :: cnaux,rnaux1,rnaux2
      real*8,save,allocatable :: caux1(:),caux2(:),raux1(:),raux2(:)
      real*8,save :: raux3(1)
#endif

      end module
