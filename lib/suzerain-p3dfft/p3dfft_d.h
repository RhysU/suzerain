/*
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
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!----------------------------------------------------------------------------
*/
#ifndef __P3DFFT_S_H____
#define __P3DFFT_S_H____

#ifdef P3DFFT_PRECISION
#error "P3DFFT_PRECISION #define collision"
#endif

#define P3DFFT_PRECISION double
#include "p3dfft_internal.h"
#undef P3DFFT_PRECISION

#endif /* __P3DFFT_S_H____ */
