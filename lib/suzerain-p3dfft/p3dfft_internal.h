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
#ifdef  __P3DFFT_INTERNAL_H____
#error "p3dfft_internal.h must not be multiply included"
#else
#define __P3DFFT_INTERNAL_H____

#include <assert.h>
#include <stdlib.h>

#ifndef P3DFFT_PRECISION
#error "P3DFFT_PRECISION must be #defined prior to including p3dfft_internal.h"
#endif


#ifdef IBM

#define FORT_MOD_NAME(NAME) __p3dfft_NMOD_##NAME
#define FORTNAME(NAME) NAME

#elif defined PGI

#define FORT_MOD_NAME(NAME) p3dfft_##NAME##_
#define FORTNAME(NAME) NAME##_

#elif defined __INTEL_COMPILER

#define FORT_MOD_NAME(NAME) p3dfft_mp_##NAME##_
#define FORTNAME(NAME) NAME##_

#elif defined __GNUC__

#define FORT_MOD_NAME(NAME) __p3dfft_MOD_##NAME
#define FORTNAME(NAME) NAME

#else

#error "Unknown compiler vendor: check p3dfft_internal.h for details"

#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Runtime access to compile time options */
extern void FORT_MOD_NAME(p3dfft_get_blocksizes)(int *nbl_x, int *nbl_y, int *nb1, int*nb);
extern void FORT_MOD_NAME(p3dfft_using_stride1)(int *used);
extern void FORT_MOD_NAME(p3dfft_get_precision)(int *prec);

/* Runtime logic */
extern void FORT_MOD_NAME(p3dfft_fftw_rigor)(int *);
extern void FORT_MOD_NAME(p3dfft_setup)(int *dims, int *nx, int *ny, int *nz, int *ow);
extern void FORT_MOD_NAME(p3dfft_get_dims)(int *, int *, int *, int *);
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B);
extern void FORT_MOD_NAME(ftran_r2c)(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B);
extern void FORT_MOD_NAME(btran_c2r)(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B);
extern void FORT_MOD_NAME(p3dfft_clean)();

/* GNU already includes an abort function */
#if defined __INTEL_COMPILER || !defined  __GNUC__
extern void FORTNAME(abort)();
#endif

inline
void p3dfft_get_blocksizes(int *nbl_x, int *nbl_y, int *nb1, int *nb)
{
    FORT_MOD_NAME(p3dfft_get_blocksizes)(nbl_x, nbl_y, nb1, nb);
}

inline
void p3dfft_fftw_rigor(unsigned rigor)
{
    int int_rigor = (int) rigor;
    FORT_MOD_NAME(p3dfft_fftw_rigor)(&int_rigor);
}

inline
int p3dfft_using_stride1()
{
    int using_stride1;
    FORT_MOD_NAME(p3dfft_using_stride1)(&using_stride1);
    return using_stride1;
}

inline
int p3dfft_get_precision()
{
    int prec;
    FORT_MOD_NAME(p3dfft_get_precision)(&prec);
    return prec;
}

inline
void p3dfft_setup(int *dims, int nx, int ny, int nz, int overwrite)
{
    /* Ensure library precision matches included header */
    if (sizeof(P3DFFT_PRECISION) == sizeof(float)) {
        assert(p3dfft_get_precision() == 1);
    } else if (sizeof(P3DFFT_PRECISION) == sizeof(double)) {
        assert(p3dfft_get_precision() == 2);
    } else {
        assert(0 /* Sanity failure */);
    }

    FORT_MOD_NAME(p3dfft_setup)(dims, &nx, &ny, &nz, &overwrite);
}

inline
void p3dfft_get_dims(int *start, int *end, int *size, int conf)
{
    FORT_MOD_NAME(p3dfft_get_dims)(start, end, size, &conf);
}

inline
void p3dfft_ftran_r2c(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B)
{
    FORT_MOD_NAME(p3dfft_ftran_r2c)(A, B);
}

inline
void p3dfft_btran_c2r(P3DFFT_PRECISION *A, P3DFFT_PRECISION *B)
{
    FORT_MOD_NAME(p3dfft_btran_c2r)(A, B);
}

inline
void p3dfft_clean()
{
    FORT_MOD_NAME(p3dfft_clean)();
}

/* GNU already includes an abort function */
#if defined __INTEL_COMPILER || !defined  __GNUC__
inline
void FORTNAME(abort)()
{
    abort();
}
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __P3DFFT_INTERNAL_H____ */
