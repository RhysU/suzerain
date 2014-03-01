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

      module p3dfft

      implicit none

      include 'mpif.h'

      private

! Set precision

#ifndef SINGLE_PREC
       integer, parameter,public :: mytype=8
       integer, parameter,public:: mpireal = MPI_DOUBLE_PRECISION
       integer,parameter,public:: mpicomplex = MPI_DOUBLE_COMPLEX
#else
       integer, parameter,public :: mytype=4
       integer, parameter,public:: mpireal = MPI_REAL
       integer,parameter,public:: mpicomplex = MPI_COMPLEX
#endif

! global variables

      integer, public :: padd
      real(8), save,public :: timers(10)

      integer,save :: NX_fft,NY_fft,NZ_fft,numtasks,iproc,jproc
      integer,save :: ipid,jpid,taskid
      integer,save :: iistart,iiend,iisize,jjstart,jjsize,jjend
      integer,save ::jistart,kjstart,jisize,kjsize,jiend,kjend

      integer,save ::  nxh,nxhp

! mpi process info
!
      logical :: mpi_set=.false.
      integer, save :: mpi_comm_cart
      integer, save :: mpi_comm_row, mpi_comm_col
      integer,save, dimension(:), allocatable :: iist,iien,iisz
      integer,save, dimension(:), allocatable :: jist,jien,jisz
      integer,save, dimension(:), allocatable :: jjst,jjen,jjsz
      integer,save, dimension(:), allocatable :: kjst,kjen,kjsz

! mpi derived data types for implementing alltoallv using send-recvs
      integer,save,dimension(:),allocatable:: IfSndCnts,IfSndStrt
      integer,save,dimension(:),allocatable:: IfRcvCnts,IfRcvStrt
      integer,save,dimension(:),allocatable:: KfSndCnts,KfSndStrt
      integer,save,dimension(:),allocatable:: KfRcvCnts,KfRcvStrt
      integer,save,dimension(:),allocatable:: JrSndCnts,JrSndStrt
      integer,save,dimension(:),allocatable:: JrRcvCnts,JrRcvStrt
      integer,save,dimension(:),allocatable:: KrSndCnts,KrSndStrt
      integer,save,dimension(:),allocatable:: KrRcvCnts,KrRcvStrt
      integer,save,dimension(:,:),allocatable:: status
      complex(mytype), save, allocatable :: buf(:),buf1(:),buf2(:)
      logical :: OW = .false.
#ifdef USE_EVEN
      integer*8,save :: IfCntMax,KfCntMax
      logical KfCntUneven
#endif

#ifdef STRIDE1
#ifdef NBL_X
      integer, parameter :: NB1 = NBL_X
#else
      integer, parameter :: NB1=4
#endif
#ifdef NBL_Y
      integer, parameter :: NB=NBL_Y
#else
      integer, parameter :: NB=1
#endif
#endif

! ghost-cell support using p3dfft_init_ghosts(), p3dfft_update_ghosts()
      logical :: ghosts_set=.false.
      integer,save, dimension(:),   allocatable :: proc_id2coords
      integer,save, dimension(:,:), allocatable :: proc_coords2id
      integer,save, dimension(:), allocatable :: gneighb_r,gneighb_c
      real(kind=mytype),save, dimension(:), allocatable :: gbuf_snd,gbuf_recv
      integer,save :: goverlap
      integer,save :: gmess_rsize(3), gmess_csize(3)
      integer,save :: gproc_rsize(3), gproc_csize(3)
      integer,save :: gproc_rstart(3), gproc_cstart(3)
      integer,save :: gproc_rend(3), gproc_cend(3)
      integer,save :: gslab_rsize(3,3), gslab_rstart(3,6), gslab_rend(3,6)
      integer,save :: gslab_csize(3,3), gslab_cstart(3,6), gslab_cend(3,6)
      integer(kind=8),save :: gmem_rstart(6), gmem_cstart(6)

      public :: p3dfft_btran_c2r
      public :: p3dfft_clean
      public :: p3dfft_fftw_rigor
      public :: p3dfft_ftran_r2c
      public :: p3dfft_get_dims
      public :: p3dfft_get_precision
      public :: p3dfft_setup
      public :: p3dfft_using_stride1

!-------------------
      contains
!-------------------

!=========================================================
! Return compile-time blocking options
!
      subroutine p3dfft_get_blocksizes(nbl_x,nbl_y,nb1,nb)
!=========================================================
      integer, intent(out) :: nbl_x, nbl_y, nb1, nb

      nbl_x = NBL_X
      nbl_y = NBL_Y
      nb1   = NB1
      nb    = NB

      end subroutine p3dfft_get_blocksizes

!=======================================================
! Return whether or not STRIDE1 was used at compile-time
! Use of {0,1} for result to simplify C interoperability
      subroutine p3dfft_using_stride1(using_stride1)
!=======================================================
      integer, intent(out) :: using_stride1

#ifdef STRIDE1
      using_stride1 = 1
#else
      using_stride1 = 0
#endif

      end subroutine p3dfft_using_stride1

!=======================================================
! Set the FFTW rigor flags to be used (e.g. FFTW_ESTIMATE).
! Must be one prior to p3dfft_setup invocation.  No effect
! when FFTW not used.
      subroutine p3dfft_fftw_rigor(flag)
!=======================================================
      use fft_spec
      implicit none
      integer, intent(in) :: flag

#ifdef FFTW
      fftw_rigor = flag
#endif

      end subroutine p3dfft_fftw_rigor

!========================================================
! Return precision used at compile-time with 1 indicating
! single precision and 2 indicating double precision.
      subroutine p3dfft_get_precision(prec)
!=======================================================
      integer, intent(out) :: prec

#ifdef SINGLE_PREC
      prec = 1
#else
      prec = 2
#endif

      end subroutine p3dfft_get_precision

!=====================================================
! Return array dimensions for either real-space (conf=1) or wavenumber-space(conf=2)
!
      subroutine p3dfft_get_dims(istart,iend,isize,conf)
!=====================================================

      integer istart(3),iend(3),isize(3),conf

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
      else

      if(conf .eq. 1) then
         istart(1) = 1
         iend(1) = NX_fft
         isize(1) = NX_fft
         istart(2) = jistart
         iend(2) = jiend
         isize(2) = jisize
         istart(3) = kjstart
         iend(3) = kjend
         isize(3) = kjsize
      else if(conf .eq. 2) then
#ifdef STRIDE1
         istart(2) = iistart
         iend(2) = iiend
         isize(2) = iisize
         istart(3) = jjstart
         iend(3) = jjend
         isize(3) = jjsize
         istart(1) = 1
         iend(1) = NZ_fft
         isize(1) = NZ_fft
#else
         istart(1) = iistart
         iend(1) = iiend
         isize(1) = iisize
         istart(2) = jjstart
         iend(2) = jjend
         isize(2) = jjsize
         istart(3) = 1
         iend(3) = NZ_fft
         isize(3) = NZ_fft
#endif
      endif

      endif
      end subroutine p3dfft_get_dims

! =========================================================
      subroutine p3dfft_setup(dims,nx,ny,nz,overwrite)
!========================================================

      implicit none

      integer i,j,k,nx,ny,nz,err
      integer ierr, dims(2),  cartid(2)
      logical periodic(2),remain_dims(2),overwrite
      integer impid, ippid, jmpid, jppid
      integer*8 nm,n1,n2,li
      real(mytype), allocatable :: R(:)

      if(nx .le. 0 .or. ny .le. 0 .or. nz .le. 0) then
         print *,'Invalid dimensions :',nx,ny,nz
         call abort
      endif

      if(mpi_set) then
         print *,'P3DFFT Setup error: the problem is already initialized. '
         print *,'Currently multiple setups not supported.'
         print *,'Quit the library using p3dfft_clean before initializing another setup'
         call abort
      endif

      OW = overwrite

      timers = 0.0

      mpi_set = .true.
      nx_fft = nx
      ny_fft = ny
      nz_fft = nz
      nxh=nx/2
      nxhp=nxh+1

      call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)

      if(dims(1) .le. 0 .or. dims(2) .le. 0 .or. &
           dims(1)*dims(2) .ne. numtasks) then
         print *,'Invalid processor geometry: ',dims,' for ',numtasks,'tasks'
         call abort
      endif

      iproc = dims(1)
      jproc = dims(2)

#ifndef DIMS_C
       i = dims(1)
       dims(1) = dims(2)
       dims(2) = i
#endif

      periodic(1) = .false.
      periodic(2) = .false.
! creating cartesian processor grid
      call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic, &
           .false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
! process with a linear id of 5 may have cartid of (3,1)

#ifdef DIMS_C
      ipid = cartid(1)
      jpid = cartid(2)
#else
      ipid = cartid(2)
      jpid = cartid(1)
#endif

! store processor-grid-informations
      cartid(1) = ipid
      cartid(2) = jpid
      allocate(proc_id2coords(0:(iproc*jproc)*2-1))
      call MPI_Allgather( cartid,        2, MPI_INTEGER, &
                          proc_id2coords, 2, MPI_INTEGER,&
                          mpi_comm_cart,ierr)
      allocate(proc_coords2id(0:iproc-1,0:jproc-1))
      do i=0,(iproc*jproc)-1
        proc_coords2id( proc_id2coords(2*i), &
                                        proc_id2coords(2*i+1)) = i
      enddo

! here i is east-west j is north-south
! impid is west neighbour ippid is east neighbour and so on
      impid = ipid - 1
      ippid = ipid + 1
      jmpid = jpid - 1
      jppid = jpid + 1
!boundary processes
      if (ipid.eq.0) impid = MPI_PROC_NULL
      if (jpid.eq.0) jmpid = MPI_PROC_NULL
      if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
      if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL
! using cart comworld create east-west(row) sub comworld

#ifdef DIMS_C
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
#else
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
#endif

      allocate (iist(0:iproc-1))
      allocate (iisz(0:iproc-1))
      allocate (iien(0:iproc-1))
      allocate (jjst(0:jproc-1))
      allocate (jjsz(0:jproc-1))
      allocate (jjen(0:jproc-1))
      allocate (jist(0:iproc-1))
      allocate (jisz(0:iproc-1))
      allocate (jien(0:iproc-1))
      allocate (kjst(0:jproc-1))
      allocate (kjsz(0:jproc-1))
      allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny,nz) => (iproc,jproc)
!
      call MapDataToProc(nxhp,iproc,iist,iien,iisz)
      call MapDataToProc(ny,iproc,jist,jien,jisz)
      call MapDataToProc(ny,jproc,jjst,jjen,jjsz)
      call MapDataToProc(nz,jproc,kjst,kjen,kjsz)

! These are local array indices for each processor

      iistart = iist(ipid)
      jjstart = jjst(jpid)
      jistart = jist(ipid)
      kjstart = kjst(jpid)
      iisize= iisz(ipid)
      jjsize= jjsz(jpid)
      jisize= jisz(ipid)
      kjsize= kjsz(jpid)
      iiend = iien(ipid)
      jjend = jjen(jpid)
      jiend = jien(ipid)
      kjend = kjen(jpid)


#ifdef USE_EVEN
      IfCntMax = iisz(iproc-1)*jisz(iproc-1)*kjsize*mytype*2
      KfCntMax = iisize * jjsz(jproc-1) * kjsz(jproc-1)*mytype*2
      if(mod(ny,jproc) .ne. 0 .or. mod(nz,jproc) .ne. 0) then
         KfCntUneven = .true.
      else
         KfCntUneven = .false.
      endif
#endif

! We may need to pad arrays due to uneven size
      padd = max(iisize*jjsize*nz_fft,iisize*ny_fft*kjsize) - nxhp*jisize*kjsize
      if(padd .le. 0) then
         padd=0
      else
         if(mod(padd,nxhp*jisize) .eq. 0) then
            padd = padd / (nxhp*jisize)
         else
            padd = padd / (nxhp*jisize)+1
         endif

      endif

! Initialize FFTW and allocate buffers for communication
      nm = nxhp * jisize * (kjsize+padd)
      allocate(buf1(nm),stat=err)
      if(err .ne. 0) then
         print *,'p3dfft_setup: Error allocating buf1 (',nm
      endif
      allocate(R(nm*2),stat=err)
      if(err .ne. 0) then
         print *,'p3dfft_setup: Error allocating R (',nm*2
      endif
      buf1 = 0.0
      R = 0.0
      call init_fft(buf1,R,nm)
      deallocate(R)

#ifdef USE_EVEN
      n1 = IfCntMax * iproc /(mytype*2)
      n2 = KfCntMax * jproc / (mytype*2)
      n1 = max(n1,n2)
      if(n1 .gt. nm) then
         deallocate(buf1)
         allocate(buf1(n1))
      endif
      n1 = max(n1,nm)
      allocate(buf2(n1))
      buf1 = 0.0
#else
      allocate(buf2(nm))
#endif

      allocate(buf(nm))
      buf = 0.0

! Displacements and buffer counts for mpi_alltoallv

      allocate (IfSndStrt(0:iproc-1))
      allocate (IfSndCnts(0:iproc-1))
      allocate (IfRcvStrt(0:iproc-1))
      allocate (IfRcvCnts(0:iproc-1))

      allocate (KfSndStrt(0:jproc-1))
      allocate (KfSndCnts(0:jproc-1))
      allocate (KfRcvStrt(0:jproc-1))
      allocate (KfRcvCnts(0:jproc-1))

      allocate (JrSndStrt(0:jproc-1))
      allocate (JrSndCnts(0:jproc-1))
      allocate (JrRcvStrt(0:jproc-1))
      allocate (JrRcvCnts(0:jproc-1))

      allocate (KrSndStrt(0:iproc-1))
      allocate (KrSndCnts(0:iproc-1))
      allocate (KrRcvStrt(0:iproc-1))
      allocate (KrRcvCnts(0:iproc-1))


!   start pointers and types of send  for the 1st forward transpose
      do i=0,iproc-1
         IfSndStrt(i) = (iist(i) -1)* jisize*kjsize*mytype*2
         IfSndCnts(i) = iisz(i) * jisize*kjsize*mytype*2

!   start pointers and types of recv for the 1st forward transpose
         IfRcvStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         IfRcvCnts(i) = jisz(i) * iisize*kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd forward transpose
      do i=0,jproc-1
         KfSndStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         KfSndCnts(i) = iisize*kjsize*jjsz(i)*mytype*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         KfRcvCnts(i) = iisize*jjsize*kjsz(i)*mytype*2
      end do

!   start pointers and types of send  for the 1st inverse transpose
      do i=0,jproc-1
         JrSndStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         JrSndCnts(i) = iisize*jjsize*kjsz(i)*mytype*2

!   start pointers and types of recv for the 1st inverse transpose
         JrRcvStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         JrRcvCnts(i) = jjsz(i) * iisize * kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd inverse transpose
      do i=0,iproc-1
         KrSndStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         KrSndCnts(i) = jisz(i) * iisize*kjsize*mytype*2

!   start pointers and types of recv for the 2nd inverse transpose
         KrRcvStrt(i) = (iist(i) -1) * jisize*kjsize*mytype*2
         KrRcvCnts(i) = jisize*iisz(i)*kjsize*mytype*2
      enddo

      end subroutine p3dfft_setup

!==================================================================
      subroutine MapDataToProc (data,proc,st,en,sz)
!========================================================
!
       implicit none
       integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
       integer i,size,nl,nu

       size=data/proc
       nu = data - size * proc
       nl = proc - nu
       st(0) = 1
       sz(0) = size
       en(0) = size
       do i=1,nl-1
         st(i) = st(i-1) + size
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      size = size + 1
      do i=nl,proc-1
         st(i) = en(i-1) + 1
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      en(proc-1)= data
      sz(proc-1)= data-st(proc-1)+1

      end subroutine


!========================================================

      subroutine init_fft(A,B,n1)

      use fft_spec
      implicit none

      integer*8 n1
      complex(mytype) A(n1)
      real(mytype) B(n1*2)


      call init_work(nx_fft,ny_fft,nz_fft)
      call plan_f_r2c(B,nx_fft,A,nxhp,nx_fft,jisize*kjsize)
      call plan_b_c2r(A,nxhp,B,nx_fft,nx_fft,jisize*kjsize)
#ifdef STRIDE1
      call plan_f_c1(A,iisize,1,A,iisize,1,ny_fft,iisize)
      call plan_b_c1(A,iisize,1,A,iisize,1,ny_fft,iisize)
         call plan_f_c2(A,1,nz_fft, &
           A,1,nz_fft,nz_fft,NB*iisize,.false.)
      if(OW) then
         call plan_b_c2(A,1,nz_fft, &
           A,1,nz_fft,nz_fft,NB*iisize,.false.)
      else
         call plan_b_c2(A,1,nz_fft, &
           B,1,nz_fft,nz_fft,NB*iisize,.false.)
      endif
#else
      call plan_f_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_b_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_f_c2(A,iisize*jjsize, 1, &
           A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      if(OW) then
         call plan_b_c2(A,iisize*jjsize, 1, &
           A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      else
         call plan_b_c2(A,iisize*jjsize, 1, &
           B,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      endif
#endif

      return
      end subroutine

!========================================================
!  3D FFT inverse transform with 2D domain decomposition
!
!  The order of array elements in memory is the same in all stages: (x,y,z)

! Input: XYZg - comlpex array, with entire Z dimension local,
!               while x and y are block-distributed among processors in 2D grid
! Output: XgYZ - an array of real, x dimension is entirely local within
! processors memory  while y and z are block-distributed among processors
! in 2D grid
!
! The arrays may occupy the same memory space
! In this case their first elements should coincide.
! Naturally, output overwrites input
!

      subroutine p3dfft_btran_c2r (XYZg,XgYZ)
!========================================================

      use fft_spec
      implicit none

      real(mytype),TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif

      integer x,y,z,i,k,nx,ny,nz
      integer(8) Nl

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif
      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

! FFT Tranform (C2C) in Z for all x and y

      if(jproc .gt. 1) then

#ifdef STRIDE1
         call init_b_c(buf, 1,nz, &
              buf, 1, nz,nz,NB*iisize)
         call bcomm1_trans(XYZg,buf2,buf,timers(3),timers(10))
#else

         if(OW) then

            if(iisize*jjsize .gt. 0) then
!SUZERAIN      call init_b_c(XYZg, iisize*jjsize, 1, &
!                   XYZg, iisize*jjsize, 1,nz,iisize*jjsize)

               timers(10) = timers(10) - MPI_Wtime()
!SUZERAIN      call exec_b_c2(XYZg, iisize*jjsize,1, &
!                   XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(10) = timers(10) + MPI_Wtime()
            endif
            call bcomm1(XYZg,buf,timers(3),timers(10))

         else
            if(iisize*jjsize .gt. 0) then
!SUZERAIN      call init_b_c(buf, iisize*jjsize, 1, &
!                buf, iisize*jjsize, 1,nz,iisize*jjsize)

               timers(10) = timers(10) - MPI_Wtime()
!SUZERAIN      call exec_b_c2(XYZg, iisize*jjsize,1, &
!                buf, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(10) = timers(10) + MPI_Wtime()
               call bcomm1(buf,buf,timers(3),timers(10))
            endif
         endif

#endif


      else
#ifdef STRIDE1
         call reorder_trans_b(XYZg,buf,buf2)
#else
         if(OW) then
            Nl = iisize*jjsize*nz
!SUZERAIN   call init_b_c(XYZg, iisize*jjsize, 1, &
!             XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
!SUZERAIN   call exec_b_c2(XYZg, iisize*jjsize, 1, &
!             XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
            call ar_copy(XYZg,buf,Nl)
         else
!SUZERAIN   call init_b_c(XYZg, iisize*jjsize, 1, &
!             buf, iisize*jjsize, 1,nz,iisize*jjsize)
!SUZERAIN   call exec_b_c2(XYZg, iisize*jjsize, 1, &
!             buf, iisize*jjsize, 1,nz,iisize*jjsize)
         endif
#endif
      endif

! Exhange in columns if needed

!
! FFT Transform (C2C) in y dimension for all x, one z-plane at a time
!
      if(iisize * kjsize .gt. 0) then

         call init_b_c(buf,iisize,1,buf,iisize,1,ny,iisize)

         timers(9) = timers(9) - MPI_Wtime()
         do z=kjstart,kjend

            call btran_y_zplane(buf,z-kjstart,iisize,kjsize,iisize,1, &
                 buf,z-kjstart,iisize,kjsize,iisize,1,ny,iisize)

         enddo
         timers(9) = timers(9) + MPI_Wtime()
      endif

      if(iproc .gt. 1) call bcomm2(buf,buf,timers(4))

! Perform Complex-to-real FFT in x dimension for all y and z
      if(jisize * kjsize .gt. 0) then

         call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(8) = timers(8) - MPI_Wtime()
         call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(8) = timers(8) + MPI_Wtime()
      endif

      return
      end subroutine

      subroutine wrap_exec_b_c2(A,strideA,B,strideB, &
           N,m,L,k)

      complex(mytype) A(L,N),B(L,N)
      integer strideA,strideB,N,m,L,k


      call exec_b_c2(A(k,1),strideA,1,B(k,1),strideB,1,N,m)

      return
      end subroutine


!========================================================
!  3D FFT Forward transform with 2D domain decomposition
!
!  The order of array elements in memory is the same in all stages: (x,y,z)

! Input: XgYZ - an array of real, x dimension is entirely local within
! processors memory  while y and z are block-distributed among processors
! in 2D grid
! Output: XYZg - comlpex array, with entire Z dimension local,
!               while x and y are block-distributed among processors in 2D grid
!
! The arrays may occupy the same memory space
! In this case their first elements should coincide.
! Naturally, output overwrites input
!
      subroutine p3dfft_ftran_r2c (XgYZ,XYZg)
!========================================================

      use fft_spec
      implicit none

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif
      integer x,y,z,i,err,nx,ny,nz
      integer(8) Nl

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft
!
! FFT transform (R2C) in X for all z and y
!
      if(jisize * kjsize .gt. 0) then
         call init_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) - MPI_Wtime()
         call exec_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)
         timers(5) = timers(5) + MPI_Wtime()
      endif

! Exchange data in rows
      if(iproc .gt. 1) call fcomm1(buf,buf,timers(1))

! FFT transform (C2C) in Y for all x and z, one Z plane at a time

      if(iisize * kjsize .gt. 0) then
         call init_f_c(buf,iisize,1,buf,iisize,1,ny,iisize)
         timers(6) = timers(6) - MPI_Wtime()
         do z=1,kjsize
            call ftran_y_zplane(buf,z-1,iisize,kjsize,iisize,1, &
                 buf,z-1,iisize,kjsize,iisize,1,ny,iisize)
         enddo
         timers(6) = timers(6) + MPI_Wtime()

      endif

! Exchange data in columns
      if(jproc .gt. 1) then
#ifdef STRIDE1
! For stride1 option combine second transpose with transform in Z
         call init_f_c(buf,1,nz, buf,1,nz,nz,iisize*NB)
         call fcomm2_trans(buf,XYZg,timers(2),timers(7))
#else


! FFT Transform (C2C) in Z for all x and y

! Transpose y-z

         call fcomm2(buf,XYZg,timers(2),timers(7))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is
! faster than done in-place

         if(iisize * jjsize .gt. 0) then
!SUZERAIN   call init_f_c(XYZg,iisize*jjsize, 1, &
!                XYZg,iisize*jjsize, 1,nz,iisize*jjsize)

            timers(7) = timers(7) - MPI_Wtime()
!SUZERAIN   call exec_f_c2(XYZg,iisize*jjsize, 1, &
!                XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
            timers(7) = timers(7) + MPI_Wtime()
         endif

#endif

      else
#ifdef STRIDE1

         call reorder_trans_f(buf,XYZg)
#else
         Nl = iisize*jjsize*nz
         call ar_copy(buf,XYZg,Nl)
!SUZERAIN call init_f_c(XYZg,iisize*jjsize, 1, &
!             XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
!SUZERAIN call exec_f_c2(XYZg,iisize*jjsize, 1, &
!             XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
#endif
      endif


      return
      end subroutine


#ifdef STRIDE1

! This routine is called only when jproc=1, and only when stride1 is used
! transform backward in Z and transpose array in memory

      subroutine reorder_trans_b(A,B,C)

      complex(mytype) B(iisize,ny_fft,nz_fft)
      complex(mytype) A(nz_fft,iisize,ny_fft)
      complex(mytype) C(nz_fft,iisize,NB*iisize*nz_fft)
      integer x,y,z,iy,iz,y2,z2

      if(OW) then
         do y=1,ny_fft,nb
            y2 = min(y+nb-1,ny_fft)
!SUZERAIN   call exec_b_c2(A(1,1,y),1,nz_fft,A(1,1,y),1,nz_fft,nz_fft,nb*iisize)
            do z=1,nz_fft,nb
               z2 = min(z+nb-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = A(iz,x,iy)
                     enddo
                  enddo
               enddo

            enddo
         enddo

      else
         do y=1,ny_fft,nb
            y2 = min(y+nb-1,ny_fft)
!SUZERAIN   call exec_b_c2(A(1,1,y),1,nz_fft,C(1,1,y),1,nz_fft,nz_fft,nb*iisize)
            do z=1,nz_fft,nb
               z2 = min(z+nb-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = C(iz,x,iy)
                     enddo
                  enddo
               enddo

            enddo
         enddo

      endif

      return
      end subroutine


! This routine is called only when jproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Z

      subroutine reorder_trans_f(A,B)

      complex(mytype) A(iisize,ny_fft,nz_fft)
      complex(mytype) B(nz_fft,iisize,ny_fft)
!      complex(mytype) C(nz_fft,iisize,ny_fft)
      integer x,y,z,iy,iz,y2,z2

      do y=1,ny_fft,nb
         y2 = min(y+nb-1,ny_fft)
         do z=1,nz_fft,nb
            z2 = min(z+nb-1,nz_fft)
            do iz=z,z2
               do iy=y,y2
                  do x=1,iisize
                     B(iz,x,iy) = A(x,iy,iz)
                  enddo
               enddo
            enddo
         enddo
!SUZERAIN call exec_f_c2(B(1,1,y),1,nz_fft,B(1,1,y),1,nz_fft,nz_fft,nb*iisize)
      enddo

      return
      end subroutine
#endif

! Communication Module
!
! Contains 4 routines for forward and backward exchanges inrows and columns
! of processors. Uses MPI_Alltoallv or MPI_Alltoall routine
!

!========================================================
! Transpose X and Y pencils

      subroutine fcomm1(source,dest,t)
!========================================================

      implicit none

      complex(mytype) source(nxhp,jisize,kjsize)
#ifdef STRIDE12
      complex(mytype) dest(ny_fft,iisize,kjsize)
#else
      complex(mytype) dest(iisize,ny_fft,kjsize)
#endif

      real(8) t
      integer x,y,i,ierr,z,xs,j,n,ix,iy
      integer*8 position,pos1

! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf

      position = 1

      do i=0,iproc-1
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  buf1(position) = source(x,y,z)
                  position = position +1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1)*IfCntMax/(mytype*2)+1
#endif
      enddo

#ifdef USE_EVEN

! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors)

      t = t - MPI_Wtime()
      call mpi_alltoall(buf1,IfCntMax, mpi_byte, &
           buf2,IfCntMax, mpi_byte,mpi_comm_row,ierr)
      t = t + MPI_Wtime()

#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
      t = t - MPI_Wtime()
      call mpi_alltoallv(buf1,IfSndCnts, IfSndStrt,mpi_byte, &
           buf2,IfRcvCnts, IfRcvStrt,mpi_byte,mpi_comm_row,ierr)
      t = MPI_Wtime() + t
#endif

! Unpack the data

#ifdef STRIDE12
      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            pos1 = position
            do y=jist(i),jien(i),nb1
               do x=1,iisize,nb1
                  do iy = y,min(y+nb1-1,jien(i))
                     position = pos1 + x-1 +(iy-jist(i))*iisize
                     do ix=x,min(x+nb1-1,iisize)
                        dest(iy,ix,z) = buf2(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            enddo
            position = pos1 + jisz(i) * iisize
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo

#else
      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            do y=jist(i),jien(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position + 1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo
#endif

      return
      end subroutine

!========================================================
! Transpose Y and Z pencils

      subroutine fcomm2(source,dest,t,tc)
!========================================================

      implicit none

#ifdef STRIDE1
      complex(mytype) source(iisize,ny_fft,kjsize)
      complex(mytype) dest(nz_fft,iisize,jjsize)
#else
      complex(mytype) source(iisize,ny_fft,kjsize)
      complex(mytype) dest(iisize,jjsize,nz_fft)
#endif
      real(8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz
      integer*8 position,pos1

! Pack send buffers for exchanging y and z for all x at once

      position = 1

      do i=0,jproc-1
         do z=1,kjsize
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo

#ifdef USE_EVEN
         position = position + (KfCntMax/(mytype*2) - jjsz(i)*iisize*kjsize)
#endif
      enddo

! Exchange y-z buffers in columns of processors

      t = t - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

      if(KfCntUneven) then

         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf2,KfCntMax, mpi_byte,mpi_comm_col,ierr)

         t = MPI_Wtime() + t

         tc = tc - MPI_Wtime()

         position = 1
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo

         tc = tc + MPI_Wtime()

      else

         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           dest,KfCntMax, mpi_byte,mpi_comm_col,ierr)
         t = MPI_Wtime() + t

      endif

#else
! Use MPI_Alltoallv

      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte, &
           dest,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = MPI_Wtime() + t

#endif
      return
      end subroutine

#ifdef STRIDE1

!========================================================
! Transpose Y and Z pencils
! Assume Stride1 data structure

      subroutine fcomm2_trans(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,ny_fft,kjsize)
!      complex(mytype) buf(nz_fft*iisize*jjsize)
      complex(mytype) dest(nz_fft,iisize,jjsize)

      real(8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,ix,x2,n
      integer*8 position,pos1

! Pack send buffers for exchanging y and z for all x at once

      position = 1

      do i=0,jproc-1
         do z=1,kjsize
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo

#ifdef USE_EVEN
         position = position + (KfCntMax/(mytype*2) - jjsz(i)*iisize*kjsize)
#endif
      enddo

! Exchange y-z buffers in columns of processors

      t = t - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

      call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf,KfCntMax, mpi_byte,mpi_comm_col,ierr)

      t = MPI_Wtime() + t
      tc = tc - MPI_Wtime()

      do y=1,jjsize,NB
         y2 = min(y+NB-1,jjsize)

         do x=1,iisize,NB1
            x2 = min(x+NB1-1,iisize)

            n = iisize - (x2-x)
            do i=0,jproc-1
               pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
               do z=kjst(i),kjen(i)
                  position = pos1
                  do iy=y,y2
                     do ix=x,x2
! Here we are sure that dest and buf are different
                        dest(z,ix,iy) = buf(position)
                        position = position +1
                     enddo
                     position = position + n
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
         enddo
!SUZERAIN call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
!              nz_fft,NB*iisize)

      enddo
      tc = tc + MPI_Wtime()

#else
! Use MPI_Alltoallv

      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte, &
           buf,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = MPI_Wtime() + t

      tc = tc - MPI_Wtime()

      if(NB .eq. 1) then

         do y=1,jjsize

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)

               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do ix=x,x2
! Here we are sure that dest and buf are different
                     dest(z,ix,y) = buf(position)
                     position = position +1
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
!SUZERAIN   call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
!                nz_fft,iisize)
         enddo

      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)

               n = iisize - (x2-x)
               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do iy=y,y2
!                  position = ((z-1)*jjsize +iy-1)*iisize + x
                     do ix=x,x2
! Here we are sure that dest and buf are different
                        dest(z,ix,iy) = buf(position)
                        position = position +1
                     enddo
                     position = position + n
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
!SUZERAIN   call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
!                nz_fft,NB*iisize)
         enddo

      endif
      tc = tc + MPI_Wtime()

#endif
      return
      end subroutine

#endif

!========================================================
! Transpose back Z to Y pencils

      subroutine bcomm1 (source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,jjsize,nz_fft)
      complex(mytype) dest(iisize,ny_fft,kjsize)
      real(8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,iz,y2,z2
      integer*8 position,pos1

!     Pack the data for sending


#ifdef USE_EVEN

      if(KfCntUneven) then
         tc = tc - MPI_Wtime()
         position = 1
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo
         tc = tc + MPI_Wtime()

         t = t - MPI_Wtime()
         call mpi_alltoall(buf1,KfCntMax,mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)

      else
         t = t - MPI_Wtime()
         call mpi_alltoall(source,KfCntMax,mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
      endif
#else

!     Exchange data in columns
      t = t - MPI_Wtime()
      call mpi_alltoallv(source,JrSndCnts, JrSndStrt,mpi_byte, &
           buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif


      t = t + MPI_Wtime()

! Unpack receive buffers into dest

      position=1
      do i=0,jproc-1
         do z=1,kjsize
#ifdef STRIDE12
            do x=1,iisize
               do y=jjst(i),jjen(i)
                  dest(y,x,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#else
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#endif
         enddo
#ifdef USE_EVEN
         position = (i+1)*KfCntMax/(mytype*2)+1
#endif
      enddo


      return
      end subroutine

#ifdef STRIDE1

!========================================================
! Transpose back Z to Y pencils
! Assumes stride1 data structure

      subroutine bcomm1_trans (source,buf3,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(nz_fft,iisize,jjsize)
      complex(mytype) buf3(nz_fft,iisize,NB)
      complex(mytype) dest(iisize,ny_fft,kjsize)

      real(8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n
      integer*8 position,pos1

!     Pack the data for sending


#ifdef USE_EVEN
! Use MPI_Alltoall

      tc = tc - MPI_Wtime()

      if(OW) then

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

!SUZERAIN   call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
!             nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)

               do i=0,jproc-1
                  pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
                  do z=kjst(i),kjen(i)
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = source(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + iisize - (x2-x)
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
            enddo
         enddo

         tc = tc + MPI_Wtime()
         t = t - MPI_Wtime()
         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
         t = t + MPI_Wtime()

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
! Here we are sure that buf is not the same as dest
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo

      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

!SUZERAIN   call exec_b_c2(source(1,1,y),1,nz_fft,buf3(1,1,y),1,nz_fft, &
!                nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)

               do i=0,jproc-1
                  pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
                  do z=kjst(i),kjen(i)
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = buf3(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + iisize - (x2-x)
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
            enddo

         enddo

         tc = tc + MPI_Wtime()
         t = t - MPI_Wtime()
         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)

         t = t + MPI_Wtime()

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo

      endif

#else
! Use MPI_Alltoallv

      tc = tc - MPI_Wtime()

      if(OW) then

         if(NB .eq.1) then

            do y=1,jjsize

!SUZERAIN      call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
!                 nz_fft,iisize)

               do x=1,iisize,NB1
                  x2 = min(x+NB1-1,iisize)

                  pos1 = (y-1)*iisize + x
                  do z=1,nz_fft
                     position = pos1
                     do ix=x,x2
                        buf1(position) = source(z,ix,y)
                        position = position+1
                     enddo
                     pos1 = pos1 + iisize*jjsize
                  enddo
               enddo
            enddo

         else

            do y=1,jjsize,NB
               y2 = min(y+NB-1,jjsize)

!SUZERAIN      call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
!                   nz_fft,NB*iisize)

               do x=1,iisize,NB1
                  x2 = min(x+NB1-1,iisize)

                  n = iisize - (x2-x)
                  pos1 = (y-1)*iisize + x
                  do z=1,nz_fft
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = source(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + n
                     enddo
                     pos1 = pos1 + iisize*jjsize
                  enddo
               enddo
            enddo

         endif

         tc = tc + MPI_Wtime()

         t = t - MPI_Wtime()
         call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, &
              buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
         t = t + MPI_Wtime()

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
! Here we are sure that buf is not the same as dest
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo

      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

!SUZERAIN   call exec_b_c2(source(1,1,y),1,nz_fft,buf3(1,1,y),1,nz_fft, &
!             nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)

               n = iisize - (x2-x)
               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do iy = y,y2
                     do ix=x,x2
                        buf1(position) = buf3(z,ix,iy)
                        position = position+1
                     enddo
                     position = position + n
                  enddo
                  pos1 = pos1 + iisize*jjsize
               enddo
            enddo
         enddo

         tc = tc + MPI_Wtime()

         t = t - MPI_Wtime()
         call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, &
              buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)

         t = t + MPI_Wtime()

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo

      endif
#endif


      return
      end subroutine

#endif

!========================================================
! Transpose back Y to X pencils

      subroutine bcomm2(source,dest,t)
!========================================================

      implicit none

      complex(mytype) dest(nxhp,jisize,kjsize)
      complex(mytype) source(iisize,ny_fft,kjsize)
      real(8) t
      integer x,y,z,i,ierr
      integer*8 position

! Pack and exchange x-z buffers in rows

      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            do y=jist(i),jien(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo

      t = t - MPI_Wtime()

#ifdef USE_EVEN
      call mpi_alltoall (buf1,IfCntMax, &
                 mpi_byte, buf2,IfCntMax,mpi_byte, &
                 mpi_comm_row,ierr)
#else
      call mpi_alltoallv (buf1,KrSndCnts, KrSndStrt, &
                 mpi_byte, buf2,KrRcvCnts,KrRcvStrt,mpi_byte, &
                 mpi_comm_row,ierr)
#endif
      t = t + MPI_Wtime()

! Unpack receive buffers into dest

      position=1
      do i=0,iproc-1
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  dest(x,y,z) = buf2(position)
                  position = position+1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo

      return
      end subroutine

!========================================================
      subroutine p3dfft_clean

!!--------------------------------------------------------------

      use fft_spec

! Clean-up FFT library
#ifdef FFTW
#ifndef SINGLE_PREC
      if (plan1_frc .ne. 0) then
          call dfftw_destroy_plan(plan1_frc)
          plan1_frc = 0
      endif
      if (plan1_bcr .ne. 0) then
          call dfftw_destroy_plan(plan1_bcr)
          plan1_bcr = 0
      endif
      if (plan1_fc .ne. 0) then
          call dfftw_destroy_plan(plan1_fc)
          plan1_fc = 0
      endif
      if (plan2_fc .ne. 0) then
          call dfftw_destroy_plan(plan2_fc)
          plan2_fc = 0
      endif
      if (plan1_bc .ne. 0) then
          call dfftw_destroy_plan(plan1_bc)
          plan1_bc = 0
      endif
      if (plan2_bc .ne. 0) then
          call dfftw_destroy_plan(plan2_bc)
          plan2_bc = 0
      endif
#else
      if (plan1_frc .ne. 0) then
          call sfftw_destroy_plan(plan1_frc)
          plan1_frc = 0
      endif
      if (plan1_bcr .ne. 0) then
          call sfftw_destroy_plan(plan1_bcr)
          plan1_bcr = 0
      endif
      if (plan1_fc .ne. 0) then
          call sfftw_destroy_plan(plan1_fc)
          plan1_fc = 0
      endif
      if (plan2_fc .ne. 0) then
          call sfftw_destroy_plan(plan2_fc)
          plan2_fc = 0
      endif
      if (plan1_bc .ne. 0) then
          call sfftw_destroy_plan(plan1_bc)
          plan1_bc = 0
      endif
      if (plan2_bc .ne. 0) then
          call sfftw_destroy_plan(plan2_bc)
          plan2_bc = 0
      endif
#endif

#elif defined ESSL
      if (allocated(caux1)) deallocate(caux1)
      if (allocated(caux2)) deallocate(caux2)
      if (allocated(raux1)) deallocate(raux1)
      if (allocated(raux2)) deallocate(raux2)
#endif

      if (allocated(buf1)) deallocate(buf1)
      if (allocated(buf2)) deallocate(buf2)
      if (allocated(buf) ) deallocate(buf)

! Clean-up what we allocated to allow re-setup
      if (allocated(proc_id2coords)) deallocate(proc_id2coords)
      if (allocated(proc_coords2id)) deallocate(proc_coords2id)

      if (allocated(iist)) deallocate(iist)
      if (allocated(iisz)) deallocate(iisz)
      if (allocated(iien)) deallocate(iien)
      if (allocated(jjst)) deallocate(jjst)
      if (allocated(jjsz)) deallocate(jjsz)
      if (allocated(jjen)) deallocate(jjen)
      if (allocated(jist)) deallocate(jist)
      if (allocated(jisz)) deallocate(jisz)
      if (allocated(jien)) deallocate(jien)
      if (allocated(kjst)) deallocate(kjst)
      if (allocated(kjsz)) deallocate(kjsz)
      if (allocated(kjen)) deallocate(kjen)

      if (allocated(IfSndStrt)) deallocate(IfSndStrt)
      if (allocated(IfSndCnts)) deallocate(IfSndCnts)
      if (allocated(IfRcvStrt)) deallocate(IfRcvStrt)
      if (allocated(IfRcvCnts)) deallocate(IfRcvCnts)
      if (allocated(KfSndStrt)) deallocate(KfSndStrt)
      if (allocated(KfSndCnts)) deallocate(KfSndCnts)
      if (allocated(KfRcvStrt)) deallocate(KfRcvStrt)
      if (allocated(KfRcvCnts)) deallocate(KfRcvCnts)
      if (allocated(JrSndStrt)) deallocate(JrSndStrt)
      if (allocated(JrSndCnts)) deallocate(JrSndCnts)
      if (allocated(JrRcvStrt)) deallocate(JrRcvStrt)
      if (allocated(JrRcvCnts)) deallocate(JrRcvCnts)
      if (allocated(KrSndStrt)) deallocate(KrSndStrt)
      if (allocated(KrSndCnts)) deallocate(KrSndCnts)
      if (allocated(KrRcvStrt)) deallocate(KrRcvStrt)
      if (allocated(KrRcvCnts)) deallocate(KrRcvCnts)

      mpi_set = .false.

      return
      end subroutine

      subroutine ar_copy(A,B,nar)

      integer(8) nar,i
      complex(mytype) A(nar,1,1),B(nar,1,1)

      do i=1,nar
         B(i,1,1)=A(i,1,1)
      enddo

      return
      end subroutine



      subroutine print_buf(A,lx,ly,lz)

      complex(mytype) A(lx,ly,lz)
      integer lx,ly,lz,i,j,k

      do k=1,lz
         do j=1,ly
            do i=1,lx
               if(abs(A(i,j,k)) .gt. 0.0000005) then
                  print *,taskid,': (',i,j,k,') =',A(i,j,k)
               endif
            enddo
         enddo
      enddo

      end subroutine

!=======================================================================
!                 ghost cell support
!               ----------------------
!
!  how to use:
! ==============
!     1) init ghost-support:
!       ==>> call p3dfft_init_ghosts(goverlap)
!       goverlap = overlap of ghost-cell-slab
!                  (limitations: only one goverlap supported)
!       (of cause we are still assuming periodic boundary conditions)
!
!     2) allocation of memory
!       each array which needs ghost-cell-support has to have
!       the min. size:  fsize(1)+2*goverlap *
!                       fsize(2)+2*goverlap *
!                       fsize(3)+2*goverlap
!
!     3) exchange ghost-cells among processors
!       to update the ghost-cells among the proceesors
!       ==>> call update_rghosts(array).
!       This will assume data represents physical space!
!
!     4) use/read ghost-cells
!       if you want to access/read ghost-cells
!       ==>> call ijk2i(i,j,k)
!       which returns the position in a 1D-array as an integer*8
!
!  array-memory:
! ==============
!     1) size if complex-type: fsize(1)+2*goverlap *
!                              fsize(2)+2*goverlap *
!                              fsize(3)+2*goverlap
!     2a) structure in memory:
!      |000000000000000000000000000|11111|22222|33333|44444|55555|66666|
!                 ^^^^              ^     ^     ^     ^     ^     ^
!          input/output of FFT      A     B     C     D     E     F
!
!      ghost-cell-slab at each pencil-side:
!      A:  -X,Y*Z       A(istart(1)-i, j,k)             => gmem_start(1)
!      B:  +X,Y*Z       B(i,j,k)                                => gmem_start(2)
!      C:  -Y,X*Z       C(i,j,k)                                => gmem_start(3)
!      D:  +Y,X*Z       D(i,j,k)                                => gmem_start(4)
!      E:  -Z,X*Y       E(i,j,k)                                => gmem_start(5)
!      F:  +Z,X*Y       F(i,j,k)                                => gmem_start(6)
!
!     2b) a more closer look at "A":
!      ....|aaaaaabbbbbbccccccdddddd|....
!           ^     ^     ^     ^
!           a1    b1    c1    d1
!
!      each ghost-cell-slab is build up by no. goverlap planes:
!      a1:  Y*Z-plane with lowest X-coord
!       |
!      d1:  Y*Z-plane with highest X-coord
!
!  important variables:
! =====================
!     1) goverlap
!       the data overlap between two processors - equal to the thickness
!       of a ghost-slab
!       This version supports only one fixed overlap value, because it simplifies
!       the coding. If different overlaps for different arrays are needed the
!       following varables must be dependent on the overlap-value and therefor
!       should have one additional dimension.
!       A good idea would be to make p3dfft_init_ghosts(..) return an id,which
!       has to get set when calling update_ghosts(data,id).
!
!     2) gmem_start(1..6) (-X,-Y,-Z,+X,+Y,+Z)
!       start of data in memory for ghost-cells of certain side
!
!     3) gmess_size(1..3) (X,Y,Z)
!       size of mpi message needed to send/recieve ghost-slab each dimension(1..3)
!
!     4a) gslab_start(1..3, 1..6)
!       relative start coord. of ghost-slab in pencil ijk(1..3), [neg./pos side(1..2) for each dimension(1..3)](1..6)
!           -X   -Y   -Z   +X   +Y   +Z
!         i
!         j
!         k
!
!     4b) gslab_end(1..3, 1..6)
!       relative end coord. of ghost-slab in pencil ijk(1..3), [neg./pos side(1..2) for each dimension(1..3)](1..6)
!           -X   -Y   -Z   +X   +Y   +Z
!         i
!         j
!         k
!
!     4c) gslab_size(1..3, 1..3)
!       size of ghost-slab for each dimension(1..3), ijk(1..3)
!
!     6) gneighb_(numtasks*(3*2))
!       extended processor grid with all processor neighbours in -/+XYZ(1..6)
!       including periodic boundary neighbours
!        example for 2 processors with gneighb_r(:) = 0,0,0,0,1,1,1,1,1,1,0,0
!         pencil-side: -X  +X  -Y  +Y  -Z  +Z
!         taskid 0:     0   0   0   0   1   1
!         taskid 1:     1   1   1   1   0   0
!=======================================================================

      subroutine p3dfft_init_ghosts(goverlap_in)
        implicit none

!       function args
        integer, intent(in) :: goverlap_in

!       other vars
        integer, dimension(:,:), allocatable :: coords2id
        integer :: pneighb_r(6), pneighb_c(6)
        integer :: pCoo(2)
        integer :: rstart(3), rend(3), cstart(3), cend(3)
        integer :: i,d,mp,mp1,msize,ierr,dmp,m,p

        if(.not. mpi_set) then
                print *,'P3DFFT error: call setup before other routines'
                ierr = 1
                call abort
        endif

        goverlap = goverlap_in

        call p3dfft_get_dims(gproc_rstart, gproc_rend, gproc_rsize, 1)
        call p3dfft_get_dims(gproc_cstart, gproc_cend, gproc_csize, 2)

        if(     goverlap .gt. gproc_csize(1) &
           .or. goverlap .gt. gproc_csize(2) &
           .or. goverlap .gt. gproc_csize(3) ) then
                        print *,'P3DFFT error: in one direction goverlap is greater then complex pencil size'
                        ierr = 1
                        call abort
        endif

!       set mpi message sizes for exchange of ghost-cells between processors
        gmess_rsize(1) = goverlap       *gproc_rsize(2) *gproc_rsize(3)
        gmess_rsize(2) = gproc_rsize(1) *goverlap       *gproc_rsize(3)
        gmess_rsize(3) = gproc_rsize(1) *gproc_rsize(2) *goverlap

        gmess_csize(1) = goverlap       *gproc_csize(2) *gproc_csize(3) *2
        gmess_csize(2) = gproc_csize(1) *goverlap       *gproc_csize(3) *2
        gmess_csize(3) = gproc_csize(1) *gproc_csize(2) *goverlap       *2

!       precalc ghost-slab start-position in memory of in/out array (-X,-Y
        gmem_rstart(1) = int(gproc_csize(1),8)  *int(gproc_csize(2),8) *int(gproc_csize(3),8) *2 +1
        gmem_rstart(2) = gmem_rstart(1) +int(goverlap      *gproc_rsize(2) *gproc_rsize(3) ,8)
        gmem_rstart(3) = gmem_rstart(2) +int(gproc_rsize(1) *goverlap      *gproc_rsize(3) ,8)
        gmem_rstart(4) = gmem_rstart(3) +int(gproc_rsize(1) *gproc_rsize(2) *goverlap      ,8)
        gmem_rstart(5) = gmem_rstart(4) +int(goverlap       *gproc_rsize(2) *gproc_rsize(3),8)
        gmem_rstart(6) = gmem_rstart(5) +int(gproc_rsize(1) *goverlap       *gproc_rsize(3),8)

        gmem_cstart(1) = int(gproc_csize(1),8)  *int(gproc_csize(2),8) *int(gproc_csize(3),8) *2 +1
        gmem_cstart(2) = gmem_cstart(1) +int(goverlap      *gproc_csize(2) *gproc_csize(3) ,8)
        gmem_cstart(3) = gmem_cstart(2) +int(gproc_csize(1) *goverlap      *gproc_csize(3) ,8)
        gmem_cstart(4) = gmem_cstart(3) +int(gproc_csize(1) *gproc_csize(2) *goverlap      ,8)
        gmem_cstart(5) = gmem_cstart(4) +int(goverlap       *gproc_csize(2) *gproc_csize(3),8)
        gmem_cstart(6) = gmem_cstart(5) +int(gproc_csize(1) *goverlap       *gproc_csize(3),8)

!       write(*,*) 'gmem_start: ', gmem_rstart(:)

!       allocate send/recieve buffer
        msize = gmess_rsize(1)
        if(gmess_rsize(2) .gt. msize) msize = gmess_rsize(2)
        if(gmess_rsize(3) .gt. msize) msize = gmess_rsize(3)
        if(gmess_csize(1) .gt. msize) msize = gmess_csize(1)
        if(gmess_csize(2) .gt. msize) msize = gmess_csize(2)
        if(gmess_csize(3) .gt. msize) msize = gmess_csize(3)
        allocate(gbuf_snd(msize),  stat=ierr)
        allocate(gbuf_recv(msize), stat=ierr)

!       precalc ghost-slab size,start,end (checked)
        gslab_rsize(:,1) = gproc_rsize(:)
        gslab_rsize(:,2) = gproc_rsize(:)
        gslab_rsize(:,3) = gproc_rsize(:)
        gslab_csize(:,1) = gproc_csize(:)
        gslab_csize(:,2) = gproc_csize(:)
        gslab_csize(:,3) = gproc_csize(:)

        gslab_rsize(1,1) = goverlap
        gslab_rsize(2,2) = goverlap
        gslab_rsize(3,3) = goverlap
        gslab_csize(1,1) = goverlap
        gslab_csize(2,2) = goverlap
        gslab_csize(3,3) = goverlap
        do d=1, 3       ! loop over X,Y,Z sides
          do mp=1,2             ! loop over direction (-/+)
                dmp = (mp-1)*3 +d
                if(mp .eq. 1) then; m=1; p=0; else; m=0; p=1; endif

                gslab_rstart(1,dmp) = m*1 +p*(gproc_rsize(1)-gslab_rsize(1,d)+1)
                gslab_rstart(2,dmp) = m*1 +p*(gproc_rsize(2)-gslab_rsize(2,d)+1)
                gslab_rstart(3,dmp) = m*1 +p*(gproc_rsize(3)-gslab_rsize(3,d)+1)

                gslab_rend(1,dmp) = gslab_rstart(1,dmp) +gslab_rsize(1,d) -1
                gslab_rend(2,dmp) = gslab_rstart(2,dmp) +gslab_rsize(2,d) -1
                gslab_rend(3,dmp) = gslab_rstart(3,dmp) +gslab_rsize(3,d) -1

                gslab_cstart(1,dmp) = m*1 +p*(gproc_csize(1)-gslab_csize(1,d)+1)
                gslab_cstart(2,dmp) = m*1 +p*(gproc_csize(2)-gslab_csize(2,d)+1)
                gslab_cstart(3,dmp) = m*1 +p*(gproc_csize(3)-gslab_csize(3,d)+1)

                gslab_cend(1,dmp) = gslab_rstart(1,dmp) +gslab_rsize(1,d) -1
                gslab_cend(2,dmp) = gslab_rstart(2,dmp) +gslab_rsize(2,d) -1
                gslab_cend(3,dmp) = gslab_rstart(3,dmp) +gslab_rsize(3,d) -1
          enddo
        enddo

!       print *, 'gproc_rstart ',taskid, gproc_rstart(:)
!       call MPI_Barrier(mpi_comm_cart,ierr)
!       print *, 'gproc_rend ',taskid, gproc_rend(:)
!       call MPI_Barrier(mpi_comm_world,ierr)
!       do d=1, 3
!               write(*,*) 'gslab_rsize ',taskid,d, gslab_rsize(d,:)
!       enddo
!       call MPI_Barrier(mpi_comm_cart,ierr)
!       do d=1, 3
!               write(*,*) 'gslab_rstart',taskid,d, gslab_rstart(d,:)
!       enddo
!       call MPI_Barrier(mpi_comm_cart,ierr)
!       do d=1, 3
!               write(*,*) 'gslab_rend  ',taskid,d, gslab_rend(d,:)
!       enddo
!       call MPI_Barrier(mpi_comm_cart,ierr)


!     ! build help-array with coordinates to proc-ids
        allocate(coords2id(-1:iproc,-1:jproc))
        coords2id(:,:) = -1
        do i=0,numtasks-1
                coords2id( proc_id2coords(2*i), proc_id2coords(2*i+1)) = i
        enddo
!     ! extend array with one line more on each side (periodicy)
        coords2id(-1,   0:jproc-1) = coords2id(iproc-1,0:jproc-1)
        coords2id(iproc,0:jproc-1) = coords2id(0,      0:jproc-1)
        coords2id(0:iproc-1,   -1) = coords2id(0:iproc-1,jproc-1)
        coords2id(0:iproc-1,jproc) = coords2id(0:iproc-1,      0)

        pCoo(1) = proc_id2coords(taskid*2)
        pCoo(2) = proc_id2coords(taskid*2+1)
        do d=0, 2       ! loop over X,Y,Z sides
          do mp=1,2             ! loop over direction (-/+)
                mp1 = mp*2-3    ! 1=>-1 , 2=>+1

!               ! build pneighb_r
                if(     d .eq. 0 .or.                         & ! no decomp. in X
                        d .eq. 1 .and. iproc .eq. 1 .or.      & ! no decomp. in Y
                        d .eq. 2 .and. jproc .eq. 1) then       ! no decomp. in Z
                                pneighb_r(d*2+mp) = taskid
                else if(d .eq. 1) then                                  ! decomp. in Y
                        pneighb_r(d*2+mp) = coords2id(pCoo(1)+mp1, pCoo(2))
                else if(d .eq. 2) then                                  ! decomp. in Z
                        pneighb_r(d*2+mp) = coords2id(pCoo(1), pCoo(2)+mp1)
                endif

!               ! build pneighb_c
                if(     d .eq. 0 .and. iproc .eq. 1 .or.      & ! no decomp. in X
                        d .eq. 1 .and. jproc .eq. 1 .or.      & ! no decomp. in Y
                        d .eq. 2) then                          ! no decomp. in Z
                                pneighb_c(d*2+mp) = taskid
                else if(d .eq. 0) then                                  ! decomp. in X
                        pneighb_c(d*2+mp) = coords2id(pCoo(1)+mp1, pCoo(2))
                else if(d .eq. 1) then                                  ! decomp. in Y
                        pneighb_c(d*2+mp) = coords2id(pCoo(1), pCoo(2)+mp1)
                endif

          enddo
        enddo

!     ! distribute pneigh_r
        allocate(gneighb_r(numtasks*6))
        call MPI_Allgather(     pneighb_r,      6, MPI_INTEGER, &
                                                gneighb_r,      6, MPI_INTEGER, &
                                                mpi_comm_cart,ierr)

!     ! distribute pneigh_c
        allocate(gneighb_c(numtasks*6))
        call MPI_Allgather(     pneighb_c,      6, MPI_INTEGER, &
                                                gneighb_c,      6, MPI_INTEGER, &
                                                mpi_comm_cart,ierr)

!       if(taskid .eq. 0) write(*,*) 'gneighb: ',gneighb_r(:)

        ghosts_set = .true.

      end subroutine p3dfft_init_ghosts

!     ------------------------------------------------------------------
!     !  p3dfft_update_ghosts(XYZ, goverlap)
!     !     # XYZ - array of real made out of two parts
!     !              1. part = input data = FFT-pencil
!     !                        [(fsize(1)*2) *(fsize(2)*2) *(fsize(3)*2)]
!     !              2. part = output data = ghost-cells
!     !                         [(fsize(1)*2) *(fsize(2)*2) *goverlap]*2        2* XY-plane
!     !                        +[(fsize(1)*2) *(fsize(3)*2) *goverlap]*2        2* XZ-plane
!     !                        +[(fsize(2)*2) *(fsize(3)*2) *goverlap]*2        2* YZ-plane
!     !     # gneighb    - neighbours of all pencils in all directions

      subroutine p3dfft_update_ghosts(XYZ, gneighb, gproc_size, &
                                     gmem_start, gslab_start, gslab_end)
        implicit none

!     ! function args
        real(mytype) ,TARGET :: XYZ(*)
        integer, intent(in) :: gneighb(*)
        integer, intent(in) :: gproc_size(*)
        integer(kind=8), intent(in) :: gmem_start(*)
        integer, intent(in) :: gslab_start(3,6)
        integer, intent(in) :: gslab_end(3,6)

!     ! usual stuff
        integer,save :: d,mp, dmp, pmd
        integer,save :: send_neighb, recv_neighb
        integer(kind=8),save :: i,j,k,jk
        integer,save :: g, buf_size
        integer(kind=8),save :: jpos,kpos
        integer, save :: request, status(mpi_status_size), ierr

!     ! exchange ghost-cells sides of pencil/slide
        do d=1, 3       ! X,Y,Z sides
                do mp=1,2
                        if(mp .eq. 1) then      ! from - to +
                                send_neighb = (taskid*6) +(d-1)*2 +2
                                recv_neighb = (taskid*6) +(d-1)*2 +1
                                dmp = d
                                pmd = 3+d
                        else                            ! from + to -
                                send_neighb = (taskid*6) +(d-1)*2 +1
                                recv_neighb = (taskid*6) +(d-1)*2 +2
                                dmp = 3+d
                                pmd = d
                        endif

!                       dimension d is contained entirely within processors memory
                        if(gneighb(send_neighb) .eq. taskid) then
!                          .and. gneighb(recv_neighb) .eq. taskid) should always be the case if gneighb(send_neighb)=taskid

!                               just copy data in local memory
                                g=0
                                do k=gslab_start(3,dmp), gslab_end(3,dmp)
                                  kpos = (k-1)*gproc_size(1)*gproc_size(2)
                                  do j=gslab_start(2,dmp), gslab_end(2,dmp)
                                        jpos = kpos +(j-1)*gproc_size(1)
                                        do i=gslab_start(1,dmp), gslab_end(1,dmp)
                                          XYZ(gmem_start(pmd) +g) = XYZ(jpos +i)
                                          g = g+1
                                        enddo
                                  enddo
                                enddo

!                       dimension d is block-distributed among processors
                        else

!                               pack send-buffer
                                g=0
                                do k=gslab_start(3,dmp), gslab_end(3,dmp)
                                  kpos = (k-1)*gproc_size(1)*gproc_size(2)
                                  do j=gslab_start(2,dmp), gslab_end(2,dmp)
                                        jpos = kpos +(j-1)*gproc_size(1)
                                        do i=gslab_start(1,dmp), gslab_end(1,dmp)
                                          g = g+1
                                          gbuf_snd(g) = XYZ(jpos +i)
                                        enddo
                                  enddo
                                enddo

!                               start a receive, send ghost-cells then wait
                                buf_size = g*mytype
                                call MPI_IRecv( gbuf_recv, buf_size, MPI_BYTE, &
                                                gneighb(send_neighb), d, mpi_comm_cart, &
                                                request, ierr)

                                call MPI_Send(  gbuf_snd, buf_size, MPI_BYTE, &
                                                gneighb(recv_neighb), d, mpi_comm_cart, ierr)

                                call MPI_Wait( request, status, ierr )

!                               unpack recieve-buffer
                                do i=1, g
                                        XYZ(gmem_start(pmd)+i-1) = gbuf_recv(i)
                                enddo

                        endif

                enddo
        enddo

      end subroutine p3dfft_update_ghosts

!     ------------------------------------------------------------------
!
!       ijk2i() - calc 1d-array-position
!
      function gr_ijk2i(i,j,k)
        implicit none

!       function args
        integer, intent(in) :: i,j,k
        integer(kind=8) :: gr_ijk2i

        if(.not. ghosts_set) then
                print *,'P3DFFT error: call p3dfft_init_ghosts before'
                return
        endif

!       ghost-cell in i-dim
        if(i .lt. gproc_rstart(1)) then
          gr_ijk2i = gmem_rstart(1) &
                +int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8) &
                +int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8) &
                +                                           int( i-gproc_rstart(1) +gslab_rsize(1,1),8)
                return
        else if(i .gt. gproc_rend(1)  ) then
          gr_ijk2i = gmem_rstart(4) &
                +int( gslab_rsize(1,1)*gslab_rsize(1,2),8) *int( k-gproc_rstart(3),8) &
                +int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8) &
                +                                           int( i-gproc_rend(1)-1,8)
                return

!       ghost-cell in j-dim
        else if(j .lt. gproc_rstart(2)) then
          gr_ijk2i = gmem_rstart(2) &
                +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8) &
                +int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,2),8) &
                +                                           int( i-gproc_rstart(1),8)
                return
        else if(j .gt. gproc_rend(2)  ) then
          gr_ijk2i = gmem_rstart(5) &
                +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8) &
                +int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rend(2)-1,8) &
                +                                           int( i-gproc_rstart(1),8)
                return

!       ghost-cell in k-dim
        else if(k .lt. gproc_rstart(3)) then
          gr_ijk2i = gmem_rstart(3) &
                +int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rstart(3) +gslab_rsize(3,3),8) &
                +int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8) &
                +                                           int( i-gproc_rstart(1),8)
        else if(k .gt. gproc_rend(3)  ) then
          gr_ijk2i = gmem_rstart(6) &
                +int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rend(3)-1,8) &
                +int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8) &
                +                                           int( i-gproc_rstart(1),8)
                return

!       no ghost-cell
        else
          gr_ijk2i = 1 &
                +int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8) &
                +int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8) &
                +                                      int(i-gproc_rstart(1),8)
        endif

      end function gr_ijk2i

!     ------------------------------------------------------------------
      subroutine update_rghosts(XgYZ)
        real(mytype) ,TARGET :: XgYZ(1,1,*)

        if(.not. ghosts_set) then
                print *,'P3DFFT error: call p3dfft_init_ghosts before'
                return
        endif
        call p3dfft_update_ghosts(XgYZ, gneighb_r, gproc_rsize, &
                                   gmem_rstart, gslab_rstart, gslab_rend)

      end subroutine update_rghosts

!     ------------------------------------------------------------------
      subroutine update_cghosts(XYZg)
        real(mytype) ,TARGET :: XYZg(1,1,*)

        print *,'P3DFFT error: support for complex ghost-cells not implemented'
        print *,'              in p3dfft_update_ghosts() and ijk2i() yet,'
        print *,'              perhaps you like to do the job (?)'
        return

        if(.not. ghosts_set) then
                print *,'P3DFFT error: call p3dfft_init_ghosts before'
                return
        endif
        call p3dfft_update_ghosts(XYZg, gneighb_c, gproc_csize, &
                                   gmem_cstart, gslab_cstart, gslab_cend)

      end subroutine update_cghosts

      end module


