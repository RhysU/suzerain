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

module largo

! Use ISO_C_BINDING to expose C-friendly API
  use, intrinsic :: iso_c_binding, only: c_associated,      &
                                         WP => c_double,    &
                                         c_int,             &
                                         c_char,            &
                                         c_loc,             &
                                         c_f_pointer,       &
                                         c_null_char,       &
                                         c_null_ptr,        &
                                         largo_ptr =>c_ptr, &
                                         largo_workspace_ptr =>c_ptr, &
                                         string_ptr =>c_ptr


! Use each largo module exporting all publicly visible symbols from each one
  use largo_workspace
  use largo_bl_spatial
  use largo_bl_spatial_first_order
  use largo_bl_spatial_mixed
  use largo_bl_spatiotemporal
  use largo_bl_temporal_chemistry
  use largo_bl_temporal
  use largo_bl_temporal_tconsistent

  implicit none ! is in effect throughout entire module

  public :: largo_allocate
  public :: largo_allocate_c

contains

  ! Generic interface routines

  ! Generic interface allocate
  subroutine largo_allocate_c(lcp, neq, ns, cmodel) bind(C, name="largo_allocate")

    ! largo workspace C pointer
    type(largo_ptr), intent(out)                :: lcp
    integer(c_int), intent(in), value           :: neq
    integer(c_int), intent(in), value           :: ns     ! number of species
    type(string_ptr), intent(in), value         :: cmodel
    character(len=255)                          :: fmodel
    logical                                     :: copy_success

    ! Copy C-style, NULL terminated string to Fortran-compatible storage
    copy_success = largo_c_f_stringcopy (cmodel, fmodel)

    ! Invoke Fortran functionality using Fortran-ready model string
    call largo_allocate(lcp, neq, ns, fmodel)

  end subroutine largo_allocate_c


  ! Generic interface allocate
  subroutine largo_allocate(lcp, neq, ns, fmodel) !bind(C)

    ! largo workspace C pointer
    type(largo_ptr), intent(out)         :: lcp
    integer(c_int), intent(in)           :: neq
    integer(c_int), intent(in)           :: ns    ! number of species
    character(*), intent(in)             :: fmodel
    type(largo_type), pointer            :: lauxp

    ! Allocate derived type variable
    allocate (lauxp)

    ! Initialize number of equations
    lauxp%neq = neq

    ! Initialize number of species
    lauxp%ns = ns

    ! Initialize number of variables
    lauxp%nvar = neq + 1

    ! Initialize according to model index
    ! FIXME: enumerate models
    select case (trim(fmodel))
    case ("bl_temporal")
      call largo_BL_temporal_allocate (lauxp%cp, lauxp%neq, lauxp%ns)
      lauxp%largo_init            => largo_BL_temporal_init
      lauxp%largo_finalize        => largo_BL_temporal_deallocate
      lauxp%largo_prestep_mean    => largo_BL_temporal_preStep_sEtaMean
      lauxp%largo_prestep_innerxz => largo_BL_temporal_preStep_sEta_innerxz
      lauxp%largo_prestep_innery  => largo_BL_temporal_preStep_sEta_innery
      lauxp%largo_prestep         => largo_BL_temporal_preStep_sEta

      lauxp%largo_continuity_mean => largo_BL_temporal_continuity_sEtaMean
      lauxp%largo_xmomentum_mean  => largo_BL_temporal_xMomentum_sEtaMean
      lauxp%largo_ymomentum_mean  => largo_BL_temporal_yMomentum_sEtaMean
      lauxp%largo_zmomentum_mean  => largo_BL_temporal_zMomentum_sEtaMean
      lauxp%largo_energy_mean     => largo_BL_temporal_energy_sEtaMean
      lauxp%largo_species_mean    => largo_BL_temporal_species_sEtaMean
      lauxp%largo_ispecies_mean   => largo_BL_temporal_ispecies_sEtaMean
      lauxp%largo_continuity_rms  => largo_BL_temporal_continuity_sEtaRms
      lauxp%largo_xmomentum_rms   => largo_BL_temporal_xMomentum_sEtaRms
      lauxp%largo_ymomentum_rms   => largo_BL_temporal_yMomentum_sEtaRms
      lauxp%largo_zmomentum_rms   => largo_BL_temporal_zMomentum_sEtaRms
      lauxp%largo_energy_rms      => largo_BL_temporal_energy_sEtaRms
      lauxp%largo_species_rms     => largo_BL_temporal_species_sEtaRms
      lauxp%largo_ispecies_rms    => largo_BL_temporal_ispecies_sEtaRms

      lauxp%largo_continuity      => largo_BL_temporal_continuity_sEta
      lauxp%largo_xmomentum       => largo_BL_temporal_xMomentum_sEta
      lauxp%largo_ymomentum       => largo_BL_temporal_yMomentum_sEta
      lauxp%largo_zmomentum       => largo_BL_temporal_zMomentum_sEta
      lauxp%largo_energy          => largo_BL_temporal_energy_sEta
      lauxp%largo_species         => largo_BL_temporal_species_sEta

      lauxp%largo_all_sources_mean => largo_BL_temporal_sEtaMean
      lauxp%largo_all_sources      => largo_BL_temporal_sEta

      lauxp%largo_prestep_baseflow => largo_BL_temporal_preStep_baseflow

    case ("bl_temporal_tensor-consistent")
      call largo_BL_temporal_tconsistent_allocate (lauxp%cp, lauxp%neq, lauxp%ns)
      lauxp%largo_init            => largo_BL_temporal_tconsistent_init
      lauxp%largo_finalize        => largo_BL_temporal_tconsistent_deallocate
      lauxp%largo_prestep_mean    => largo_BL_temporal_tconsistent_preStep_sEtaMean
      lauxp%largo_prestep_innerxz => largo_BL_temporal_tconsistent_preStep_sEta_innerxz
      lauxp%largo_prestep_innery  => largo_BL_temporal_tconsistent_preStep_sEta_innery
      lauxp%largo_prestep         => largo_BL_temporal_tconsistent_preStep_sEta

      lauxp%largo_continuity_mean => largo_BL_temporal_tconsistent_continuity_sEtaMean
      lauxp%largo_xmomentum_mean  => largo_BL_temporal_tconsistent_xMomentum_sEtaMean
      lauxp%largo_ymomentum_mean  => largo_BL_temporal_tconsistent_yMomentum_sEtaMean
      lauxp%largo_zmomentum_mean  => largo_BL_temporal_tconsistent_zMomentum_sEtaMean
      lauxp%largo_energy_mean     => largo_BL_temporal_tconsistent_energy_sEtaMean
      lauxp%largo_species_mean    => largo_BL_temporal_tconsistent_species_sEtaMean

      lauxp%largo_continuity      => largo_BL_temporal_tconsistent_continuity_sEta_
      lauxp%largo_xmomentum       => largo_BL_temporal_tconsistent_xMomentum_sEta_
      lauxp%largo_ymomentum       => largo_BL_temporal_tconsistent_yMomentum_sEta_
      lauxp%largo_zmomentum       => largo_BL_temporal_tconsistent_zMomentum_sEta_
      lauxp%largo_energy          => largo_BL_temporal_tconsistent_energy_sEta_
      lauxp%largo_species         => largo_BL_temporal_tconsistent_species_sEta_

      lauxp%largo_all_sources     => largo_BL_temporal_tconsistent_sEta

      lauxp%largo_prestep_baseflow => largo_BL_temporal_tconsistent_preStep_baseflow

    case ("bl_spatiotemporal")
      call largo_BL_spatiotemporal_allocate (lauxp%cp, lauxp%neq, lauxp%ns)
      lauxp%largo_init            => largo_BL_spatiotemporal_init
      lauxp%largo_finalize        => largo_BL_spatiotemporal_deallocate
      lauxp%largo_prestep_mean    => largo_BL_spatiotemporal_preStep_sEtaMean
      lauxp%largo_prestep_innerxz => largo_BL_spatiotemporal_preStep_sEta_innerxz
      lauxp%largo_prestep_innery  => largo_BL_spatiotemporal_preStep_sEta_innery
      lauxp%largo_prestep         => largo_BL_spatiotemporal_preStep_sEta

      lauxp%largo_continuity_mean => largo_BL_spatiotemporal_continuity_sEtaMean
      lauxp%largo_xmomentum_mean  => largo_BL_spatiotemporal_xMomentum_sEtaMean
      lauxp%largo_ymomentum_mean  => largo_BL_spatiotemporal_yMomentum_sEtaMean
      lauxp%largo_zmomentum_mean  => largo_BL_spatiotemporal_zMomentum_sEtaMean
      lauxp%largo_energy_mean     => largo_BL_spatiotemporal_energy_sEtaMean
      lauxp%largo_species_mean    => largo_BL_spatiotemporal_species_sEtaMean
      lauxp%largo_ispecies_mean   => largo_BL_spatiotemporal_ispecies_sEtaMean
      lauxp%largo_continuity_rms  => largo_BL_spatiotemporal_continuity_sEtaRms
      lauxp%largo_xmomentum_rms   => largo_BL_spatiotemporal_xMomentum_sEtaRms
      lauxp%largo_ymomentum_rms   => largo_BL_spatiotemporal_yMomentum_sEtaRms
      lauxp%largo_zmomentum_rms   => largo_BL_spatiotemporal_zMomentum_sEtaRms
      lauxp%largo_energy_rms      => largo_BL_spatiotemporal_energy_sEtaRms
      lauxp%largo_species_rms     => largo_BL_spatiotemporal_species_sEtaRms
      lauxp%largo_ispecies_rms    => largo_BL_spatiotemporal_ispecies_sEtaRms

      lauxp%largo_continuity      => largo_BL_spatiotemporal_continuity_sEta
      lauxp%largo_xmomentum       => largo_BL_spatiotemporal_xMomentum_sEta
      lauxp%largo_ymomentum       => largo_BL_spatiotemporal_yMomentum_sEta
      lauxp%largo_zmomentum       => largo_BL_spatiotemporal_zMomentum_sEta
      lauxp%largo_energy          => largo_BL_spatiotemporal_energy_sEta
      lauxp%largo_species         => largo_BL_spatiotemporal_species_sEta

      lauxp%largo_all_sources_mean => largo_BL_spatiotemporal_sEtaMean
      lauxp%largo_all_sources      => largo_BL_spatiotemporal_sEta

      lauxp%largo_init_wall_baseflow => largo_BL_spatiotemporal_init_wall_baseflow
      lauxp%largo_prestep_baseflow   => largo_BL_spatiotemporal_preStep_baseflow

    case default
      ! FIXME: throw an error, model not declared
    end select

    ! Get C pointer from Fortran pointer
    lcp = c_loc(lauxp)

  end subroutine largo_allocate


  ! Generic interface allocate
  subroutine largo_deallocate(lcp) bind(C)

    ! largo C pointer
    type(largo_ptr), intent(out) :: lcp
    type(largo_type), pointer    :: lauxp

    call c_f_pointer(lcp, lauxp)

    ! Deallocate the model in a model-specific fashion
    if (associated(lauxp%largo_finalize) .and. c_associated(lauxp%cp)) then
      call lauxp%largo_finalize(lauxp%cp)
    end if

    ! Nullify procedure pointers as a defensive precaution
    if (associated(lauxp%largo_init            )) nullify(lauxp%largo_init            )
    if (associated(lauxp%largo_finalize        )) nullify(lauxp%largo_finalize        )
    if (associated(lauxp%largo_prestep_mean    )) nullify(lauxp%largo_prestep_mean    )
    if (associated(lauxp%largo_prestep_innerxz )) nullify(lauxp%largo_prestep_innerxz )
    if (associated(lauxp%largo_prestep_innery  )) nullify(lauxp%largo_prestep_innery  )
    if (associated(lauxp%largo_prestep         )) nullify(lauxp%largo_prestep         )
    if (associated(lauxp%largo_continuity_mean )) nullify(lauxp%largo_continuity_mean )
    if (associated(lauxp%largo_xmomentum_mean  )) nullify(lauxp%largo_xmomentum_mean  )
    if (associated(lauxp%largo_ymomentum_mean  )) nullify(lauxp%largo_ymomentum_mean  )
    if (associated(lauxp%largo_zmomentum_mean  )) nullify(lauxp%largo_zmomentum_mean  )
    if (associated(lauxp%largo_energy_mean     )) nullify(lauxp%largo_energy_mean     )
    if (associated(lauxp%largo_species_mean    )) nullify(lauxp%largo_species_mean    )
    if (associated(lauxp%largo_ispecies_mean   )) nullify(lauxp%largo_ispecies_mean   )
    if (associated(lauxp%largo_continuity_rms  )) nullify(lauxp%largo_continuity_rms  )
    if (associated(lauxp%largo_xmomentum_rms   )) nullify(lauxp%largo_xmomentum_rms   )
    if (associated(lauxp%largo_ymomentum_rms   )) nullify(lauxp%largo_ymomentum_rms   )
    if (associated(lauxp%largo_zmomentum_rms   )) nullify(lauxp%largo_zmomentum_rms   )
    if (associated(lauxp%largo_energy_rms      )) nullify(lauxp%largo_energy_rms      )
    if (associated(lauxp%largo_species_rms     )) nullify(lauxp%largo_species_rms     )
    if (associated(lauxp%largo_ispecies_rms    )) nullify(lauxp%largo_ispecies_rms    )

    if (associated(lauxp%largo_continuity      )) nullify(lauxp%largo_continuity      )
    if (associated(lauxp%largo_xmomentum       )) nullify(lauxp%largo_xmomentum       )
    if (associated(lauxp%largo_ymomentum       )) nullify(lauxp%largo_ymomentum       )
    if (associated(lauxp%largo_zmomentum       )) nullify(lauxp%largo_zmomentum       )
    if (associated(lauxp%largo_energy          )) nullify(lauxp%largo_energy          )
    if (associated(lauxp%largo_species         )) nullify(lauxp%largo_species         )

    if (associated(lauxp%largo_all_sources_mean)) nullify(lauxp%largo_all_sources_mean)
    if (associated(lauxp%largo_all_sources     )) nullify(lauxp%largo_all_sources     )

    if (associated(lauxp%largo_prestep_baseflow)) nullify(lauxp%largo_prestep_baseflow)
    if (associated(lauxp%largo_init_wall_baseflow)) &
                                                  nullify(lauxp%largo_init_wall_baseflow)

    ! Deallocate derived type variable
    deallocate(lauxp)

    ! Nullify C pointer
    lcp = c_null_ptr

  end subroutine largo_deallocate


  ! Generic interface prestep subroutine
  subroutine largo_init(lcp, grDelta, grDA, grDArms) bind(C)

    real(WP), intent(in), value         :: grDelta
    real(WP), dimension(*), intent(in)  :: grDA
    real(WP), dimension(*), intent(in)  :: grDArms
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_init(lauxp%cp, grDelta,                 &
                                    grDA     (1:lauxp%nvar), &
                                    grDArms  (1:lauxp%nvar))

  end subroutine largo_init


  ! Generic interface prestep subroutine
  subroutine largo_preStep_sEtaMean(lcp, y, mean, ddy_mean) bind(C)

    real(WP), intent(in), value         :: y
    real(WP), dimension(*), intent(in)  :: mean
    real(WP), dimension(*), intent(in)  :: ddy_mean
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_prestep_mean(lauxp%cp, y, &
                             mean         (1:lauxp%nvar),  &
                             ddy_mean     (1:lauxp%nvar))

  end subroutine largo_preStep_sEtaMean


  ! Generic interface prestep subroutine
  subroutine largo_preStep_sEta_innery(lcp, y,                          &
                                           mean,     rms,     mean_rqq, &
                                       ddy_mean, ddy_rms, ddy_mean_rqq) bind(C)

    real(WP), intent(in), value         :: y
    real(WP), dimension(*), intent(in)  :: mean
    real(WP), dimension(*), intent(in)  :: rms
    real(WP), dimension(*), intent(in)  :: mean_rqq
    real(WP), dimension(*), intent(in)  :: ddy_mean
    real(WP), dimension(*), intent(in)  :: ddy_rms
    real(WP), dimension(*), intent(in)  :: ddy_mean_rqq
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_prestep_innery(lauxp%cp, y, &
                             mean         (1:lauxp%nvar), &
                             rms          (1:lauxp%nvar), &
                             mean_rqq     (1:lauxp%nvar), &
                             ddy_mean     (1:lauxp%nvar), &
                             ddy_rms      (1:lauxp%nvar), &
                             ddy_mean_rqq (1:lauxp%nvar))

  end subroutine largo_preStep_sEta_innery


  ! Generic interface prestep subroutine
  subroutine largo_preStep_sEta_innerxz(lcp, qflow) bind(C)

    real(WP), dimension(*), intent(in)  :: qflow
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_prestep_innerxz(lauxp%cp,   &
                             qflow    (1:lauxp%nvar))

  end subroutine largo_preStep_sEta_innerxz


  ! Generic interface prestep subroutine
  subroutine largo_preStep_sEta(lcp, y, qflow,                   &
                                    mean,     rms,     mean_rqq, &
                                ddy_mean, ddy_rms, ddy_mean_rqq) bind(C)

    real(WP), intent(in), value         :: y
    real(WP), dimension(*), intent(in)  :: qflow
    real(WP), dimension(*), intent(in)  :: mean
    real(WP), dimension(*), intent(in)  :: rms
    real(WP), dimension(*), intent(in)  :: mean_rqq
    real(WP), dimension(*), intent(in)  :: ddy_mean
    real(WP), dimension(*), intent(in)  :: ddy_rms
    real(WP), dimension(*), intent(in)  :: ddy_mean_rqq
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_prestep(lauxp%cp, y,       &
                             qflow        (1:lauxp%nvar), &
                             mean         (1:lauxp%nvar), &
                             rms          (1:lauxp%nvar), &
                             mean_rqq     (1:lauxp%nvar), &
                             ddy_mean     (1:lauxp%nvar), &
                             ddy_rms      (1:lauxp%nvar), &
                             ddy_mean_rqq (1:lauxp%nvar))

  end subroutine largo_preStep_sEta


  ! Generic interface prestep subroutine
  subroutine largo_preStep_baseflow(lcp,     base, ddy_base, &
                                         ddt_base, ddx_base, &
                                         src_base) bind(C)

    real(WP), dimension(*), intent(in)  ::     base
    real(WP), dimension(*), intent(in)  :: ddy_base
    real(WP), dimension(*), intent(in)  :: ddt_base
    real(WP), dimension(*), intent(in)  :: ddx_base
    real(WP), dimension(*), intent(in)  :: src_base
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    if (associated(lauxp%largo_prestep_baseflow)) then
      call lauxp%largo_prestep_baseflow(lauxp%cp,               &
                                            base (1:lauxp%nvar), &
                                        ddy_base (1:lauxp%nvar), &
                                        ddt_base (1:lauxp%nvar), &
                                        ddx_base (1:lauxp%nvar), &
                                        src_base (1:lauxp%neq))
    end if

  end subroutine largo_preStep_baseflow


  ! Generic interface prestep subroutine
  subroutine largo_init_wall_baseflow(lcp, wall_base, wall_ddy_base, &
                                           wall_ddt_base, wall_ddx_base, &
                                           wall_src_base) bind(C)

    real(WP), dimension(*), intent(in)  ::     wall_base
    real(WP), dimension(*), intent(in)  :: wall_ddy_base
    real(WP), dimension(*), intent(in)  :: wall_ddt_base
    real(WP), dimension(*), intent(in)  :: wall_ddx_base
    real(WP), dimension(*), intent(in)  :: wall_src_base
    type(largo_ptr), value              :: lcp
    type(largo_type), pointer           :: lauxp

    call c_f_pointer(lcp, lauxp)
    if (associated(lauxp%largo_init_wall_baseflow)) then
      call lauxp%largo_init_wall_baseflow(lauxp%cp,               &
                                            wall_base (1:lauxp%nvar), &
                                        wall_ddy_base (1:lauxp%nvar), &
                                        wall_ddt_base (1:lauxp%nvar), &
                                        wall_ddx_base (1:lauxp%nvar), &
                                        wall_src_base (1:lauxp%neq))
    end if

  end subroutine largo_init_wall_baseflow


#define DECLARE_GENERIC_SUBROUTINE(token)token (lcp, A, B, src) bind(C);\
  type(largo_ptr), value                 :: lcp;\
  real(WP)       , intent(in), value     :: A, B;\
  real(WP)       , intent(inout)         :: src;\
  type(largo_type), pointer              :: lauxp;\
  call c_f_pointer(lcp, lauxp)

  ! Generic interface continuity_sEtaMean subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_continuity_sEtaMean)
    call lauxp%largo_continuity_mean(lauxp%cp, A, B, src)
  end subroutine largo_continuity_sEtaMean

  ! Generic interface xMomentum_sEtaMean subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_xMomentum_sEtaMean)
    call lauxp%largo_xmomentum_mean(lauxp%cp, A, B, src)
  end subroutine largo_xMomentum_sEtaMean

  ! Generic interface yMomentum_sEtaMean subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_yMomentum_sEtaMean)
    call lauxp%largo_ymomentum_mean(lauxp%cp, A, B, src)
  end subroutine largo_yMomentum_sEtaMean

  ! Generic interface zMomentum_sEtaMean subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_zMomentum_sEtaMean)
    call lauxp%largo_zmomentum_mean(lauxp%cp, A, B, src)
  end subroutine largo_zMomentum_sEtaMean

  ! Generic interface energy_sEtaMean subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_energy_sEtaMean)
    call lauxp%largo_energy_mean(lauxp%cp, A, B, src)
  end subroutine largo_energy_sEtaMean

  ! Generic interface ispecies_sEtaMean subroutine
  subroutine largo_ispecies_sEtaMean (lcp, A, B, src, is) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), intent(inout)                :: src
    integer(c_int), intent(in), value      :: is
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_ispecies_mean(lauxp%cp, A, B, src, is)
  end subroutine largo_ispecies_sEtaMean

  ! Generic interface species_sEtaMean subroutine
  subroutine largo_species_sEtaMean (lcp, A, B, srcvec) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), dimension(*), intent(inout)  :: srcvec
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_species_mean(lauxp%cp, A, B, srcvec(1:lauxp%ns))
  end subroutine largo_species_sEtaMean

  ! Generic interface continuity_sEtaRms subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_continuity_sEtaRms)
    call lauxp%largo_continuity_rms(lauxp%cp, A, B, src)
  end subroutine largo_continuity_sEtaRms

  ! Generic interface xMomentum_sEtaRms subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_xMomentum_sEtaRms)
    call lauxp%largo_xmomentum_rms(lauxp%cp, A, B, src)
  end subroutine largo_xMomentum_sEtaRms

  ! Generic interface yMomentum_sEtaRms subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_yMomentum_sEtaRms)
    call lauxp%largo_ymomentum_rms(lauxp%cp, A, B, src)
  end subroutine largo_yMomentum_sEtaRms

  ! Generic interface zMomentum_sEtaRms subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_zMomentum_sEtaRms)
    call lauxp%largo_zmomentum_rms(lauxp%cp, A, B, src)
  end subroutine largo_zMomentum_sEtaRms

  ! Generic interface energy_sEtaRms subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_energy_sEtaRms)
    call lauxp%largo_energy_rms(lauxp%cp, A, B, src)
  end subroutine largo_energy_sEtaRms

  ! Generic interface ispecies_sEtaRms subroutine
  subroutine largo_ispecies_sEtaRms (lcp, A, B, src, is) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), intent(inout)                :: src
    integer(c_int), intent(in), value      :: is
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_ispecies_rms(lauxp%cp, A, B, src, is)
  end subroutine largo_ispecies_sEtaRms

  ! Generic interface species_sEtaRms subroutine
  subroutine largo_species_sEtaRms (lcp, A, B, srcvec) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), dimension(*), intent(inout)  :: srcvec
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_species_rms(lauxp%cp, A, B, srcvec(1:lauxp%ns))
  end subroutine largo_species_sEtaRms

  ! Generic interface continuity_sEta wrapper subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_continuity_sEta)
    call lauxp%largo_continuity(lauxp%cp, A, B, src)
  end subroutine largo_continuity_sEta

  ! Generic interface xMomentum_sEta wrapper subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_xMomentum_sEta)
    call lauxp%largo_xmomentum(lauxp%cp, A, B, src)
  end subroutine largo_xMomentum_sEta

  ! Generic interface yMomentum_sEta wrapper subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_yMomentum_sEta)
    call lauxp%largo_ymomentum(lauxp%cp, A, B, src)
  end subroutine largo_yMomentum_sEta

  ! Generic interface zMomentum_sEta wrapper subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_zMomentum_sEta)
    call lauxp%largo_zmomentum(lauxp%cp, A, B, src)
  end subroutine largo_zMomentum_sEta

  ! Generic interface energy_sEta wrapper subroutine
  subroutine DECLARE_GENERIC_SUBROUTINE(largo_energy_sEta)
    call lauxp%largo_energy(lauxp%cp, A, B, src)
  end subroutine largo_energy_sEta

  ! Generic interface species_sEta wrapper subroutine
  subroutine largo_species_sEta (lcp, A, B, srcvec) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), dimension(*), intent(inout)  :: srcvec
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_species(lauxp%cp, A, B, srcvec(1:lauxp%ns))
  end subroutine largo_species_sEta

  ! Generic interface sEtaMean wrapper subroutine
  subroutine largo_sEtaMean (lcp, A, B, srcvec) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), dimension(*), intent(inout)  :: srcvec
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_all_sources_mean(lauxp%cp, A, B, srcvec(1:lauxp%neq))
  end subroutine largo_sEtaMean

  ! Generic interface sEta wrapper subroutine
  subroutine largo_sEta (lcp, A, B, srcvec) bind(C)
    type(largo_ptr), value                 :: lcp
    real(WP), intent(in), value            :: A, B
    real(WP), dimension(*), intent(inout)  :: srcvec
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    call lauxp%largo_all_sources(lauxp%cp, A, B, srcvec(1:lauxp%neq))
  end subroutine largo_sEta


  ! Convenience method to get the data workspace pointer from the generic
  ! pointer
  subroutine largo_get_workspace_pointer(lcp, lwcp)

    type(largo_ptr), intent(in)            :: lcp
    type(largo_workspace_ptr), intent(out) :: lwcp
    type(largo_type), pointer              :: lauxp

    call c_f_pointer(lcp, lauxp)
    lwcp = lauxp%cp
  end subroutine largo_get_workspace_pointer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>Copy the contents of a C-style string to a Fortran-style string
!!On non-NULL \c c_string, copy as much of its contents as possible into
!!\c f_string.  Return success if-and-only-if the entire copy succeeded.
!!Otherwise return failure.  In addition to returning failure, copying
!!a NULL pointer clears \c f_string.
!NOTE: Taken from the esio library
  function largo_c_f_stringcopy (c_string, f_string) result (success)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_associated, c_f_pointer, &
                                           c_char, c_null_char

    implicit none

    logical                               :: success
    type(c_ptr),      intent(in)          :: c_string
    character(len=*), intent(inout)       :: f_string
    character(len=1,kind=c_char), pointer :: tmp_str(:)
    integer                               :: i, n(1)

    success = .false.
    if (c_associated(c_string)) then
      n(1) = len(f_string)
      call c_f_pointer(c_string, tmp_str, n)
      do i = 1, n(1)
        if (tmp_str(i) == c_null_char) then
          f_string(i:) = ''   ! Clear any remaining Fortran string storage
          success = .true.    ! Success because full copy succeeded
          exit
        end if
        f_string(i:i) = tmp_str(i)
      end do
    else
      f_string = ''           ! Clear entire Fortran string on NULL pointer
    end if

  end function largo_c_f_stringcopy

end module largo
