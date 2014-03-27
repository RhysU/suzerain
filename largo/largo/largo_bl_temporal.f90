module largo_BL_temporal

! Use ISO_C_BINDING to expose C-friendly API
  use, intrinsic :: iso_c_binding, only: c_associated,   &
                                         WP => c_double, &
                                         c_int,          &
                                         c_loc,          &
                                         c_f_pointer,    &
                                         c_null_ptr,     &
                                         largo_workspace_ptr =>c_ptr

  implicit none ! is in effect throughout entire module

  private


  ! interface for rans prestep function
  abstract interface
    subroutine prestep_setamean_rans(cp, y, mean, ddy_mean)
      import
      real(WP), intent(in)                  :: y
      real(WP), dimension(*), intent(in)    :: mean
      real(WP), dimension(*), intent(in)    :: ddy_mean
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_setamean_rans
  end interface


  ! interface for rans source function
  abstract interface
    subroutine sourcevec_rans(cp, A, B, srcvec)
      import
      type(largo_workspace_ptr), intent(in)  :: cp
      real(WP), intent(in)                   :: A, B
      real(WP), dimension(*), intent(inout)  :: srcvec
    end subroutine sourcevec_rans
  end interface


  ! largo workspace type declaration
  type :: largo_BL_temporal_workspace_type

    ! Threshold chosen from experience with mean densities
    ! like 0.0025 requiring cutoffs like 1e-10_WP.
    real(WP) :: cutoff_wrt_rho = 4e-8_WP

    real(WP) :: grt_delta   = 1.0_WP

    real(WP) :: Ts_rho      = 0.0_WP
    real(WP) :: Ts_rhoU     = 0.0_WP
    real(WP) :: Ts_rhoV     = 0.0_WP
    real(WP) :: Ts_rhoW     = 0.0_WP
    real(WP) :: Ts_rhoE     = 0.0_WP

    real(WP) :: ygrms_rho   = 0.0_WP
    real(WP) :: ygrms_rhoU  = 0.0_WP
    real(WP) :: ygrms_rhoV  = 0.0_WP
    real(WP) :: ygrms_rhoW  = 0.0_WP
    real(WP) :: ygrms_rhoE  = 0.0_WP

    real(WP) :: mean_rho    = 0.0_WP
    real(WP) :: mean_rhoU   = 0.0_WP
    real(WP) :: mean_rhoV   = 0.0_WP
    real(WP) :: mean_rhoW   = 0.0_WP
    real(WP) :: mean_rhoE   = 0.0_WP

    real(WP) :: fluc_rho    = 0.0_WP
    real(WP) :: fluc_rhoU   = 0.0_WP
    real(WP) :: fluc_rhoV   = 0.0_WP
    real(WP) :: fluc_rhoW   = 0.0_WP
    real(WP) :: fluc_rhoE   = 0.0_WP

    real(WP) :: TsRms_rho   = 0.0_WP
    real(WP) :: TsRms_rhoU  = 0.0_WP
    real(WP) :: TsRms_rhoV  = 0.0_WP
    real(WP) :: TsRms_rhoW  = 0.0_WP
    real(WP) :: TsRms_rhoE  = 0.0_WP

    real(WP), allocatable, dimension(:) :: Ts_rhos
    real(WP), allocatable, dimension(:) :: ygrms_rhos
    real(WP), allocatable, dimension(:) :: mean_rhos
    real(WP), allocatable, dimension(:) :: fluc_rhos
    real(WP), allocatable, dimension(:) :: TsRms_rhos

    ! baseflow variables
    real(WP) :: base_rho        = 0.0_WP
    real(WP) :: base_rhoU       = 0.0_WP
    real(WP) :: base_rhoV       = 0.0_WP
    real(WP) :: base_rhoW       = 0.0_WP
    real(WP) :: base_rhoE       = 0.0_WP

    real(WP) :: ddy_base_rho    = 0.0_WP
    real(WP) :: ddy_base_rhoU   = 0.0_WP
    real(WP) :: ddy_base_rhoV   = 0.0_WP
    real(WP) :: ddy_base_rhoW   = 0.0_WP
    real(WP) :: ddy_base_rhoE   = 0.0_WP

    real(WP) :: ddt_base_rho    = 0.0_WP
    real(WP) :: ddt_base_rhoU   = 0.0_WP
    real(WP) :: ddt_base_rhoV   = 0.0_WP
    real(WP) :: ddt_base_rhoW   = 0.0_WP
    real(WP) :: ddt_base_rhoE   = 0.0_WP

    real(WP) :: src_base_rho    = 0.0_WP
    real(WP) :: src_base_rhoU   = 0.0_WP
    real(WP) :: src_base_rhoV   = 0.0_WP
    real(WP) :: src_base_rhoW   = 0.0_WP
    real(WP) :: src_base_rhoE   = 0.0_WP

    real(WP) :: grt_DA_rho      = 0.0_WP
    real(WP) :: grt_DA_rhoU     = 0.0_WP
    real(WP) :: grt_DA_rhoV     = 0.0_WP
    real(WP) :: grt_DA_rhoW     = 0.0_WP
    real(WP) :: grt_DA_rhoE     = 0.0_WP

    real(WP) :: grt_DA_rms_rho  = 0.0_WP
    real(WP) :: grt_DA_rms_rhoU = 0.0_WP
    real(WP) :: grt_DA_rms_rhoV = 0.0_WP
    real(WP) :: grt_DA_rms_rhoW = 0.0_WP
    real(WP) :: grt_DA_rms_rhoE = 0.0_WP

    real(WP), allocatable, dimension(:) :: base_rhos
    real(WP), allocatable, dimension(:) :: ddy_base_rhos
    real(WP), allocatable, dimension(:) :: ddt_base_rhos
    real(WP), allocatable, dimension(:) :: src_base_rhos
    real(WP), allocatable, dimension(:) :: grt_DA_rhos
    real(WP), allocatable, dimension(:) :: grt_DA_rms_rhos

    ! RANS variables
    real(WP), allocatable, dimension(:) :: Ts_tvar
    real(WP), allocatable, dimension(:) :: mean_tvar
    real(WP), allocatable, dimension(:) :: grt_DA_tvar
    real(WP), allocatable, dimension(:) :: grt_DA_rms_tvar

  end type largo_BL_temporal_workspace_type

  procedure(prestep_setamean_rans), pointer :: largo_BL_temporal_prestep_sEtaMean_rans => NULL()
  procedure(sourcevec_rans),        pointer :: largo_BL_temporal_sEta_rans             => NULL()

  ! Indices
  integer(c_int), parameter :: irho  = 1
  integer(c_int), parameter :: irhoU = 2
  integer(c_int), parameter :: irhoV = 3
  integer(c_int), parameter :: irhoW = 4
  integer(c_int), parameter :: irhoE = 5

  ! Number of equations
  integer(c_int) :: neq_ = 0

  ! Number of species
  integer(c_int) :: ns_  = 0

  ! Number of tubulence variables
  integer(c_int) :: ntvar_  = 0

  ! Index of first tubulence variable
  integer(c_int) :: itvar0_  = 0

  public  :: largo_BL_temporal_allocate
  public  :: largo_BL_temporal_deallocate
  public  :: largo_BL_temporal_init
  public  :: largo_BL_temporal_preStep_sEta
  public  :: largo_BL_temporal_preStep_sEta_innery
  public  :: largo_BL_temporal_preStep_sEta_innerxz
  public  :: largo_BL_temporal_preStep_sEtaMean
  public  :: largo_BL_temporal_continuity_sEtaMean
  public  :: largo_BL_temporal_xMomentum_sEtaMean
  public  :: largo_BL_temporal_yMomentum_sEtaMean
  public  :: largo_BL_temporal_zMomentum_sEtaMean
  public  :: largo_BL_temporal_energy_sEtaMean
  public  :: largo_BL_temporal_ispecies_sEtaMean
  public  :: largo_BL_temporal_species_sEtaMean
  public  :: largo_BL_temporal_continuity_sEtaRms
  public  :: largo_BL_temporal_xMomentum_sEtaRms
  public  :: largo_BL_temporal_yMomentum_sEtaRms
  public  :: largo_BL_temporal_zMomentum_sEtaRms
  public  :: largo_BL_temporal_energy_sEtaRms
  public  :: largo_BL_temporal_ispecies_sEtaRms
  public  :: largo_BL_temporal_species_sEtaRms

  private :: largo_BL_temporal_preStep_sEtaRms

  public  :: largo_BL_temporal_continuity_sEta
  public  :: largo_BL_temporal_xMomentum_sEta
  public  :: largo_BL_temporal_yMomentum_sEta
  public  :: largo_BL_temporal_zMomentum_sEta
  public  :: largo_BL_temporal_energy_sEta
  public  :: largo_BL_temporal_species_sEta

  public  :: largo_BL_temporal_sEtaMean
  public  :: largo_BL_temporal_sEta

  public  :: largo_BL_temporal_preStep_baseflow

  public  :: largo_BL_temporal_get_ntvar_rans
  public  :: largo_BL_temporal_init_rans
  public  :: largo_BL_temporal_prestep_sEtaMean_rans
  public  :: largo_BL_temporal_sEta_rans

contains

  subroutine largo_BL_temporal_allocate(cp, neq, ns, ntvar, ransmodel)

    type(largo_workspace_ptr), intent(out) :: cp
    ! neq=number of equations, might be needed later
    integer(c_int), intent(in)   :: neq
    integer(c_int), intent(in)   :: ns    ! number of species
    integer(c_int), intent(in)   :: ntvar ! number of turbulence variables
    character(len=*)                                :: ransmodel
    type(largo_BL_temporal_workspace_type), pointer :: auxp

    ! Allocate derived type variable
    allocate(auxp)

    ! Initialize values
    neq_    = neq
    ns_     = ns
    ntvar_  = ntvar
    itvar0_ = 5 + ns_ + 1

    ! Allocate arrays for species
    if (ns_ > 0) then
      allocate(auxp%Ts_rhos        (1:ns_))
      allocate(auxp%ygrms_rhos     (1:ns_))
      allocate(auxp%mean_rhos      (1:ns_))
      allocate(auxp%fluc_rhos      (1:ns_))
      allocate(auxp%TsRms_rhos     (1:ns_))

      allocate(auxp%base_rhos      (1:ns_))
      allocate(auxp%ddy_base_rhos  (1:ns_))
      allocate(auxp%ddt_base_rhos  (1:ns_))
      allocate(auxp%src_base_rhos  (1:ns_))
      allocate(auxp%grt_DA_rhos    (1:ns_))
      allocate(auxp%grt_DA_rms_rhos(1:ns_))

      auxp%base_rhos     = 0.0_WP
      auxp%ddy_base_rhos = 0.0_WP
      auxp%ddt_base_rhos = 0.0_WP
      auxp%src_base_rhos = 0.0_WP
      auxp%grt_DA_rhos    = 0.0_WP
      auxp%grt_DA_rms_rhos= 0.0_WP
    end if

    ! Check number of turbulence variables
    select case (trim(ransmodel))
    case ("laminar", "dns")
      ! FIXME: if (ntvar_ /= 0) point to an error
    case ("turbulent_viscosity")
      ! FIXME: if (ntvar_ /= 1) point to an error
    case ("k_epsilon")
      ! FIXME: if (ntvar_ /= 2) point to an error
    case ("k_omega")
      ! FIXME: if (ntvar_ /= 2) point to an error
    case default
      ! FIXME: RANS model not defined
    end select

    ! Allocate arrays for turbulent variables
    if (ntvar_ > 0) then
      allocate(auxp%Ts_tvar         (1:ntvar_))
      allocate(auxp%mean_tvar       (1:ntvar_))
      allocate(auxp%grt_DA_tvar     (1:ntvar_))
      allocate(auxp%grt_DA_rms_tvar (1:ntvar_))

      auxp%Ts_tvar         = 0.0_WP
      auxp%mean_tvar       = 0.0_WP
      auxp%grt_DA_tvar     = 0.0_WP
      auxp%grt_DA_rms_tvar = 0.0_WP
    end if

    ! Initialize turbulence variables function pointers
    largo_BL_temporal_prestep_sEtaMean_rans => largo_bl_temporal_prestep_sEtaMean_rans_generic
    largo_BL_temporal_sEta_rans             => largo_bl_temporal_sEta_rans_generic

    ! Get C pointer from Fortran pointer
    cp = c_loc(auxp)

  end subroutine largo_BL_temporal_allocate


  subroutine largo_BL_temporal_deallocate(cp)

    type(largo_workspace_ptr), intent(inout)  :: cp
    type(largo_BL_temporal_workspace_type), pointer :: auxp

    call c_f_pointer(cp, auxp)

    ! Deallocate arrays for species
    if (allocated(auxp%Ts_rhos         ))  deallocate(auxp%Ts_rhos        )
    if (allocated(auxp%ygrms_rhos      ))  deallocate(auxp%ygrms_rhos     )
    if (allocated(auxp%mean_rhos       ))  deallocate(auxp%mean_rhos      )
    if (allocated(auxp%fluc_rhos       ))  deallocate(auxp%fluc_rhos      )
    if (allocated(auxp%TsRms_rhos      ))  deallocate(auxp%TsRms_rhos     )

    if (allocated(auxp%base_rhos       ))  deallocate(auxp%base_rhos      )
    if (allocated(auxp%ddy_base_rhos   ))  deallocate(auxp%ddy_base_rhos  )
    if (allocated(auxp%ddt_base_rhos   ))  deallocate(auxp%ddt_base_rhos  )
    if (allocated(auxp%src_base_rhos   ))  deallocate(auxp%src_base_rhos  )
    if (allocated(auxp%grt_DA_rhos     ))  deallocate(auxp%grt_DA_rhos    )
    if (allocated(auxp%grt_DA_rms_rhos ))  deallocate(auxp%grt_DA_rms_rhos)

    ! Deallocate arrays for RANS turbulence variables
    if (allocated(auxp%Ts_tvar         ))  deallocate(auxp%Ts_tvar        )
    if (allocated(auxp%mean_tvar       ))  deallocate(auxp%mean_tvar      )
    if (allocated(auxp%grt_DA_tvar     ))  deallocate(auxp%grt_DA_tvar    )
    if (allocated(auxp%grt_DA_rms_tvar ))  deallocate(auxp%grt_DA_rms_tvar)

    ! Deallocate array of derived types
    deallocate(auxp)

    ! Nullify C pointer
    cp = c_null_ptr

  end subroutine largo_BL_temporal_deallocate


  ! get number of turbulence variables
  subroutine largo_BL_temporal_get_ntvar_rans(cp, ntvar)

    ! largo workspace C pointer
    type(largo_workspace_ptr), intent(in)            :: cp
    integer(c_int)           , intent(out)           :: ntvar

!!$     ! Get Fortran pointer from C pointer
!!$     call c_f_pointer(cp, auxp)

    ntvar = ntvar_

  end subroutine largo_BL_temporal_get_ntvar_rans


  subroutine largo_BL_temporal_init(cp, grt_delta, grt_DA, grt_DA_rms)

    real(WP), intent(in)                  :: grt_delta
    real(WP), dimension(*), intent(in)    :: grt_DA
    real(WP), dimension(*), intent(in)    :: grt_DA_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Set growth rates
    auxp%grt_delta   = grt_delta

    auxp%grt_DA_rho  = grt_DA(irho )
    auxp%grt_DA_rhoU = grt_DA(irhoU)
    auxp%grt_DA_rhoV = grt_DA(irhoV)
    auxp%grt_DA_rhoW = grt_DA(irhoW)
    auxp%grt_DA_rhoE = grt_DA(irhoE)

    auxp%grt_DA_rms_rho  = grt_DA_rms(irho )
    auxp%grt_DA_rms_rhoU = grt_DA_rms(irhoU)
    auxp%grt_DA_rms_rhoV = grt_DA_rms(irhoV)
    auxp%grt_DA_rms_rhoW = grt_DA_rms(irhoW)
    auxp%grt_DA_rms_rhoE = grt_DA_rms(irhoE)

    do is=1, ns_
      auxp%grt_DA_rhos    (is) = grt_DA    (5+is)
      auxp%grt_DA_rms_rhos(is) = grt_DA_rms(5+is)
    end do

  end subroutine largo_BL_temporal_init


  ! RANS initialization
  subroutine largo_BL_temporal_init_rans(cp, grt_delta, grt_DA, grt_DA_rms)

    ! largo workspace C pointer
    real(WP), intent(in)                  :: grt_delta
    real(WP), dimension(*), intent(in)    :: grt_DA
    real(WP), dimension(*), intent(in)    :: grt_DA_rms
    type(largo_workspace_ptr), intent(in)            :: cp
    type(largo_BL_temporal_workspace_type), pointer  :: auxp
    integer(c_int) :: it

!!$     call largo_BL_temporal_init(cp, grt_delta          , &
!!$                                     grt_DA    (1:neq_), &
!!$                                     grt_DA_rms(1:neq_))

    call largo_BL_temporal_init(cp, grt_delta, grt_DA, grt_DA_rms)

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    do it=1, ntvar_
      auxp%grt_DA_tvar     (it) = grt_DA     (itvar0_ + it -1)
      auxp%grt_DA_rms_tvar (it) = grt_DA_rms (itvar0_ + it -1)
    end do

  end subroutine largo_BL_temporal_init_rans


  subroutine largo_BL_temporal_preStep_baseflow(cp,     base, ddy_base, &
                                          ddt_base, ddx_base, src_base)

    real(WP), dimension(*), intent(in)    :: base
    real(WP), dimension(*), intent(in)    :: ddy_base
    real(WP), dimension(*), intent(in)    :: ddt_base
    real(WP), dimension(*), intent(in)    :: ddx_base
    real(WP), dimension(*), intent(in)    :: src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store baseflow information
    auxp%base_rho  = base(irho )
    auxp%base_rhoU = base(irhoU)
    auxp%base_rhoV = base(irhoV)
    auxp%base_rhoW = base(irhoW)
    auxp%base_rhoE = base(irhoE)

    auxp%ddy_base_rho  = ddy_base(irho )
    auxp%ddy_base_rhoU = ddy_base(irhoU)
    auxp%ddy_base_rhoV = ddy_base(irhoV)
    auxp%ddy_base_rhoW = ddy_base(irhoW)
    auxp%ddy_base_rhoE = ddy_base(irhoE)

    auxp%ddt_base_rho  = ddt_base(irho )
    auxp%ddt_base_rhoU = ddt_base(irhoU)
    auxp%ddt_base_rhoV = ddt_base(irhoV)
    auxp%ddt_base_rhoW = ddt_base(irhoW)
    auxp%ddt_base_rhoE = ddt_base(irhoE)

    auxp%src_base_rho  = src_base(irho )
    auxp%src_base_rhoU = src_base(irhoU)
    auxp%src_base_rhoV = src_base(irhoV)
    auxp%src_base_rhoW = src_base(irhoW)
    auxp%src_base_rhoE = src_base(irhoE)

    do is=1, ns_
      auxp%base_rhos    (is) =     base(5+is)
      auxp%ddy_base_rhos(is) = ddy_base(5+is)
      auxp%ddt_base_rhos(is) = ddt_base(5+is)
      auxp%src_base_rhos(is) = src_base(5+is)
    end do

  end subroutine largo_BL_temporal_preStep_baseflow


  subroutine largo_BL_temporal_preStep_sEtaMean(cp, y, mean, ddy_mean)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: ddy_mean
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! These ones depend on y only
    auxp%mean_rho  = mean(irho )
    auxp%mean_rhoU = mean(irhoU)
    auxp%mean_rhoV = mean(irhoV)
    auxp%mean_rhoW = mean(irhoW)
    auxp%mean_rhoE = mean(irhoE)

    auxp%Ts_rho  = - auxp%ddt_base_rho  - auxp%grt_DA_rho  * (mean(irho )-auxp%base_rho ) + y * auxp%grt_delta * (ddy_mean(irho ) - auxp%ddy_base_rho ) + auxp%src_base_rho
    auxp%Ts_rhoU = - auxp%ddt_base_rhoU - auxp%grt_DA_rhoU * (mean(irhoU)-auxp%base_rhoU) + y * auxp%grt_delta * (ddy_mean(irhoU) - auxp%ddy_base_rhoU) + auxp%src_base_rhoU
    auxp%Ts_rhoV = - auxp%ddt_base_rhoV - auxp%grt_DA_rhoV * (mean(irhoV)-auxp%base_rhoV) + y * auxp%grt_delta * (ddy_mean(irhoV) - auxp%ddy_base_rhoV) + auxp%src_base_rhoV
    auxp%Ts_rhoW = - auxp%ddt_base_rhoW - auxp%grt_DA_rhoW * (mean(irhoW)-auxp%base_rhoW) + y * auxp%grt_delta * (ddy_mean(irhoW) - auxp%ddy_base_rhoW) + auxp%src_base_rhoW
    auxp%Ts_rhoE = - auxp%ddt_base_rhoE - auxp%grt_DA_rhoE * (mean(irhoE)-auxp%base_rhoE) + y * auxp%grt_delta * (ddy_mean(irhoE) - auxp%ddy_base_rhoE) + auxp%src_base_rhoE

    do is=1, ns_
      auxp%mean_rhos(is) = mean(5+is)
      auxp%Ts_rhos(is)  = - auxp%ddt_base_rhos(is)  - auxp%grt_DA_rhos(is)  * (mean(5+is)-auxp%base_rhos(is)) + y * auxp%grt_delta * (ddy_mean(5+is) - auxp%ddy_base_rhos(is)) + auxp%src_base_rhos(is)
    end do

  end subroutine largo_BL_temporal_preStep_sEtaMean


  subroutine largo_BL_temporal_preStep_sEtaRms(cp, y, rms, ddy_rms)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: rms
    real(WP), dimension(*), intent(in)    :: ddy_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp
    real(WP) :: eps

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Obtain epsilon-like cutoff scaled by local mean density
    eps = auxp%mean_rho * auxp%cutoff_wrt_rho

    ! These ones depend on y only
    auxp%ygrms_rho  = 0.0_WP
    if (rms(irho ) > eps)  &
       auxp%ygrms_rho  = y * auxp%grt_delta * ddy_rms(irho )/(rms(irho ))

    auxp%ygrms_rhoU = 0.0_WP
    if (rms(irhoU) > eps)  &
       auxp%ygrms_rhoU = y * auxp%grt_delta * ddy_rms(irhoU)/(rms(irhoU))

    auxp%ygrms_rhoV = 0.0_WP
    if (rms(irhoV) > eps)  &
       auxp%ygrms_rhoV = y * auxp%grt_delta * ddy_rms(irhoV)/(rms(irhoV))

    auxp%ygrms_rhoW = 0.0_WP
    if (rms(irhoW) > eps)  &
       auxp%ygrms_rhoW = y * auxp%grt_delta * ddy_rms(irhoW)/(rms(irhoW))

    auxp%ygrms_rhoE = 0.0_WP
    if (rms(irhoE) > eps)  &
       auxp%ygrms_rhoE = y * auxp%grt_delta * ddy_rms(irhoE)/(rms(irhoE))

    do is=1, ns_
      auxp%ygrms_rhos(is) = 0.0_WP
      if (rms(5+is) > eps)  &
         auxp%ygrms_rhos(is) = y * auxp%grt_delta * ddy_rms(5+is)/(rms(5+is))
    end do

  end subroutine largo_BL_temporal_preStep_sEtaRms


  subroutine largo_BL_temporal_preStep_sEta_innerxz(cp, qflow)

    real(WP), dimension(*), intent(in)    :: qflow
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    auxp%fluc_rho  = qflow(irho ) - auxp%mean_rho
    auxp%fluc_rhoU = qflow(irhoU) - auxp%mean_rhoU
    auxp%fluc_rhoV = qflow(irhoV) - auxp%mean_rhoV
    auxp%fluc_rhoW = qflow(irhoW) - auxp%mean_rhoW
    auxp%fluc_rhoE = qflow(irhoE) - auxp%mean_rhoE

    auxp%TsRms_rho  = auxp%fluc_rho  * (- auxp%grt_DA_rms_rho  + auxp%ygrms_rho )
    auxp%TsRms_rhoU = auxp%fluc_rhoU * (- auxp%grt_DA_rms_rhoU + auxp%ygrms_rhoU)
    auxp%TsRms_rhoV = auxp%fluc_rhoV * (- auxp%grt_DA_rms_rhoV + auxp%ygrms_rhoV)
    auxp%TsRms_rhoW = auxp%fluc_rhoW * (- auxp%grt_DA_rms_rhoW + auxp%ygrms_rhoW)
    auxp%TsRms_rhoE = auxp%fluc_rhoE * (- auxp%grt_DA_rms_rhoE + auxp%ygrms_rhoE)

    do is=1, ns_
      auxp%fluc_rhos(is)  = qflow(5+is) - auxp%mean_rhos(is)
      auxp%TsRms_rhos(is) = auxp%fluc_rhos(is) * (- auxp%grt_DA_rms_rhos(is) + auxp%ygrms_rhos(is))
    end do

  end subroutine largo_BL_temporal_preStep_sEta_innerxz


  subroutine largo_BL_temporal_preStep_sEta_innery(cp, y,                  &
                                              mean,     rms,     mean_rqq, &
                                          ddy_mean, ddy_rms, ddy_mean_rqq)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: rms
    real(WP), dimension(*), intent(in)    :: mean_rqq        ! not used
    real(WP), dimension(*), intent(in)    :: ddy_mean
    real(WP), dimension(*), intent(in)    :: ddy_rms
    real(WP), dimension(*), intent(in)    :: ddy_mean_rqq    ! not used
    type(largo_workspace_ptr), intent(in) :: cp

    call largo_BL_temporal_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_temporal_preStep_sEtaRms(cp, y, rms, ddy_rms)

  end subroutine largo_BL_temporal_preStep_sEta_innery


  subroutine largo_BL_temporal_preStep_sEta(cp, y, qflow,                    &
                                                mean,     rms,     mean_rqq, &
                                            ddy_mean, ddy_rms, ddy_mean_rqq)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: qflow
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: rms
    real(WP), dimension(*), intent(in)    :: mean_rqq        ! not used
    real(WP), dimension(*), intent(in)    :: ddy_mean
    real(WP), dimension(*), intent(in)    :: ddy_rms
    real(WP), dimension(*), intent(in)    :: ddy_mean_rqq    ! not used
    type(largo_workspace_ptr), intent(in) :: cp

    ! Delegate to Y then XZ routines to ensure their behavior reproduced
    call largo_BL_temporal_preStep_sEta_innery(cp, y,                           &
                                                   mean,     rms,     mean_rqq, &
                                               ddy_mean, ddy_rms, ddy_mean_rqq)
    call largo_BL_temporal_preStep_sEta_innerxz(cp, qflow)

  end subroutine largo_BL_temporal_preStep_sEta


  subroutine largo_BL_temporal_preStep_sEtaMean_rans_generic(cp, y, mean, ddy_mean)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: ddy_mean
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_workspace_type), pointer   :: auxp
    integer(c_int) :: it

    call largo_BL_temporal_preStep_sEtaMean(cp, y, mean, ddy_mean)

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    do it=1, ntvar_
      auxp%mean_tvar(it) = mean(itvar0_+it-1)
      auxp%Ts_tvar(it)  = - auxp%grt_DA_tvar(it)  * auxp%mean_tvar(it) + y * auxp%grt_delta * ddy_mean(itvar0_+it-1)
    end do

  end subroutine largo_BL_temporal_preStep_sEtaMean_rans_generic


#define DECLARE_SUBROUTINE(token)token (cp, A, B, src);\
  type(largo_workspace_ptr), intent(in)  :: cp;\
  real(WP)       , intent(in)            :: A, B;\
  real(WP)       , intent(inout)         :: src;\
  type(largo_BL_temporal_workspace_type), pointer    :: auxp;\
  call c_f_pointer(cp, auxp);\


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_continuity_sEtaMean)
    src = A * src + B * auxp%Ts_rho
  end subroutine largo_BL_temporal_continuity_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_xMomentum_sEtaMean)
    src = A * src + B * auxp%Ts_rhoU
  end subroutine largo_BL_temporal_xMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_yMomentum_sEtaMean)
    src = A * src + B * auxp%Ts_rhoV
  end subroutine largo_BL_temporal_yMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_zMomentum_sEtaMean)
    src = A * src + B * auxp%Ts_rhoW
  end subroutine largo_BL_temporal_zMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_energy_sEtaMean)
    src = A * src + B * auxp%Ts_rhoE
  end subroutine largo_BL_temporal_energy_sEtaMean


  subroutine largo_BL_temporal_ispecies_sEtaMean (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)   :: cp
    real(WP)       , intent(in)             :: A, B
    real(WP)       , intent(inout)          :: src
    type(largo_BL_temporal_workspace_type), pointer     :: auxp
    integer(c_int), intent(in)              :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * auxp%Ts_rhos(is)
  end subroutine largo_BL_temporal_ispecies_sEtaMean


  subroutine largo_BL_temporal_species_sEtaMean (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_temporal_ispecies_sEtaMean (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_species_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_continuity_sEtaRms)
    src = A * src + B * auxp%TsRms_rho
  end subroutine largo_BL_temporal_continuity_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_xMomentum_sEtaRms)
    src = A * src + B * auxp%TsRms_rhoU
  end subroutine largo_BL_temporal_xMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_yMomentum_sEtaRms)
    src = A * src + B * auxp%TsRms_rhoV
  end subroutine largo_BL_temporal_yMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_zMomentum_sEtaRms)
    src = A * src + B * auxp%TsRms_rhoW
  end subroutine largo_BL_temporal_zMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_energy_sEtaRms)
    src = A * src + B * auxp%TsRms_rhoE
  end subroutine largo_BL_temporal_energy_sEtaRms


  subroutine largo_BL_temporal_ispecies_sEtaRms (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP)       , intent(inout)            :: src
    type(largo_BL_temporal_workspace_type), pointer       :: auxp
    integer(c_int), intent(in)                :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * auxp%TsRms_rhos(is)
  end subroutine largo_BL_temporal_ispecies_sEtaRms


  subroutine largo_BL_temporal_species_sEtaRms (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_temporal_ispecies_sEtaRms (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_species_sEtaRms


  subroutine largo_BL_temporal_sEta_rans_generic (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec
    type(largo_BL_temporal_workspace_type), pointer       :: auxp
    integer(c_int)                            :: it

    call largo_BL_temporal_sEtaMean (cp, A, B, srcvec) 

    call c_f_pointer(cp, auxp)
    do it = 1, ntvar_
       srcvec(itvar0_+it-1) = A * srcvec(itvar0_+it-1) + B * auxp%Ts_tvar(it)
    end do

  end subroutine largo_BL_temporal_sEta_rans_generic


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_continuity_sEta)
    call largo_BL_temporal_continuity_sEtaMean (cp,      A, B, src)
    call largo_BL_temporal_continuity_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_temporal_continuity_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_xMomentum_sEta)
    call largo_BL_temporal_xMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_temporal_xMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_temporal_xMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_yMomentum_sEta)
    call largo_BL_temporal_yMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_temporal_yMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_temporal_yMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_zMomentum_sEta)
    call largo_BL_temporal_zMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_temporal_zMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_temporal_zMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_energy_sEta)
    call largo_BL_temporal_energy_sEtaMean (cp,      A, B, src)
    call largo_BL_temporal_energy_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_temporal_energy_sEta


  subroutine largo_BL_temporal_species_sEta (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_temporal_ispecies_sEtaMean (cp,      A, B, srcvec(is), is)
      call largo_BL_temporal_ispecies_sEtaRms  (cp, 1.0_WP, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_species_sEta


  subroutine largo_BL_temporal_sEtaMean (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_temporal_continuity_sEtaMean (cp,      A, B, srcvec(irho ))
    call largo_BL_temporal_xMomentum_sEtaMean  (cp,      A, B, srcvec(irhou))
    call largo_BL_temporal_yMomentum_sEtaMean  (cp,      A, B, srcvec(irhov))
    call largo_BL_temporal_zMomentum_sEtaMean  (cp,      A, B, srcvec(irhow))
    call largo_BL_temporal_energy_sEtaMean     (cp,      A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_temporal_ispecies_sEtaMean (cp,      A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_temporal_sEtaMean


  subroutine largo_BL_temporal_sEta (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_temporal_continuity_sEtaMean (cp,      A, B, srcvec(irho ))
    call largo_BL_temporal_continuity_sEtaRms  (cp, 1.0_WP, B, srcvec(irho ))

    call largo_BL_temporal_xMomentum_sEtaMean  (cp,      A, B, srcvec(irhou))
    call largo_BL_temporal_xMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhou))

    call largo_BL_temporal_yMomentum_sEtaMean  (cp,      A, B, srcvec(irhov))
    call largo_BL_temporal_yMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhov))

    call largo_BL_temporal_zMomentum_sEtaMean  (cp,      A, B, srcvec(irhow))
    call largo_BL_temporal_zMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhow))

    call largo_BL_temporal_energy_sEtaMean     (cp,      A, B, srcvec(irhoE))
    call largo_BL_temporal_energy_sEtaRms      (cp, 1.0_WP, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_temporal_ispecies_sEtaMean (cp,      A, B, srcvec(5+is), is)
      call largo_BL_temporal_ispecies_sEtaRms  (cp, 1.0_WP, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_temporal_sEta

end module largo_BL_temporal
