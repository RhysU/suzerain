module largo_BL_spatiotemporal_consistent

  ! Use ISO_C_BINDING to expose C-friendly API through generic interface
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


  type :: largo_BL_spatiotemporal_consistent_workspace_type

    real(WP) :: gr_delta   = 1.0_WP
    real(WP) :: grt_delta   = 1.0_WP
    real(WP) :: grx_delta   = 1.0_WP

    real(WP) :: Xs_rho = 0.0_WP
    real(WP) :: Xs_U   = 0.0_WP
    real(WP) :: Xs_V   = 0.0_WP
    real(WP) :: Xs_W   = 0.0_WP
    real(WP) :: Xs_E   = 0.0_WP
    real(WP) :: Xs_P   = 0.0_WP

    real(WP) :: XsArms_rho = 0.0_WP
    real(WP) :: XsArms_U   = 0.0_WP
    real(WP) :: XsArms_V   = 0.0_WP
    real(WP) :: XsArms_W   = 0.0_WP
    real(WP) :: XsArms_E   = 0.0_WP

!!$     real(WP) :: XsFull_rho = 0.0_WP
!!$     real(WP) :: XsFull_U   = 0.0_WP
!!$     real(WP) :: XsFull_V   = 0.0_WP
!!$     real(WP) :: XsFull_W   = 0.0_WP
!!$     real(WP) :: XsFull_E   = 0.0_WP

    real(WP) :: S_rho = 0.0_WP
    real(WP) :: S_U   = 0.0_WP
    real(WP) :: S_V   = 0.0_WP
    real(WP) :: S_W   = 0.0_WP
    real(WP) :: S_E   = 0.0_WP

    real(WP) :: SArms_rho = 0.0_WP
    real(WP) :: SArms_U   = 0.0_WP
    real(WP) :: SArms_V   = 0.0_WP
    real(WP) :: SArms_W   = 0.0_WP
    real(WP) :: SArms_E   = 0.0_WP

    real(WP) :: SFull_rho = 0.0_WP
    real(WP) :: SFull_U   = 0.0_WP
    real(WP) :: SFull_V   = 0.0_WP
    real(WP) :: SFull_W   = 0.0_WP
    real(WP) :: SFull_E   = 0.0_WP

    real(WP), allocatable, dimension(:) :: Xs_cs
    real(WP), allocatable, dimension(:) :: ffluc_cs
    real(WP), allocatable, dimension(:) :: XsArms_cs
!!$     real(WP), allocatable, dimension(:) :: XsFull_cs
    real(WP), allocatable, dimension(:) :: S_cs
    real(WP), allocatable, dimension(:) :: SArms_cs
    real(WP), allocatable, dimension(:) :: SFull_cs

    ! TC:
    real(WP) :: fluc_rho = 0.0_WP
    real(WP) :: ffluc_U  = 0.0_WP
    real(WP) :: ffluc_V  = 0.0_WP
    real(WP) :: ffluc_W  = 0.0_WP
    real(WP) :: ffluc_E  = 0.0_WP

    real(WP) :: mean_rho = 0.0_WP
    real(WP) :: fav_U    = 0.0_WP
    real(WP) :: fav_V    = 0.0_WP
    real(WP) :: fav_W    = 0.0_WP
    real(WP) :: fav_E    = 0.0_WP
    real(WP) :: mean_p   = 0.0_WP

    real(WP) :: field_rho = 0.0_WP
    real(WP) :: field_U   = 0.0_WP
    real(WP) :: field_V   = 0.0_WP
    real(WP) :: field_W   = 0.0_WP
    real(WP) :: field_E   = 0.0_WP

    real(WP) :: dmean_rho = 0.0_WP
    real(WP) :: dfav_U    = 0.0_WP
    real(WP) :: dfav_V    = 0.0_WP
    real(WP) :: dfav_W    = 0.0_WP
    real(WP) :: dfav_E    = 0.0_WP

    real(WP) :: rhoupup = 0.0_WP
    real(WP) :: rhovpvp = 0.0_WP
    real(WP) :: rhowpwp = 0.0_WP
    real(WP) :: rhoEpEp = 0.0_WP

    real(WP) :: drhoupup = 0.0_WP
    real(WP) :: drhovpvp = 0.0_WP
    real(WP) :: drhowpwp = 0.0_WP
    real(WP) :: drhoEpEp = 0.0_WP

    real(WP) :: Arms_rho = 0.0_WP
    real(WP) :: Arms_U   = 0.0_WP
    real(WP) :: Arms_V   = 0.0_WP
    real(WP) :: Arms_W   = 0.0_WP
    real(WP) :: Arms_E   = 0.0_WP

    real(WP) :: dArms_rho = 0.0_WP
    real(WP) :: dArms_U   = 0.0_WP
    real(WP) :: dArms_V   = 0.0_WP
    real(WP) :: dArms_W   = 0.0_WP
    real(WP) :: dArms_E   = 0.0_WP

    real(WP) :: grx_DA_rho  = 0.0_WP
    real(WP) :: grx_DA_U    = 0.0_WP
    real(WP) :: grx_DA_V    = 0.0_WP
    real(WP) :: grx_DA_W    = 0.0_WP
    real(WP) :: grx_DA_E    = 0.0_WP
    real(WP) :: grx_DA_p    = 0.0_WP

    real(WP) :: grt_DA_rms_rho = 0.0_WP
    real(WP) :: grt_DA_rms_U   = 0.0_WP
    real(WP) :: grt_DA_rms_V   = 0.0_WP
    real(WP) :: grt_DA_rms_W   = 0.0_WP
    real(WP) :: grt_DA_rms_E   = 0.0_WP
    real(WP) :: grt_DA_rms_p   = 0.0_WP

    real(WP) :: grx_DA_rms_rho = 0.0_WP
    real(WP) :: grx_DA_rms_U   = 0.0_WP
    real(WP) :: grx_DA_rms_V   = 0.0_WP
    real(WP) :: grx_DA_rms_W   = 0.0_WP
    real(WP) :: grx_DA_rms_E   = 0.0_WP
    real(WP) :: grx_DA_rms_p   = 0.0_WP

    real(WP) :: ygArms_rho = 0.0_WP
    real(WP) :: ygArms_U   = 0.0_WP
    real(WP) :: ygArms_V   = 0.0_WP
    real(WP) :: ygArms_W   = 0.0_WP
    real(WP) :: ygArms_E   = 0.0_WP

    real(WP), allocatable, dimension(:) ::      field_cs
    real(WP), allocatable, dimension(:) ::        fav_cs
    real(WP), allocatable, dimension(:) ::       dfav_cs
    real(WP), allocatable, dimension(:) ::       Arms_cs
    real(WP), allocatable, dimension(:) ::      dArms_cs
    real(WP), allocatable, dimension(:) ::     ygArms_cs
    real(WP), allocatable, dimension(:) ::     grx_DA_cs
    real(WP), allocatable, dimension(:) :: grt_DA_rms_cs
    real(WP), allocatable, dimension(:) :: grx_DA_rms_cs

    ! baseflow variables
    real(WP) :: base_rho   = 0.0_WP
    real(WP) :: base_U     = 0.0_WP
    real(WP) :: base_V     = 0.0_WP
    real(WP) :: base_W     = 0.0_WP
    real(WP) :: base_E     = 0.0_WP
    real(WP) :: base_p     = 0.0_WP

    real(WP) :: ddy_base_rho = 0.0_WP
    real(WP) :: ddy_base_U   = 0.0_WP
    real(WP) :: ddy_base_V   = 0.0_WP
    real(WP) :: ddy_base_W   = 0.0_WP
    real(WP) :: ddy_base_E   = 0.0_WP
    real(WP) :: ddy_base_p   = 0.0_WP

    real(WP) :: ddx_base_rho = 0.0_WP
    real(WP) :: ddx_base_U   = 0.0_WP
    real(WP) :: ddx_base_V   = 0.0_WP
    real(WP) :: ddx_base_W   = 0.0_WP
    real(WP) :: ddx_base_E   = 0.0_WP
    real(WP) :: ddx_base_p   = 0.0_WP

    real(WP) :: src_base_rho  = 0.0_WP
    real(WP) :: src_base_rhoU = 0.0_WP
    real(WP) :: src_base_rhoV = 0.0_WP
    real(WP) :: src_base_rhoW = 0.0_WP
    real(WP) :: src_base_rhoE = 0.0_WP

    real(WP) :: wall_base_u   = 0.0_WP

    real(WP), allocatable, dimension(:) ::     base_cs
    real(WP), allocatable, dimension(:) :: ddy_base_cs
    real(WP), allocatable, dimension(:) :: ddx_base_cs
    real(WP), allocatable, dimension(:) :: src_base_rhos

    ! RANS variables
    real(WP), allocatable, dimension(:) ::         Xs_tvar
    real(WP), allocatable, dimension(:) ::        fav_tvar
    real(WP), allocatable, dimension(:) ::       dfav_tvar
    real(WP), allocatable, dimension(:) ::     grx_DA_tvar
    real(WP), allocatable, dimension(:) :: grx_DA_rms_tvar
    real(WP), allocatable, dimension(:) :: grt_DA_rms_tvar
    real(WP), allocatable, dimension(:) ::          S_tvar

    integer(c_int) :: ip
    integer(c_int) :: ip_base    ! there's no baseflow info for tvars

  end type largo_BL_spatiotemporal_consistent_workspace_type

  procedure(prestep_setamean_rans), pointer :: largo_BL_spatiotemporal_consistent_prestep_sEtaMean_rans => NULL()
  procedure(sourcevec_rans),        pointer :: largo_BL_spatiotemporal_consistent_sEta_rans             => NULL()

  integer(c_int), parameter :: irho  = 1
  integer(c_int), parameter :: irhoU = 2
  integer(c_int), parameter :: irhoV = 3
  integer(c_int), parameter :: irhoW = 4
  integer(c_int), parameter :: irhoE = 5
!!$   real(WP), parameter :: eps = 1.0E-10_WP
  real(WP), parameter :: eps = 1.0E-7_WP

  ! Number of equations
  integer(c_int) :: neq_ = 0

  ! Number of species
  integer(c_int) :: ns_  = 0

  ! Number of tubulence variables
  integer(c_int) :: ntvar_  = 0

  ! Index of first tubulence variable
  integer(c_int) :: itvar0_  = 0

  public  :: largo_BL_spatiotemporal_consistent_allocate
  public  :: largo_BL_spatiotemporal_consistent_deallocate
  public  :: largo_BL_spatiotemporal_consistent_init
  public  :: largo_BL_spatiotemporal_consistent_preStep_sEta
  public  :: largo_BL_spatiotemporal_consistent_preStep_sEta_innery
  public  :: largo_BL_spatiotemporal_consistent_preStep_sEta_innerxz
  public  :: largo_BL_spatiotemporal_consistent_preStep_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_continuity_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_xMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_yMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_zMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_energy_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_ispecies_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_species_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_continuity_sEta_
  public  :: largo_BL_spatiotemporal_consistent_xMomentum_sEta_
  public  :: largo_BL_spatiotemporal_consistent_yMomentum_sEta_
  public  :: largo_BL_spatiotemporal_consistent_zMomentum_sEta_
  public  :: largo_BL_spatiotemporal_consistent_energy_sEta_
  public  :: largo_BL_spatiotemporal_consistent_ispecies_sEta_
  public  :: largo_BL_spatiotemporal_consistent_species_sEta_

  public  :: largo_BL_spatiotemporal_consistent_sEtaMean
  public  :: largo_BL_spatiotemporal_consistent_sEta

  private :: largo_BL_spatiotemporal_consistent_preStep_sEta_

  public  :: largo_BL_spatiotemporal_consistent_init_wall_baseflow
  public  :: largo_BL_spatiotemporal_consistent_preStep_baseflow

  public  :: largo_BL_spatiotemporal_consistent_get_ntvar_rans
  public  :: largo_BL_spatiotemporal_consistent_init_rans
  public  :: largo_BL_spatiotemporal_consistent_prestep_sEtaMean_rans
  public  :: largo_BL_spatiotemporal_consistent_sEta_rans

contains

  subroutine largo_BL_spatiotemporal_consistent_allocate(cp, neq, ns, ntvar, &
                                                         ransmodel)

    type(largo_workspace_ptr), intent(out) :: cp
    ! neq=number of equations, might be needed later
    integer(c_int), intent(in)   :: neq
    integer(c_int), intent(in)   :: ns    ! number of species
    integer(c_int), intent(in)   :: ntvar ! number of turbulence variables
    character(len=*)                                                 :: ransmodel
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp

    ! Allocate derived type variable
    allocate(auxp)

    ! Initialize values
    neq_    = neq
    ns_     = ns
    ntvar_  = ntvar
    itvar0_ = 5 + ns_ + 1

    ! Allocate arrays for species
    if (ns_ > 0) then
      allocate(auxp%Xs_cs         (1:ns_))
      allocate(auxp%ffluc_cs      (1:ns_))
      allocate(auxp%fav_cs        (1:ns_))
      allocate(auxp%dfav_cs       (1:ns_))
      allocate(auxp%field_cs      (1:ns_))
      allocate(auxp%XsArms_cs     (1:ns_))
      allocate(auxp%Arms_cs       (1:ns_))
      allocate(auxp%dArms_cs      (1:ns_))
      allocate(auxp%ygArms_cs     (1:ns_))

      allocate(auxp%S_cs          (1:ns_))
      allocate(auxp%SArms_cs      (1:ns_))
      allocate(auxp%SFull_cs      (1:ns_))

      allocate(auxp%grx_DA_cs     (1:ns_))
      allocate(auxp%grt_DA_rms_cs (1:ns_))
      allocate(auxp%grx_DA_rms_cs (1:ns_))

      allocate(auxp%base_cs       (1:ns_))
      allocate(auxp%ddy_base_cs   (1:ns_))
      allocate(auxp%ddx_base_cs   (1:ns_))
      allocate(auxp%src_base_rhos (1:ns_))

      auxp%Xs_cs         = 0.0_WP
      auxp%ffluc_cs      = 0.0_WP
      auxp%fav_cs        = 0.0_WP
      auxp%dfav_cs       = 0.0_WP
      auxp%field_cs      = 0.0_WP
      auxp%XsArms_cs     = 0.0_WP
      auxp%Arms_cs       = 0.0_WP
      auxp%dArms_cs      = 0.0_WP
      auxp%ygArms_cs     = 0.0_WP

      auxp%S_cs          = 0.0_WP
      auxp%SArms_cs      = 0.0_WP
      auxp%SFull_cs      = 0.0_WP

      auxp%grx_DA_cs     = 0.0_WP
      auxp%grt_DA_rms_cs = 0.0_WP
      auxp%grx_DA_rms_cs = 0.0_WP
      auxp%base_cs       = 0.0_WP
      auxp%ddy_base_cs   = 0.0_WP
      auxp%ddx_base_cs   = 0.0_WP
      auxp%src_base_rhos = 0.0_WP
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
    case ("self-similar")
      ! ntvar can have any value, all variables are
      ! computed assuming self-similar evolution
    case default
      ! FIXME: RANS model not defined
    end select

    ! Allocate arrays for turbulent variables
    if (ntvar_ > 0) then
      allocate(auxp%Xs_tvar         (1:ntvar_))
      allocate(auxp%fav_tvar        (1:ntvar_))
      allocate(auxp%dfav_tvar       (1:ntvar_))
      allocate(auxp%grx_DA_tvar     (1:ntvar_))
      allocate(auxp%grx_DA_rms_tvar (1:ntvar_))
      allocate(auxp%grt_DA_rms_tvar (1:ntvar_))
      allocate(auxp%S_tvar          (1:ntvar_))

      auxp%Xs_tvar         = 0.0_WP
      auxp%fav_tvar        = 0.0_WP
      auxp%dfav_tvar       = 0.0_WP
      auxp%grx_DA_tvar     = 0.0_WP
      auxp%grx_DA_rms_tvar = 0.0_WP
      auxp%grt_DA_rms_tvar = 0.0_WP
      auxp%S_tvar          = 0.0_WP
    end if

    ! Initialize turbulence variables function pointers
    largo_BL_spatiotemporal_consistent_prestep_sEtaMean_rans &
      & => largo_bl_spatiotemporal_consistent_prestep_sEtaMean_rans_gen
    largo_BL_spatiotemporal_consistent_sEta_rans             &
      & => largo_bl_spatiotemporal_consistent_sEta_rans_gen

    ! Pressure variable indices
    auxp%ip      = 5 + ns_ + ntvar_ + 1
    auxp%ip_base = 5 + ns_ + 1

    ! Get C pointer from Fortran pointer
    cp = c_loc(auxp)

  end subroutine largo_BL_spatiotemporal_consistent_allocate


  subroutine largo_BL_spatiotemporal_consistent_deallocate(cp)

    type(largo_workspace_ptr), intent(inout)  :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp

    call c_f_pointer(cp, auxp)

    ! Deallocate arrays for species
    if (allocated(auxp%Xs_cs             )) deallocate(auxp%Xs_cs            )
    if (allocated(auxp%ffluc_cs          )) deallocate(auxp%ffluc_cs         )
    if (allocated(auxp%fav_cs            )) deallocate(auxp%fav_cs           )
    if (allocated(auxp%dfav_cs           )) deallocate(auxp%dfav_cs          )
    if (allocated(auxp%field_cs          )) deallocate(auxp%field_cs         )
    if (allocated(auxp%XsArms_cs         )) deallocate(auxp%XsArms_cs        )
    if (allocated(auxp%Arms_cs           )) deallocate(auxp%Arms_cs          )
    if (allocated(auxp%dArms_cs          )) deallocate(auxp%dArms_cs         )
    if (allocated(auxp%ygArms_cs         )) deallocate(auxp%ygArms_cs        )
    if (allocated(auxp%S_cs              )) deallocate(auxp%S_cs             )
    if (allocated(auxp%SFull_cs          )) deallocate(auxp%SFull_cs         )

    if (allocated(auxp%grx_DA_cs         )) deallocate(auxp%grx_DA_cs        )
    if (allocated(auxp%grt_DA_rms_cs     )) deallocate(auxp%grt_DA_rms_cs    )
    if (allocated(auxp%grx_DA_rms_cs     )) deallocate(auxp%grx_DA_rms_cs    )

    if (allocated(auxp%base_cs           )) deallocate(auxp%base_cs          )
    if (allocated(auxp%ddy_base_cs       )) deallocate(auxp%ddy_base_cs      )
    if (allocated(auxp%ddx_base_cs       )) deallocate(auxp%ddx_base_cs      )
    if (allocated(auxp%src_base_rhos     )) deallocate(auxp%src_base_rhos    )

    ! Deallocate arrays for RANS turbulence variables
    if (allocated(auxp%Xs_tvar           ))  deallocate(auxp%Xs_tvar         )
    if (allocated(auxp%fav_tvar          ))  deallocate(auxp%fav_tvar        )
    if (allocated(auxp%dfav_tvar         ))  deallocate(auxp%dfav_tvar       )
    if (allocated(auxp%grx_DA_tvar       ))  deallocate(auxp%grx_DA_tvar     )
    if (allocated(auxp%grx_DA_rms_tvar   ))  deallocate(auxp%grx_DA_rms_tvar )
    if (allocated(auxp%grt_DA_rms_tvar   ))  deallocate(auxp%grt_DA_rms_tvar )
    if (allocated(auxp%S_tvar            ))  deallocate(auxp%S_tvar          )

    ! Deallocate array of derived types
    deallocate(auxp)

    ! Nullify C pointer
    cp = c_null_ptr

  end subroutine largo_BL_spatiotemporal_consistent_deallocate


  ! get number of turbulence variables
  subroutine largo_BL_spatiotemporal_consistent_get_ntvar_rans(cp, ntvar)

    ! largo workspace C pointer
    type(largo_workspace_ptr), intent(in)            :: cp
    integer(c_int)           , intent(out)           :: ntvar

!!$     ! Get Fortran pointer from C pointer
!!$     call c_f_pointer(cp, auxp)

    ntvar = ntvar_

  end subroutine largo_BL_spatiotemporal_consistent_get_ntvar_rans


  subroutine largo_BL_spatiotemporal_consistent_init(cp, grx_delta, grx_DA, grx_DA_rms)

    real(WP), intent(in)                  :: grx_delta
    real(WP), dimension(*), intent(in)    :: grx_DA
    real(WP), dimension(*), intent(in)    :: grx_DA_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Set growth rates
    auxp%grx_delta   = grx_delta

    auxp%grx_DA_rho  = grx_DA(irho )
    auxp%grx_DA_U    = grx_DA(irhoU)
    auxp%grx_DA_V    = grx_DA(irhoV)
    auxp%grx_DA_W    = grx_DA(irhoW)
    auxp%grx_DA_E    = grx_DA(irhoE)
    auxp%grx_DA_p    = grx_DA(auxp%ip)

    auxp%grx_DA_rms_rho  = grx_DA_rms(irho )
    auxp%grx_DA_rms_U    = grx_DA_rms(irhoU)
    auxp%grx_DA_rms_V    = grx_DA_rms(irhoV)
    auxp%grx_DA_rms_W    = grx_DA_rms(irhoW)
    auxp%grx_DA_rms_E    = grx_DA_rms(irhoE)
    auxp%grx_DA_rms_p    = grx_DA_rms(auxp%ip)

    do is=1, ns_
      auxp%grx_DA_cs    (is) = grx_DA    (5+is)
      auxp%grx_DA_rms_cs(is) = grx_DA_rms(5+is)
    end do

    ! Compute grt_delta using as velocity scale the 
    ! inviscid streamwise velocity at the wall
    ! NOTE: added the computation of grt_delta here as 
    ! well to avoid having as a requirement to call this method
    ! before the init_wall_baseflow one
    if (auxp%wall_base_u /= 0.0_WP) then
      auxp%grt_delta      = auxp%grx_delta      * auxp%wall_base_u 
      auxp%grt_DA_rms_rho = grx_DA_rms(irho )   * auxp%wall_base_u  
      auxp%grt_DA_rms_U   = grx_DA_rms(irhoU)   * auxp%wall_base_u
      auxp%grt_DA_rms_V   = grx_DA_rms(irhoV)   * auxp%wall_base_u
      auxp%grt_DA_rms_W   = grx_DA_rms(irhoW)   * auxp%wall_base_u
      auxp%grt_DA_rms_E   = grx_DA_rms(irhoE)   * auxp%wall_base_u
      auxp%grt_DA_rms_p   = grx_DA_rms(auxp%ip) * auxp%wall_base_u
      do is=1, ns_
        auxp%grt_DA_rms_cs(is) = grx_DA_rms(5+is) * auxp%wall_base_u
      end do
    end if

  end subroutine largo_BL_spatiotemporal_consistent_init


  ! RANS initialization
  subroutine largo_BL_spatiotemporal_consistent_init_rans(cp, grx_delta, grx_DA, grx_DA_rms)

    real(WP), intent(in)                             :: grx_delta
    real(WP), dimension(*), intent(in)               :: grx_DA
    real(WP), dimension(*), intent(in)               :: grx_DA_rms
    type(largo_workspace_ptr), intent(in)            :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer  :: auxp
    integer(c_int) :: it

    call largo_BL_spatiotemporal_consistent_init(cp, grx_delta, grx_DA, grx_DA_rms)

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    do it=1, ntvar_
      auxp%grx_DA_tvar     (it) = grx_DA     (itvar0_ + it - 1)
      auxp%grx_DA_rms_tvar (it) = grx_DA_rms (itvar0_ + it - 1)
    end do

    if (auxp%wall_base_u /= 0.0_WP) then
      do it=1, ntvar_
        auxp%grt_DA_rms_tvar(it) = grx_DA_rms(itvar0_+it-1) * auxp%wall_base_u
      end do
    end if

  end subroutine largo_BL_spatiotemporal_consistent_init_rans


  subroutine largo_BL_spatiotemporal_consistent_init_wall_baseflow(cp,  &   
               wall_base, wall_ddy_base, wall_ddt_base,      &
                          wall_ddx_base, wall_src_base)

    real(WP), dimension(*), intent(in)   :: wall_base
    real(WP), dimension(*), intent(in)   :: wall_ddy_base
    real(WP), dimension(*), intent(in)   :: wall_ddt_base
    real(WP), dimension(*), intent(in)   :: wall_ddx_base
    real(WP), dimension(*), intent(in)   :: wall_src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store relevant wall baseflow information
    auxp%wall_base_u  = wall_base(irhoU) / wall_base(irho)

    ! Compute grt_* using as velocity scale the 
    ! inviscid streamwise velocity at the wall 
    auxp%grt_delta       = auxp%grx_delta      * auxp%wall_base_u 
    auxp%grt_DA_rms_rho  = auxp%grx_DA_rms_rho * auxp%wall_base_u  
    auxp%grt_DA_rms_U    = auxp%grx_DA_rms_U   * auxp%wall_base_u
    auxp%grt_DA_rms_V    = auxp%grx_DA_rms_V   * auxp%wall_base_u
    auxp%grt_DA_rms_W    = auxp%grx_DA_rms_W   * auxp%wall_base_u
    auxp%grt_DA_rms_E    = auxp%grx_DA_rms_E   * auxp%wall_base_u
    auxp%grt_DA_rms_p    = auxp%grx_DA_rms_p   * auxp%wall_base_u
    do is=1, ns_
      auxp%grt_DA_rms_cs(is) = auxp%grx_DA_rms_cs(is) * auxp%wall_base_u
    end do

  end subroutine largo_BL_spatiotemporal_consistent_init_wall_baseflow


  subroutine largo_BL_spatiotemporal_consistent_preStep_baseflow(cp, base, &
                              ddy_base, ddt_base, ddx_base, src_base)

    real(WP), dimension(*), intent(in)    :: base
    real(WP), dimension(*), intent(in)    :: ddy_base
    real(WP), dimension(*), intent(in)    :: ddt_base
    real(WP), dimension(*), intent(in)    :: ddx_base
    real(WP), dimension(*), intent(in)    :: src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store baseflow information
    if (base(irho) > 0.0_WP) then
      auxp%base_rho = base(irho )
      auxp%base_U   = base(irhoU)/base(irho )
      auxp%base_V   = base(irhoV)/base(irho )
      auxp%base_W   = base(irhoW)/base(irho )
      auxp%base_E   = base(irhoE)/base(irho )
      auxp%base_p   = base(auxp%ip_base)

      auxp%ddy_base_rho = ddy_base(irho )
      auxp%ddy_base_U   = ddy_base(irhoU)/base(irho ) - auxp%base_U/base(irho ) * ddy_base(irho )
      auxp%ddy_base_V   = ddy_base(irhoV)/base(irho ) - auxp%base_V/base(irho ) * ddy_base(irho )
      auxp%ddy_base_W   = ddy_base(irhoW)/base(irho ) - auxp%base_W/base(irho ) * ddy_base(irho )
      auxp%ddy_base_E   = ddy_base(irhoE)/base(irho ) - auxp%base_E/base(irho ) * ddy_base(irho )
      auxp%ddy_base_p   = ddy_base(auxp%ip_base)

      auxp%ddx_base_rho = ddx_base(irho )
      auxp%ddx_base_U   = ddx_base(irhoU)/base(irho ) - auxp%base_U/base(irho ) * ddx_base(irho )
      auxp%ddx_base_V   = ddx_base(irhoV)/base(irho ) - auxp%base_V/base(irho ) * ddx_base(irho )
      auxp%ddx_base_W   = ddx_base(irhoW)/base(irho ) - auxp%base_W/base(irho ) * ddx_base(irho )
      auxp%ddx_base_E   = ddx_base(irhoE)/base(irho ) - auxp%base_E/base(irho ) * ddx_base(irho )
      auxp%ddx_base_p   = ddx_base(auxp%ip_base)
    end if

    auxp%src_base_rho  = src_base(irho )
    auxp%src_base_rhoU = src_base(irhoU)
    auxp%src_base_rhoV = src_base(irhoV)
    auxp%src_base_rhoW = src_base(irhoW)
    auxp%src_base_rhoE = src_base(irhoE)

    do is=1, ns_
      auxp%base_cs      (is) =     base(5+is)/base(irho )
      auxp%ddy_base_cs  (is) = ddy_base(5+is)/base(irho ) - auxp%base_cs(is)/base(irho ) * ddy_base(irho )
      auxp%ddx_base_cs  (is) = ddx_base(5+is)/base(irho ) - auxp%base_cs(is)/base(irho ) * ddx_base(irho )
      auxp%src_base_rhos(is) = src_base(5+is)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_preStep_baseflow


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: ddy_mean
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Compute/get mean and Favre averages
    auxp%mean_rho = mean(irho)
    auxp%fav_U    = mean(irhoU)/mean(irho)
    auxp%fav_V    = mean(irhoV)/mean(irho)
    auxp%fav_W    = mean(irhoW)/mean(irho)
    auxp%fav_E    = mean(irhoE)/mean(irho)
    do is=1, ns_
      auxp%fav_cs(is) = mean(5+is)/mean(irho)
    end do
    auxp%mean_p    = mean(auxp%ip)

    ! Compute/get y-derivative of Favre averages
    auxp%dmean_rho = ddy_mean(irho)
    auxp%dfav_U    = ddy_mean(irhoU)/mean(irho) &
      &             - auxp%fav_U/mean(irho) * ddy_mean(irho)
    auxp%dfav_V    = ddy_mean(irhoV)/mean(irho) &
      &             - auxp%fav_V/mean(irho) * ddy_mean(irho)
    auxp%dfav_W    = ddy_mean(irhoW)/mean(irho) &
      &             - auxp%fav_W/mean(irho) * ddy_mean(irho)
    auxp%dfav_E    = ddy_mean(irhoE)/mean(irho) &
      &             - auxp%fav_E/mean(irho) * ddy_mean(irho)

    do is=1, ns_
      auxp%dfav_cs(is)  =  ddy_mean(5+is)/mean(irho) &
      &             - auxp%fav_cs(is)/mean(irho) * ddy_mean(irho)
    end do

    ! These ones depend on y only
    auxp%Xs_rho =  - auxp%ddx_base_rho - auxp%grx_DA_rho * (auxp%mean_rho-auxp%base_rho) + y * auxp%grx_delta * (ddy_mean(irho )  -auxp%ddy_base_rho)
    auxp%Xs_U   =  - auxp%ddx_base_U   - auxp%grx_DA_U   * (auxp%fav_U   -auxp%base_U  ) + y * auxp%grx_delta * (auxp%dfav_U      -auxp%ddy_base_U  )
    auxp%Xs_V   =  - auxp%ddx_base_V   - auxp%grx_DA_V   * (auxp%fav_V   -auxp%base_V  ) + y * auxp%grx_delta * (auxp%dfav_V      -auxp%ddy_base_V  )
    auxp%Xs_W   =  - auxp%ddx_base_W   - auxp%grx_DA_W   * (auxp%fav_W   -auxp%base_W  ) + y * auxp%grx_delta * (auxp%dfav_W      -auxp%ddy_base_W  )
    auxp%Xs_E   =  - auxp%ddx_base_E   - auxp%grx_DA_E   * (auxp%fav_E   -auxp%base_E  ) + y * auxp%grx_delta * (auxp%dfav_E      -auxp%ddy_base_E  )

    auxp%Xs_p    = - auxp%ddx_base_p   - auxp%grx_DA_p   * (mean(auxp%ip)-auxp%base_p )  + y * auxp%grx_delta * (ddy_mean(auxp%ip)-auxp%ddy_base_p  )

    do is=1, ns_
      auxp%Xs_cs(is)  = - auxp%ddx_base_cs(is) - auxp%grx_DA_cs(is) * (auxp%fav_cs(is)-auxp%base_cs(is)) + y * auxp%grx_delta * (auxp%dfav_cs(is)-auxp%ddy_base_cs(is)) 
    end do

    ! Mean sources for primitive variables
    auxp%S_rho =  auxp%fav_U * auxp%Xs_rho +             auxp%mean_rho * auxp%Xs_U
    auxp%S_U   =  auxp%fav_U * auxp%Xs_U   +                                             1.0_WP/auxp%mean_rho * auxp%Xs_P
    auxp%S_V   =  auxp%fav_U * auxp%Xs_V
    auxp%S_W   =  auxp%fav_U * auxp%Xs_W
    auxp%S_E   =  auxp%fav_U * auxp%Xs_E   + auxp%mean_p/auxp%mean_rho * auxp%Xs_U + auxp%fav_u/auxp%mean_rho * auxp%Xs_P

    do is=1, ns_
      auxp%S_cs(is)  = auxp%fav_U * auxp%Xs_cs(is)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEtaMean


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: dmean_rqq
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in)          :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Compute \mean{ru"u"} for each u component and for energy
    ! In a continuous setting, these must be nonnegative so
    ! enforce that constraint discretely to increase robustness.
    auxp%rhoupup  = max(0.0_WP, mean_rqq(irhoU) - auxp%mean_rho * auxp%fav_u**2)

    auxp%rhovpvp  = max(0.0_WP, mean_rqq(irhoV) - auxp%mean_rho * auxp%fav_v**2)

    auxp%rhowpwp  = max(0.0_WP, mean_rqq(irhoW) - auxp%mean_rho * auxp%fav_w**2)

    auxp%rhoEpEp  = max(0.0_WP, mean_rqq(irhoE) - auxp%mean_rho * auxp%fav_E**2)

    ! Compute d(\mean{ru"u"})/dy for each u component and for enegry
    auxp%drhoupup = dmean_rqq(irhoU) &
      &          - auxp%dmean_rho * auxp%fav_u**2 &
      &          - auxp%mean_rho  * 2.0_WP*auxp%fav_u*auxp%dfav_u

    auxp%drhovpvp = dmean_rqq(irhoV) &
      &          - auxp%dmean_rho * auxp%fav_v**2 &
      &          - auxp%mean_rho  * 2.0_WP*auxp%fav_v*auxp%dfav_v

    auxp%drhowpwp = dmean_rqq(irhoW) &
      &          - auxp%dmean_rho * auxp%fav_w**2 &
      &          - auxp%mean_rho  * 2.0_WP*auxp%fav_w*auxp%dfav_w

    auxp%drhoEpEp = dmean_rqq(irhoE) &
      &          - auxp%dmean_rho * auxp%fav_E**2 &
      &          - auxp%mean_rho  * 2.0_WP*auxp%fav_E*auxp%dfav_E

    ! Assign Arms_rho
    !NOTE: Arms terms for density are based on that for the mean
    !      no need to compute them
    !auxp%Arms_rho = auxp%mean_rho

    ! Compute Arms_U=sqrt{2*k}
    auxp%Arms_U   = sqrt((auxp%rhoupup + auxp%rhovpvp + auxp%rhowpwp) / &
     &                     auxp%mean_rho)
    auxp%Arms_V   = auxp%Arms_U
    auxp%Arms_W   = auxp%Arms_U
    auxp%Arms_E   = sqrt(auxp%rhoEpEp/ auxp%mean_rho)

    ! Assign dArms_rho
    !NOTE: Arms terms for density are based on that for the mean
    !      no need to compute them
    !auxp%dArms_rho = auxp%dmean_rho

    ! Compute d(Armsu)/dy=d(sqrt{2*tke})/dy
    auxp%dArms_U  = 0.0_WP
    auxp%dArms_V  = 0.0_WP
    auxp%dArms_W  = 0.0_WP

    if (auxp%Arms_U > eps) then
      auxp%dArms_U  = 0.5_WP / auxp%Arms_U * ( &
        &   (auxp%drhoupup + auxp%drhovpvp + auxp%drhowpwp) / auxp%mean_rho &
        & - (auxp%rhoupup + auxp%rhovpvp + auxp%rhowpwp) / auxp%mean_rho**2 * auxp%dmean_rho )
      auxp%dArms_V  = auxp%dArms_U
      auxp%dArms_W  = auxp%dArms_U
    end if

    auxp%dArms_E  = 0.0_WP
    if (auxp%Arms_E > eps) &
      & auxp%dArms_E  = 0.5_WP / auxp%Arms_E * ( &
          & auxp%drhoEpEp / auxp%mean_rho - auxp%rhoEpEp / auxp%mean_rho**2 * auxp%dmean_rho )

    do is=1, ns_
      auxp%Arms_cs(is)  = auxp%mean_rho
      auxp%dArms_cs(is) = auxp%dmean_rho
    end do

    ! These ones depend on y only
    !NOTE: Arms terms for density are based on that for the mean
    !      no need to compute them
    !auxp%ygArms_rho  = 0.0_WP
    !if (auxp%Arms_rho > eps) auxp%ygArms_rho = y * auxp%gr_delta * auxp%dArms_rho/auxp%Arms_rho

    ! FIXME: Fix growth rates
    auxp%ygArms_U = 0.0_WP
    if (auxp%Arms_U   > eps) auxp%ygArms_U   = y * auxp%grt_delta * auxp%dArms_U/auxp%Arms_U

    auxp%ygArms_V = 0.0_WP
    if (auxp%Arms_V   > eps) auxp%ygArms_V   = y * auxp%grt_delta * auxp%dArms_V/auxp%Arms_V

    auxp%ygArms_W = 0.0_WP
    if (auxp%Arms_W   > eps) auxp%ygArms_W   = y * auxp%grt_delta * auxp%dArms_W/auxp%Arms_W

    auxp%ygArms_E = 0.0_WP
    if (auxp%Arms_E   > eps) auxp%ygArms_E   = y * auxp%grt_delta * auxp%dArms_E/auxp%Arms_E

    do is=1, ns_
      auxp%ygArms_cs(is) = 0.0_WP
      if (auxp%Arms_cs(is) > eps)  auxp%ygArms_cs(is) = y * auxp%grt_delta * auxp%dArms_cs(is)/auxp%Arms_cs(is)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_innerxz(cp, qflow)

    real(WP), dimension(*), intent(in) :: qflow
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    auxp%field_rho = qflow(irho)
    auxp%field_U   = qflow(irhoU)/qflow(irho)
    auxp%field_V   = qflow(irhoV)/qflow(irho)
    auxp%field_W   = qflow(irhoW)/qflow(irho)
    auxp%field_E   = qflow(irhoE)/qflow(irho)

    auxp%fluc_rho  = auxp%field_rho - auxp%mean_rho
    auxp%ffluc_U   = auxp%field_U   - auxp%fav_U
    auxp%ffluc_V   = auxp%field_V   - auxp%fav_V
    auxp%ffluc_W   = auxp%field_W   - auxp%fav_W
    auxp%ffluc_E   = auxp%field_E   - auxp%fav_E

    auxp%XsArms_rho = auxp%fluc_rho / auxp%mean_rho * auxp%Xs_rho 
    auxp%XsArms_U   = auxp%ffluc_U  * (- auxp%grt_DA_rms_U   + auxp%ygArms_U   )
    auxp%XsArms_V   = auxp%ffluc_V  * (- auxp%grt_DA_rms_V   + auxp%ygArms_V   )
    auxp%XsArms_W   = auxp%ffluc_W  * (- auxp%grt_DA_rms_W   + auxp%ygArms_W   )
    auxp%XsArms_E   = auxp%ffluc_E  * (- auxp%grt_DA_rms_E   + auxp%ygArms_E   )

    auxp%SArms_rho  = auxp%fav_U * auxp%XsArms_rho + auxp%fluc_rho * auxp%Xs_U
    auxp%SArms_U    = auxp%XsArms_U   
    auxp%SArms_V    = auxp%XsArms_V   
    auxp%SArms_W    = auxp%XsArms_W   
    auxp%SArms_E    = auxp%XsArms_E   

    do is=1, ns_
      auxp%field_cs (is)  = qflow(5+is)/qflow(irho)
      auxp%ffluc_cs (is)  = auxp%field_cs(is) - auxp%fav_cs   (is)
      auxp%XsArms_cs(is)  = auxp%ffluc_cs(is) * (- auxp%grt_DA_rms_cs(is) + auxp%ygArms_cs(is))
      auxp%SArms_cs (is)  = auxp%XsArms_cs(is)
    end do

    ! Compute mean plus fluctuations primitive sources
    auxp%SFull_rho = auxp%S_rho + auxp%SArms_rho
    auxp%SFull_U   = auxp%S_U   + auxp%SArms_U
    auxp%SFull_V   = auxp%S_V   + auxp%SArms_V
    auxp%SFull_W   = auxp%S_W   + auxp%SArms_W
    auxp%SFull_E   = auxp%S_E   + auxp%SArms_E

    do is=1, ns_
      auxp%SFull_cs(is) = auxp%S_cs(is) + auxp%SArms_cs(is) 
    end do

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_innerxz


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_innery(cp, y, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: ddy_mean
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: dmean_rqq
    type(largo_workspace_ptr), intent(in)        :: cp

    call largo_BL_spatiotemporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_spatiotemporal_consistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEta_innery


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEta(cp, y, qflow, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in) :: qflow
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: ddy_mean
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: dmean_rqq
    type(largo_workspace_ptr), intent(in)        :: cp

    ! Delegate to Y then XZ routines to ensure their behavior reproduced
    call largo_BL_spatiotemporal_consistent_preStep_sEta_innery(cp, y, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)
    call largo_BL_spatiotemporal_consistent_preStep_sEta_innerxz(cp, qflow)

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEta


  subroutine largo_BL_spatiotemporal_consistent_preStep_sEtaMean_rans_gen(cp, y, mean, ddy_mean)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: ddy_mean
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer   :: auxp
    integer(c_int) :: it

    call largo_BL_spatiotemporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! 
    do it=1, ntvar_
      auxp%fav_tvar (it) = mean(itvar0_+it-1)/auxp%mean_rho
      auxp%dfav_tvar(it) = ddy_mean(itvar0_+it-1)/auxp%mean_rho &
      &                     - auxp%fav_tvar(it)/auxp%mean_rho * auxp%dmean_rho
      auxp%Xs_tvar  (it) = -0.0_WP - auxp%grx_DA_tvar(it) * (auxp%fav_tvar(it)-0.0_WP) + y * auxp%grx_delta * (auxp%dfav_tvar(it)-0.0_WP) 
      auxp%S_tvar   (it) = auxp%fav_U * auxp%Xs_tvar(it)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_preStep_sEtaMean_rans_gen


#define DECLARE_SUBROUTINE(token)token (cp, A, B, src);\
  type(largo_workspace_ptr), intent(in)   :: cp;\
  real(WP)       , intent(in)             :: A, B;\
  real(WP)       , intent(inout)          :: src;\
  type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp;\
  call c_f_pointer(cp, auxp)\

  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_continuity_sEtaMean)
    src = A * src + B * (auxp%S_rho + auxp%src_base_rho)
  end subroutine largo_BL_spatiotemporal_consistent_continuity_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_xMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%S_U &
      &                  + auxp%fav_U    * auxp%S_rho + auxp%src_base_rhoU)
  end subroutine largo_BL_spatiotemporal_consistent_xMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_yMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%S_V &
      &                  + auxp%fav_V    * auxp%S_rho + auxp%src_base_rhoV)
  end subroutine largo_BL_spatiotemporal_consistent_yMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_zMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%S_W &
      &                  + auxp%fav_W    * auxp%S_rho + auxp%src_base_rhoW)
  end subroutine largo_BL_spatiotemporal_consistent_zMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_energy_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%S_E &
      &                  + auxp%fav_E    * auxp%S_rho + auxp%src_base_rhoE)
  end subroutine largo_BL_spatiotemporal_consistent_energy_sEtaMean


  subroutine largo_BL_spatiotemporal_consistent_ispecies_sEtaMean (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * (  auxp%mean_rho   * auxp%S_cs(is) &
      &                  + auxp%fav_cs(is) * auxp%S_rho  + auxp%src_base_rhos(is))
  end subroutine largo_BL_spatiotemporal_consistent_ispecies_sEtaMean


  subroutine largo_BL_spatiotemporal_consistent_species_sEtaMean (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_spatiotemporal_consistent_ispecies_sEtaMean (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_spatiotemporal_consistent_species_sEtaMean


  subroutine largo_BL_spatiotemporal_consistent_sEtaMean (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_spatiotemporal_consistent_continuity_sEtaMean (cp,      A, B, srcvec(irho ))
    call largo_BL_spatiotemporal_consistent_xMomentum_sEtaMean  (cp,      A, B, srcvec(irhou))
    call largo_BL_spatiotemporal_consistent_yMomentum_sEtaMean  (cp,      A, B, srcvec(irhov))
    call largo_BL_spatiotemporal_consistent_zMomentum_sEtaMean  (cp,      A, B, srcvec(irhow))
    call largo_BL_spatiotemporal_consistent_energy_sEtaMean     (cp,      A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_spatiotemporal_consistent_ispecies_sEtaMean (cp,      A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_sEtaMean


  subroutine largo_BL_spatiotemporal_consistent_sEta_rans_gen (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer       :: auxp
    integer(c_int)                            :: it

    call largo_BL_spatiotemporal_consistent_sEtaMean (cp, A, B, srcvec) 

    call c_f_pointer(cp, auxp)
    do it = 1, ntvar_
       srcvec(itvar0_+it-1) = A * srcvec(itvar0_+it-1) + B * (auxp%mean_rho     * auxp%S_tvar(it) &
      &                                                     + auxp%fav_tvar(it) * auxp%S_rho)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_sEta_rans_gen


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_continuity_sEta_)
    src = A * src + B * (auxp%SFull_rho + auxp%src_base_rho)
  end subroutine largo_BL_spatiotemporal_consistent_continuity_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_xMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%SFull_U &
      &                  + auxp%field_U   * auxp%SFull_rho + auxp%src_base_rhoU)
  end subroutine largo_BL_spatiotemporal_consistent_xMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_yMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%SFull_V &
      &                  + auxp%field_V   * auxp%SFull_rho + auxp%src_base_rhoV)
  end subroutine largo_BL_spatiotemporal_consistent_yMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_zMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%SFull_W &
      &                  + auxp%field_W   * auxp%SFull_rho + auxp%src_base_rhoW)
  end subroutine largo_BL_spatiotemporal_consistent_zMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_consistent_energy_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%SFull_E &
      &                  + auxp%field_E   * auxp%SFull_rho + auxp%src_base_rhoE)
  end subroutine largo_BL_spatiotemporal_consistent_energy_sEta_


  subroutine largo_BL_spatiotemporal_consistent_ispecies_sEta_ (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_spatiotemporal_consistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * (  auxp%field_rho    * auxp%SFull_cs(is) &
      &                  + auxp%field_cs(is) * auxp%SFull_rho + auxp%src_base_rhos(is))
  end subroutine largo_BL_spatiotemporal_consistent_ispecies_sEta_


  subroutine largo_BL_spatiotemporal_consistent_species_sEta_ (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)   :: cp
    real(WP)       , intent(in)             :: A, B
    real(WP), dimension(*), intent(inout)   :: srcvec ! "*" = ns_
    integer(c_int)                          :: is

    do is = 1, ns_
      call largo_BL_spatiotemporal_consistent_ispecies_sEta_ (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_spatiotemporal_consistent_species_sEta_


  subroutine largo_BL_spatiotemporal_consistent_sEta (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_spatiotemporal_consistent_continuity_sEta_ (cp, A, B, srcvec(irho ))
    call largo_BL_spatiotemporal_consistent_xMomentum_sEta_  (cp, A, B, srcvec(irhou))
    call largo_BL_spatiotemporal_consistent_yMomentum_sEta_  (cp, A, B, srcvec(irhov))
    call largo_BL_spatiotemporal_consistent_zMomentum_sEta_  (cp, A, B, srcvec(irhow))
    call largo_BL_spatiotemporal_consistent_energy_sEta_     (cp, A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_spatiotemporal_consistent_ispecies_sEta_ (cp, A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_spatiotemporal_consistent_sEta

end module largo_BL_spatiotemporal_consistent
