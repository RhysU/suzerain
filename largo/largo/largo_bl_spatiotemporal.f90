module largo_BL_spatiotemporal

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

  ! largo workspace type declaration
  type :: largo_BL_spatiotemporal_workspace_type

    real(WP) :: grt_delta   = 1.0_WP
    real(WP) :: grx_delta   = 1.0_WP

    real(WP) :: Xs_rho    = 0.0_WP
    real(WP) :: Xs_rhoU   = 0.0_WP
    real(WP) :: Xs_rhoV   = 0.0_WP
    real(WP) :: Xs_rhoW   = 0.0_WP
    real(WP) :: Xs_rhoE   = 0.0_WP

    real(WP) :: Xs_p      = 0.0_WP
    real(WP) :: Xs_rhoH   = 0.0_WP

    real(WP) :: ygrms_rho   = 0.0_WP
    real(WP) :: ygrms_rhoU  = 0.0_WP
    real(WP) :: ygrms_rhoV  = 0.0_WP
    real(WP) :: ygrms_rhoW  = 0.0_WP
    real(WP) :: ygrms_rhoE  = 0.0_WP

    real(WP) :: mean_rho   = 0.0_WP
    real(WP) :: mean_rhoU  = 0.0_WP
    real(WP) :: mean_rhoV  = 0.0_WP
    real(WP) :: mean_rhoW  = 0.0_WP
    real(WP) :: mean_rhoE  = 0.0_WP

    real(WP) :: fav_U      = 0.0_WP
    real(WP) :: fav_V      = 0.0_WP
    real(WP) :: fav_W      = 0.0_WP
    real(WP) :: fav_E      = 0.0_WP

    real(WP) :: mean_p     = 0.0_WP
    real(WP) :: fav_H      = 0.0_WP

    real(WP) :: fluc_rho   = 0.0_WP
    real(WP) :: fluc_rhoU  = 0.0_WP
    real(WP) :: fluc_rhoV  = 0.0_WP
    real(WP) :: fluc_rhoW  = 0.0_WP
    real(WP) :: fluc_rhoE  = 0.0_WP

    real(WP) :: dtsRms_rho   = 0.0_WP
    real(WP) :: dtsRms_rhoU  = 0.0_WP
    real(WP) :: dtsRms_rhoV  = 0.0_WP
    real(WP) :: dtsRms_rhoW  = 0.0_WP
    real(WP) :: dtsRms_rhoE  = 0.0_WP

    real(WP), allocatable, dimension(:) :: Xs_rhos
    real(WP), allocatable, dimension(:) :: ygrms_rhos
    real(WP), allocatable, dimension(:) :: mean_rhos
    real(WP), allocatable, dimension(:) :: fav_cs
    real(WP), allocatable, dimension(:) :: fluc_rhos
    real(WP), allocatable, dimension(:) :: dtsRms_rhos

    ! baseflow variables
    real(WP) :: base_rho      = 0.0_WP
    real(WP) :: base_rhoU     = 0.0_WP
    real(WP) :: base_rhoV     = 0.0_WP
    real(WP) :: base_rhoW     = 0.0_WP
    real(WP) :: base_rhoE     = 0.0_WP
    real(WP) :: base_p        = 0.0_WP

    real(WP) :: ddy_base_rho  = 0.0_WP
    real(WP) :: ddy_base_rhoU = 0.0_WP
    real(WP) :: ddy_base_rhoV = 0.0_WP
    real(WP) :: ddy_base_rhoW = 0.0_WP
    real(WP) :: ddy_base_rhoE = 0.0_WP
    real(WP) :: ddy_base_p    = 0.0_WP

    real(WP) :: ddx_base_rho  = 0.0_WP
    real(WP) :: ddx_base_rhoU = 0.0_WP
    real(WP) :: ddx_base_rhoV = 0.0_WP
    real(WP) :: ddx_base_rhoW = 0.0_WP
    real(WP) :: ddx_base_rhoE = 0.0_WP
    real(WP) :: ddx_base_p    = 0.0_WP

    real(WP) :: src_base_rho  = 0.0_WP
    real(WP) :: src_base_rhoU = 0.0_WP
    real(WP) :: src_base_rhoV = 0.0_WP
    real(WP) :: src_base_rhoW = 0.0_WP
    real(WP) :: src_base_rhoE = 0.0_WP

    real(WP) :: wall_base_u   = 0.0_WP

    real(WP) :: gr_DA_rho     = 0.0_WP
    real(WP) :: gr_DA_rhoU    = 0.0_WP
    real(WP) :: gr_DA_rhoV    = 0.0_WP
    real(WP) :: gr_DA_rhoW    = 0.0_WP
    real(WP) :: gr_DA_rhoE    = 0.0_WP
    real(WP) :: gr_DA_p       = 0.0_WP

    real(WP) :: gr_DA_rms_rho  = 0.0_WP
    real(WP) :: gr_DA_rms_rhoU = 0.0_WP
    real(WP) :: gr_DA_rms_rhoV = 0.0_WP
    real(WP) :: gr_DA_rms_rhoW = 0.0_WP
    real(WP) :: gr_DA_rms_rhoE = 0.0_WP
    real(WP) :: gr_DA_rms_p    = 0.0_WP

    real(WP), allocatable, dimension(:) :: base_rhos
    real(WP), allocatable, dimension(:) :: ddy_base_rhos
    real(WP), allocatable, dimension(:) :: ddx_base_rhos
    real(WP), allocatable, dimension(:) :: src_base_rhos
    real(WP), allocatable, dimension(:) :: gr_DA_rhos
    real(WP), allocatable, dimension(:) :: gr_DA_rms_rhos

    integer(c_int) :: ip

  end type largo_BL_spatiotemporal_workspace_type

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

  ! Tolerance to consider rms = 0
  real(WP), parameter :: eps = 1.0E-10_WP

  public  :: largo_BL_spatiotemporal_allocate
  public  :: largo_BL_spatiotemporal_deallocate
  public  :: largo_BL_spatiotemporal_init
  public  :: largo_BL_spatiotemporal_preStep_sEta
  public  :: largo_BL_spatiotemporal_preStep_sEta_innery
  public  :: largo_BL_spatiotemporal_preStep_sEta_innerxz
  public  :: largo_BL_spatiotemporal_preStep_sEtaMean
  public  :: largo_BL_spatiotemporal_continuity_sEtaMean
  public  :: largo_BL_spatiotemporal_xMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_yMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_zMomentum_sEtaMean
  public  :: largo_BL_spatiotemporal_energy_sEtaMean
  public  :: largo_BL_spatiotemporal_ispecies_sEtaMean
  public  :: largo_BL_spatiotemporal_species_sEtaMean
  public  :: largo_BL_spatiotemporal_continuity_sEtaRms
  public  :: largo_BL_spatiotemporal_xMomentum_sEtaRms
  public  :: largo_BL_spatiotemporal_yMomentum_sEtaRms
  public  :: largo_BL_spatiotemporal_zMomentum_sEtaRms
  public  :: largo_BL_spatiotemporal_energy_sEtaRms
  public  :: largo_BL_spatiotemporal_ispecies_sEtaRms
  public  :: largo_BL_spatiotemporal_species_sEtaRms

  private :: largo_BL_spatiotemporal_preStep_sEtaRms

  public  :: largo_BL_spatiotemporal_continuity_sEta
  public  :: largo_BL_spatiotemporal_xMomentum_sEta
  public  :: largo_BL_spatiotemporal_yMomentum_sEta
  public  :: largo_BL_spatiotemporal_zMomentum_sEta
  public  :: largo_BL_spatiotemporal_energy_sEta
  public  :: largo_BL_spatiotemporal_species_sEta

  public  :: largo_BL_spatiotemporal_sEtaMean
  public  :: largo_BL_spatiotemporal_sEta

  public  :: largo_BL_spatiotemporal_init_wall_baseflow
  public  :: largo_BL_spatiotemporal_preStep_baseflow

contains

  subroutine largo_BL_spatiotemporal_allocate(cp, neq, ns)

    type(largo_workspace_ptr), intent(out) :: cp
    ! neq=number of equations, might be needed later
    integer(c_int), intent(in)   :: neq
    integer(c_int), intent(in)   :: ns    ! number of species
    type(largo_BL_spatiotemporal_workspace_type), pointer :: auxp

    ! Allocate derived type variable
    allocate(auxp)

    ! Initialize values
    neq_ = neq
    ns_  = ns

    ! Allocate arrays for species
    if (ns_ > 0) then
      allocate(auxp%Xs_rhos    (1:ns_))
      allocate(auxp%ygrms_rhos (1:ns_))
      allocate(auxp%mean_rhos  (1:ns_))
      allocate(auxp%fav_cs     (1:ns_))
      allocate(auxp%fluc_rhos  (1:ns_))
      allocate(auxp%dtsRms_rhos(1:ns_))

      allocate(auxp%base_rhos     (1:ns_))
      allocate(auxp%ddy_base_rhos (1:ns_))
      allocate(auxp%ddx_base_rhos (1:ns_))
      allocate(auxp%src_base_rhos (1:ns_))
      allocate(auxp%gr_DA_rhos    (1:ns_))
      allocate(auxp%gr_DA_rms_rhos(1:ns_))

      auxp%base_rhos     = 0.0_WP
      auxp%ddy_base_rhos = 0.0_WP
      auxp%ddx_base_rhos = 0.0_WP
      auxp%src_base_rhos = 0.0_WP
      auxp%gr_DA_rhos    = 0.0_WP
      auxp%gr_DA_rms_rhos= 0.0_WP
    end if

    ! Pressure variable index
    auxp%ip = 5 + ns_ + 1

    ! Get C pointer from Fortran pointer
    cp = c_loc(auxp)

  end subroutine largo_BL_spatiotemporal_allocate


  subroutine largo_BL_spatiotemporal_deallocate(cp)

    type(largo_workspace_ptr), intent(inout)  :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer :: auxp

    call c_f_pointer(cp, auxp)

    ! Deallocate arrays for species
    if (allocated(auxp%Xs_rhos       ))  deallocate(auxp%Xs_rhos       )
    if (allocated(auxp%ygrms_rhos    ))  deallocate(auxp%ygrms_rhos    )
    if (allocated(auxp%mean_rhos     ))  deallocate(auxp%mean_rhos     )
    if (allocated(auxp%fav_cs        ))  deallocate(auxp%fav_cs        )
    if (allocated(auxp%fluc_rhos     ))  deallocate(auxp%fluc_rhos     )
    if (allocated(auxp%dtsRms_rhos   ))  deallocate(auxp%dtsRms_rhos   )

    if (allocated(auxp%base_rhos     ))  deallocate(auxp%base_rhos     )
    if (allocated(auxp%ddy_base_rhos ))  deallocate(auxp%ddy_base_rhos )
    if (allocated(auxp%ddx_base_rhos ))  deallocate(auxp%ddx_base_rhos )
    if (allocated(auxp%src_base_rhos ))  deallocate(auxp%src_base_rhos )
    if (allocated(auxp%gr_DA_rhos    ))  deallocate(auxp%gr_DA_rhos    )
    if (allocated(auxp%gr_DA_rms_rhos))  deallocate(auxp%gr_DA_rms_rhos)

    ! Deallocate array of derived types
    deallocate(auxp)

    ! Nullify C pointer
    cp = c_null_ptr

  end subroutine largo_BL_spatiotemporal_deallocate


  subroutine largo_BL_spatiotemporal_init(cp, grt_delta, gr_DA, gr_DA_rms)

    real(WP), intent(in)                  :: grt_delta
    real(WP), dimension(*), intent(in)    :: gr_DA
    real(WP), dimension(*), intent(in)    :: gr_DA_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Set growth rates
    auxp%grt_delta   = grt_delta

    auxp%gr_DA_rho  = gr_DA(irho )
    auxp%gr_DA_rhoU = gr_DA(irhoU)
    auxp%gr_DA_rhoV = gr_DA(irhoV)
    auxp%gr_DA_rhoW = gr_DA(irhoW)
    auxp%gr_DA_rhoE = gr_DA(irhoE)
    auxp%gr_DA_p    = gr_DA(auxp%ip)

    auxp%gr_DA_rms_rho  = gr_DA_rms(irho )
    auxp%gr_DA_rms_rhoU = gr_DA_rms(irhoU)
    auxp%gr_DA_rms_rhoV = gr_DA_rms(irhoV)
    auxp%gr_DA_rms_rhoW = gr_DA_rms(irhoW)
    auxp%gr_DA_rms_rhoE = gr_DA_rms(irhoE)
    auxp%gr_DA_rms_p    = gr_DA_rms(auxp%ip)

    do is=1, ns_
      auxp%gr_DA_rhos    (is) = gr_DA    (5+is)
      auxp%gr_DA_rms_rhos(is) = gr_DA_rms(5+is)
    end do

    ! Compute grx_delta using as velocity scale the 
    ! inviscid streamwise velocity at the wall
    ! NOTE: added the computation of grx_delta here as 
    ! well to avoid having as a requirement to call this method
    ! before the init_wall_baseflow one
    if (auxp%wall_base_u /= 0.0_WP) then
      auxp%grx_delta    = auxp%grt_delta   / auxp%wall_base_u 
    end if

  end subroutine largo_BL_spatiotemporal_init


  subroutine largo_BL_spatiotemporal_init_wall_baseflow(cp,  &   
               wall_base, wall_ddy_base, wall_ddt_base, wall_ddx_base)

    real(WP), dimension(*), intent(in)   :: wall_base
    real(WP), dimension(*), intent(in)   :: wall_ddy_base
    real(WP), dimension(*), intent(in)   :: wall_ddt_base
    real(WP), dimension(*), intent(in)   :: wall_ddx_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store relevant wall baseflow information
    auxp%wall_base_u  = wall_base(irhoU) / wall_base(irho)

    ! Compute grx_delta using as velocity scale the 
    ! inviscid streamwise velocity at the wall 
    auxp%grx_delta    = auxp%grt_delta   / auxp%wall_base_u 

  end subroutine largo_BL_spatiotemporal_init_wall_baseflow


  subroutine largo_BL_spatiotemporal_preStep_baseflow(cp,   &  
               base, ddy_base, ddt_base, ddx_base, src_base)

    real(WP), dimension(*), intent(in)    :: base
    real(WP), dimension(*), intent(in)    :: ddy_base
    real(WP), dimension(*), intent(in)    :: ddt_base
    real(WP), dimension(*), intent(in)    :: ddx_base
    real(WP), dimension(*), intent(in)    :: src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store baseflow information
    auxp%base_rho  = base(irho )
    auxp%base_rhoU = base(irhoU)
    auxp%base_rhoV = base(irhoV)
    auxp%base_rhoW = base(irhoW)
    auxp%base_rhoE = base(irhoE)
    auxp%base_p    = base(auxp%ip)

    auxp%ddy_base_rho  = ddy_base(irho )
    auxp%ddy_base_rhoU = ddy_base(irhoU)
    auxp%ddy_base_rhoV = ddy_base(irhoV)
    auxp%ddy_base_rhoW = ddy_base(irhoW)
    auxp%ddy_base_rhoE = ddy_base(irhoE)
    auxp%ddy_base_p    = ddy_base(auxp%ip)

    auxp%ddx_base_rho  = ddx_base(irho )
    auxp%ddx_base_rhoU = ddx_base(irhoU)
    auxp%ddx_base_rhoV = ddx_base(irhoV)
    auxp%ddx_base_rhoW = ddx_base(irhoW)
    auxp%ddx_base_rhoE = ddx_base(irhoE)
    auxp%ddx_base_p    = ddx_base(auxp%ip)

    auxp%src_base_rho  = src_base(irho )
    auxp%src_base_rhoU = src_base(irhoU)
    auxp%src_base_rhoV = src_base(irhoV)
    auxp%src_base_rhoW = src_base(irhoW)
    auxp%src_base_rhoE = src_base(irhoE)

    do is=1, ns_
      auxp%base_rhos    (is) =     base(5+is)
      auxp%ddy_base_rhos(is) = ddy_base(5+is)
      auxp%ddx_base_rhos(is) = ddx_base(5+is)
      auxp%src_base_rhos(is) = src_base(5+is)
    end do

  end subroutine largo_BL_spatiotemporal_preStep_baseflow


  subroutine largo_BL_spatiotemporal_preStep_sEtaMean(cp, y, mean, ddy_mean)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: mean
    real(WP), dimension(*), intent(in)    :: ddy_mean
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! These ones depend on y only
    auxp%mean_rho  = mean(irho )
    auxp%mean_rhoU = mean(irhoU)
    auxp%mean_rhoV = mean(irhoV)
    auxp%mean_rhoW = mean(irhoW)
    auxp%mean_rhoE = mean(irhoE)
    auxp%mean_p    = mean(auxp%ip)

    ! Compute/get mean and Favre averages
    auxp%fav_U    = mean(irhoU)/mean(irho)
    auxp%fav_V    = mean(irhoV)/mean(irho)
    auxp%fav_W    = mean(irhoW)/mean(irho)
    auxp%fav_E    = mean(irhoE)/mean(irho)
    do is=1, ns_
      auxp%mean_rhos(is) = mean(5+is)
      auxp%fav_cs   (is) = mean(5+is)/mean(irho)
    end do
    auxp%fav_H    = auxp%fav_E + auxp%mean_p/mean(irho)

    auxp%Xs_rho  = - auxp%ddx_base_rho  - auxp%gr_DA_rho  * (mean(irho )-auxp%base_rho ) + y * auxp%grx_delta * (ddy_mean(irho ) - auxp%ddy_base_rho ) !+ auxp%src_base_rho
    auxp%Xs_rhoU = - auxp%ddx_base_rhoU - auxp%gr_DA_rhoU * (mean(irhoU)-auxp%base_rhoU) + y * auxp%grx_delta * (ddy_mean(irhoU) - auxp%ddy_base_rhoU) !+ auxp%src_base_rhoU
    auxp%Xs_rhoV = - auxp%ddx_base_rhoV - auxp%gr_DA_rhoV * (mean(irhoV)-auxp%base_rhoV) + y * auxp%grx_delta * (ddy_mean(irhoV) - auxp%ddy_base_rhoV) !+ auxp%src_base_rhoV
    auxp%Xs_rhoW = - auxp%ddx_base_rhoW - auxp%gr_DA_rhoW * (mean(irhoW)-auxp%base_rhoW) + y * auxp%grx_delta * (ddy_mean(irhoW) - auxp%ddy_base_rhoW) !+ auxp%src_base_rhoW
    auxp%Xs_rhoE = - auxp%ddx_base_rhoE - auxp%gr_DA_rhoE * (mean(irhoE)-auxp%base_rhoE) + y * auxp%grx_delta * (ddy_mean(irhoE) - auxp%ddy_base_rhoE) !+ auxp%src_base_rhoE

    auxp%Xs_p    = - auxp%ddx_base_p    - auxp%gr_DA_p    * (mean(auxp%ip)-auxp%base_p ) + y * auxp%grx_delta * (ddy_mean(auxp%ip)- auxp%ddy_base_p)
    auxp%Xs_rhoH =   auxp%Xs_rhoE + auxp%Xs_p 

    do is=1, ns_
      auxp%Xs_rhos(is)  = - auxp%ddx_base_rhos(is) - auxp%gr_DA_rhos(is)  * (mean(5+is)-auxp%base_rhos(is)) + y * auxp%grx_delta * (ddy_mean(5+is) - auxp%ddy_base_rhos(is))
    end do

  end subroutine largo_BL_spatiotemporal_preStep_sEtaMean


  subroutine largo_BL_spatiotemporal_preStep_sEtaRms(cp, y, rms, ddy_rms)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in)    :: rms
    real(WP), dimension(*), intent(in)    :: ddy_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

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

  end subroutine largo_BL_spatiotemporal_preStep_sEtaRms


  subroutine largo_BL_spatiotemporal_preStep_sEta_innerxz(cp, qflow)

    real(WP), dimension(*), intent(in)    :: qflow
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_spatiotemporal_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    auxp%fluc_rho  = qflow(irho ) - auxp%mean_rho
    auxp%fluc_rhoU = qflow(irhoU) - auxp%mean_rhoU
    auxp%fluc_rhoV = qflow(irhoV) - auxp%mean_rhoV
    auxp%fluc_rhoW = qflow(irhoW) - auxp%mean_rhoW
    auxp%fluc_rhoE = qflow(irhoE) - auxp%mean_rhoE

    auxp%dtsRms_rho  = auxp%fluc_rho  * (- auxp%gr_DA_rms_rho  + auxp%ygrms_rho ) 
    auxp%dtsRms_rhoU = auxp%fluc_rhoU * (- auxp%gr_DA_rms_rhoU + auxp%ygrms_rhoU)
    auxp%dtsRms_rhoV = auxp%fluc_rhoV * (- auxp%gr_DA_rms_rhoV + auxp%ygrms_rhoV)
    auxp%dtsRms_rhoW = auxp%fluc_rhoW * (- auxp%gr_DA_rms_rhoW + auxp%ygrms_rhoW)
    auxp%dtsRms_rhoE = auxp%fluc_rhoE * (- auxp%gr_DA_rms_rhoE + auxp%ygrms_rhoE)

    do is=1, ns_
      auxp%fluc_rhos(is)  = qflow(5+is) - auxp%mean_rhos(is)
      auxp%dtsRms_rhos(is) = auxp%fluc_rhos(is) * (- auxp%gr_DA_rms_rhos(is) + auxp%ygrms_rhos(is))
    end do

  end subroutine largo_BL_spatiotemporal_preStep_sEta_innerxz


  subroutine largo_BL_spatiotemporal_preStep_sEta_innery(cp, y,                  &
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

    call largo_BL_spatiotemporal_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_spatiotemporal_preStep_sEtaRms(cp, y, rms, ddy_rms)

  end subroutine largo_BL_spatiotemporal_preStep_sEta_innery


  subroutine largo_BL_spatiotemporal_preStep_sEta(cp, y, qflow,                    &
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


    call largo_BL_spatiotemporal_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_spatiotemporal_preStep_sEtaRms(cp, y, rms, ddy_rms)
    call largo_BL_spatiotemporal_preStep_sEta_innerxz(cp, qflow)

  end subroutine largo_BL_spatiotemporal_preStep_sEta


#define DECLARE_SUBROUTINE(token)token (cp, A, B, src);\
  type(largo_workspace_ptr), intent(in)  :: cp;\
  real(WP)       , intent(in)            :: A, B;\
  real(WP)       , intent(inout)         :: src;\
  type(largo_BL_spatiotemporal_workspace_type), pointer    :: auxp;\
  call c_f_pointer(cp, auxp);\


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_continuity_sEtaMean)
    src = A * src + B * auxp%Xs_rhoU
  end subroutine largo_BL_spatiotemporal_continuity_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_xMomentum_sEtaMean)
    src = A * src + B * (2.0_wp * auxp%fav_u * auxp%Xs_rhoU - auxp%fav_u**2 * auxp%Xs_rho + auxp%Xs_p)
  end subroutine largo_BL_spatiotemporal_xMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_yMomentum_sEtaMean)
    src = A * src + B * (auxp%fav_u * auxp%Xs_rhoV + auxp%fav_v * auxp%Xs_rhoU - auxp%fav_u*auxp%fav_v * auxp%Xs_rho )
  end subroutine largo_BL_spatiotemporal_yMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_zMomentum_sEtaMean)
    src = A * src + B * (auxp%fav_u * auxp%Xs_rhoW + auxp%fav_w * auxp%Xs_rhoU - auxp%fav_u*auxp%fav_w * auxp%Xs_rho )
  end subroutine largo_BL_spatiotemporal_zMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_energy_sEtaMean)
    src = A * src + B * (auxp%fav_u * auxp%Xs_rhoH + auxp%fav_H * auxp%Xs_rhoU - auxp%fav_u*auxp%fav_H*auxp%Xs_rho)
  end subroutine largo_BL_spatiotemporal_energy_sEtaMean


  subroutine largo_BL_spatiotemporal_ispecies_sEtaMean (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)   :: cp
    real(WP)       , intent(in)             :: A, B
    real(WP)       , intent(inout)          :: src
    type(largo_BL_spatiotemporal_workspace_type), pointer     :: auxp
    integer(c_int), intent(in)              :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * (auxp%fav_u * auxp%Xs_rhos(is) + auxp%fav_cs(is) * auxp%Xs_rhoU - auxp%fav_u * auxp%fav_cs(is) * auxp%Xs_rho )
  end subroutine largo_BL_spatiotemporal_ispecies_sEtaMean


  subroutine largo_BL_spatiotemporal_species_sEtaMean (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_spatiotemporal_ispecies_sEtaMean (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_spatiotemporal_species_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_continuity_sEtaRms)
    src = A * src + B * auxp%dtsRms_rho
  end subroutine largo_BL_spatiotemporal_continuity_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_xMomentum_sEtaRms)
    src = A * src + B * auxp%dtsRms_rhoU
  end subroutine largo_BL_spatiotemporal_xMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_yMomentum_sEtaRms)
    src = A * src + B * auxp%dtsRms_rhoV
  end subroutine largo_BL_spatiotemporal_yMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_zMomentum_sEtaRms)
    src = A * src + B * auxp%dtsRms_rhoW
  end subroutine largo_BL_spatiotemporal_zMomentum_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_energy_sEtaRms)
    src = A * src + B * auxp%dtsRms_rhoE
  end subroutine largo_BL_spatiotemporal_energy_sEtaRms


  subroutine largo_BL_spatiotemporal_ispecies_sEtaRms (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP)       , intent(inout)            :: src
    type(largo_BL_spatiotemporal_workspace_type), pointer       :: auxp
    integer(c_int), intent(in)                :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * auxp%dtsRms_rhos(is)
  end subroutine largo_BL_spatiotemporal_ispecies_sEtaRms


  subroutine largo_BL_spatiotemporal_species_sEtaRms (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_spatiotemporal_ispecies_sEtaRms (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_spatiotemporal_species_sEtaRms


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_continuity_sEta)
    call largo_BL_spatiotemporal_continuity_sEtaMean (cp,      A, B, src)
    call largo_BL_spatiotemporal_continuity_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_spatiotemporal_continuity_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_xMomentum_sEta)
    call largo_BL_spatiotemporal_xMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_spatiotemporal_xMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_spatiotemporal_xMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_yMomentum_sEta)
    call largo_BL_spatiotemporal_yMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_spatiotemporal_yMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_spatiotemporal_yMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_zMomentum_sEta)
    call largo_BL_spatiotemporal_zMomentum_sEtaMean (cp,      A, B, src)
    call largo_BL_spatiotemporal_zMomentum_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_spatiotemporal_zMomentum_sEta


  subroutine DECLARE_SUBROUTINE(largo_BL_spatiotemporal_energy_sEta)
    call largo_BL_spatiotemporal_energy_sEtaMean (cp,      A, B, src)
    call largo_BL_spatiotemporal_energy_sEtaRms  (cp, 1.0_WP, B, src)
  end subroutine largo_BL_spatiotemporal_energy_sEta


  subroutine largo_BL_spatiotemporal_species_sEta (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_spatiotemporal_ispecies_sEtaMean (cp,      A, B, srcvec(is), is)
      call largo_BL_spatiotemporal_ispecies_sEtaRms  (cp, 1.0_WP, B, srcvec(is), is)
    end do
  end subroutine largo_BL_spatiotemporal_species_sEta


  subroutine largo_BL_spatiotemporal_sEtaMean (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_spatiotemporal_continuity_sEtaMean (cp,      A, B, srcvec(irho ))
    call largo_BL_spatiotemporal_xMomentum_sEtaMean  (cp,      A, B, srcvec(irhou))
    call largo_BL_spatiotemporal_yMomentum_sEtaMean  (cp,      A, B, srcvec(irhov))
    call largo_BL_spatiotemporal_zMomentum_sEtaMean  (cp,      A, B, srcvec(irhow))
    call largo_BL_spatiotemporal_energy_sEtaMean     (cp,      A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_spatiotemporal_ispecies_sEtaMean (cp,      A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_spatiotemporal_sEtaMean


  subroutine largo_BL_spatiotemporal_sEta (cp, A, B, srcvec) bind(C)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_spatiotemporal_continuity_sEtaMean (cp,      A, B, srcvec(irho ))
    call largo_BL_spatiotemporal_continuity_sEtaRms  (cp, 1.0_WP, B, srcvec(irho ))

    call largo_BL_spatiotemporal_xMomentum_sEtaMean  (cp,      A, B, srcvec(irhou))
    call largo_BL_spatiotemporal_xMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhou))

    call largo_BL_spatiotemporal_yMomentum_sEtaMean  (cp,      A, B, srcvec(irhov))
    call largo_BL_spatiotemporal_yMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhov))

    call largo_BL_spatiotemporal_zMomentum_sEtaMean  (cp,      A, B, srcvec(irhow))
    call largo_BL_spatiotemporal_zMomentum_sEtaRms   (cp, 1.0_WP, B, srcvec(irhow))

    call largo_BL_spatiotemporal_energy_sEtaMean     (cp,      A, B, srcvec(irhoE))
    call largo_BL_spatiotemporal_energy_sEtaRms      (cp, 1.0_WP, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_spatiotemporal_ispecies_sEtaMean (cp,      A, B, srcvec(5+is), is)
      call largo_BL_spatiotemporal_ispecies_sEtaRms  (cp, 1.0_WP, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_spatiotemporal_sEta

end module largo_BL_spatiotemporal
