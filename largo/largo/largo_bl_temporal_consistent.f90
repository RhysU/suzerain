module largo_BL_temporal_consistent

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

  type :: largo_BL_temporal_consistent_workspace_type

    real(WP) :: gr_delta   = 1.0_WP

    real(WP) :: Ts_rho = 0.0_WP
    real(WP) :: Ts_U   = 0.0_WP
    real(WP) :: Ts_V   = 0.0_WP
    real(WP) :: Ts_W   = 0.0_WP
    real(WP) :: Ts_E   = 0.0_WP

    real(WP) :: TsArms_rho = 0.0_WP
    real(WP) :: TsArms_U   = 0.0_WP
    real(WP) :: TsArms_V   = 0.0_WP
    real(WP) :: TsArms_W   = 0.0_WP
    real(WP) :: TsArms_E   = 0.0_WP

    real(WP) :: TsFull_rho = 0.0_WP
    real(WP) :: TsFull_U   = 0.0_WP
    real(WP) :: TsFull_V   = 0.0_WP
    real(WP) :: TsFull_W   = 0.0_WP
    real(WP) :: TsFull_E   = 0.0_WP

    real(WP), allocatable, dimension(:) :: Ts_cs
    real(WP), allocatable, dimension(:) :: ffluc_cs
    real(WP), allocatable, dimension(:) :: TsArms_cs
    real(WP), allocatable, dimension(:) :: TsFull_cs


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

    real(WP) :: gr_DA_rho = 0.0_WP
    real(WP) :: gr_DA_U   = 0.0_WP
    real(WP) :: gr_DA_V   = 0.0_WP
    real(WP) :: gr_DA_W   = 0.0_WP
    real(WP) :: gr_DA_E   = 0.0_WP

    real(WP) :: gr_DA_rms_rho = 0.0_WP
    real(WP) :: gr_DA_rms_U   = 0.0_WP
    real(WP) :: gr_DA_rms_V   = 0.0_WP
    real(WP) :: gr_DA_rms_W   = 0.0_WP
    real(WP) :: gr_DA_rms_E   = 0.0_WP

    real(WP) :: ygArms_rho = 0.0_WP
    real(WP) :: ygArms_U   = 0.0_WP
    real(WP) :: ygArms_V   = 0.0_WP
    real(WP) :: ygArms_W   = 0.0_WP
    real(WP) :: ygArms_E   = 0.0_WP

    real(WP), allocatable, dimension(:) ::  field_cs
    real(WP), allocatable, dimension(:) ::    fav_cs
    real(WP), allocatable, dimension(:) ::   dfav_cs
    real(WP), allocatable, dimension(:) ::   Arms_cs
    real(WP), allocatable, dimension(:) ::  dArms_cs
    real(WP), allocatable, dimension(:) :: ygArms_cs
    real(WP), allocatable, dimension(:) ::  gr_DA_cs
    real(WP), allocatable, dimension(:) :: gr_DA_rms_cs

    ! baseflow variables
    real(WP) :: base_rho   = 0.0_WP
    real(WP) :: base_U     = 0.0_WP
    real(WP) :: base_V     = 0.0_WP
    real(WP) :: base_W     = 0.0_WP
    real(WP) :: base_E     = 0.0_WP

    real(WP) :: ddy_base_rho = 0.0_WP
    real(WP) :: ddy_base_U   = 0.0_WP
    real(WP) :: ddy_base_V   = 0.0_WP
    real(WP) :: ddy_base_W   = 0.0_WP
    real(WP) :: ddy_base_E   = 0.0_WP

    real(WP) :: ddt_base_rho = 0.0_WP
    real(WP) :: ddt_base_U   = 0.0_WP
    real(WP) :: ddt_base_V   = 0.0_WP
    real(WP) :: ddt_base_W   = 0.0_WP
    real(WP) :: ddt_base_E   = 0.0_WP

    ! FIXME: Flow volumetric sources should be added by
    !        the user, they are not being added here.
    !real(WP) :: src_base_rho = 0.0_WP
    !real(WP) :: src_base_U   = 0.0_WP
    !real(WP) :: src_base_V   = 0.0_WP
    !real(WP) :: src_base_W   = 0.0_WP
    !real(WP) :: src_base_E   = 0.0_WP

    real(WP), allocatable, dimension(:) :: base_cs
    real(WP), allocatable, dimension(:) :: ddy_base_cs
    real(WP), allocatable, dimension(:) :: ddt_base_cs
    real(WP), allocatable, dimension(:) :: src_base_cs

  end type largo_BL_temporal_consistent_workspace_type

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


  public  :: largo_BL_temporal_consistent_allocate
  public  :: largo_BL_temporal_consistent_deallocate
  public  :: largo_BL_temporal_consistent_init
  public  :: largo_BL_temporal_consistent_preStep_sEta
  public  :: largo_BL_temporal_consistent_preStep_sEta_innery
  public  :: largo_BL_temporal_consistent_preStep_sEta_innerxz
  public  :: largo_BL_temporal_consistent_preStep_sEtaMean
  public  :: largo_BL_temporal_consistent_continuity_sEtaMean
  public  :: largo_BL_temporal_consistent_xMomentum_sEtaMean
  public  :: largo_BL_temporal_consistent_yMomentum_sEtaMean
  public  :: largo_BL_temporal_consistent_zMomentum_sEtaMean
  public  :: largo_BL_temporal_consistent_energy_sEtaMean
  public  :: largo_BL_temporal_consistent_ispecies_sEtaMean
  public  :: largo_BL_temporal_consistent_species_sEtaMean
  public  :: largo_BL_temporal_consistent_continuity_sEta_
  public  :: largo_BL_temporal_consistent_xMomentum_sEta_
  public  :: largo_BL_temporal_consistent_yMomentum_sEta_
  public  :: largo_BL_temporal_consistent_zMomentum_sEta_
  public  :: largo_BL_temporal_consistent_energy_sEta_
  public  :: largo_BL_temporal_consistent_ispecies_sEta_
  public  :: largo_BL_temporal_consistent_species_sEta_
  public  :: largo_BL_temporal_consistent_sEta

  private :: largo_BL_temporal_consistent_preStep_sEta_

  public  :: largo_BL_temporal_consistent_preStep_baseflow

contains

  subroutine largo_BL_temporal_consistent_allocate(cp, neq, ns)

    type(largo_workspace_ptr), intent(out) :: cp
    ! neq=number of equations, might be needed later
    integer(c_int), intent(in)   :: neq
    integer(c_int), intent(in)   :: ns    ! number of species
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp

    ! Allocate derived type variable
    allocate(auxp)

    ! Initialize values
    neq_ = neq
    ns_  = ns

    ! Allocate arrays for species
    if (ns_ > 0) then
      allocate(auxp%Ts_cs        (1:ns_))
      allocate(auxp%ffluc_cs     (1:ns_))
      allocate(auxp%fav_cs       (1:ns_))
      allocate(auxp%dfav_cs      (1:ns_))
      allocate(auxp%field_cs     (1:ns_))
      allocate(auxp%TsArms_cs    (1:ns_))
      allocate(auxp%Arms_cs      (1:ns_))
      allocate(auxp%dArms_cs     (1:ns_))
      allocate(auxp%ygArms_cs    (1:ns_))
      allocate(auxp%TsFull_cs    (1:ns_))
      allocate(auxp%gr_DA_cs     (1:ns_))
      allocate(auxp%gr_DA_rms_cs (1:ns_))

      allocate(auxp%base_cs      (1:ns_))
      allocate(auxp%ddy_base_cs  (1:ns_))
      allocate(auxp%ddt_base_cs  (1:ns_))
      allocate(auxp%src_base_cs  (1:ns_))

      auxp%Arms_cs      = 0.0_WP
      auxp%dArms_cs     = 0.0_WP
      auxp%gr_DA_cs     = 0.0_WP
      auxp%gr_DA_rms_cs = 0.0_WP
      auxp%base_cs      = 0.0_WP
      auxp%ddy_base_cs  = 0.0_WP
      auxp%ddt_base_cs  = 0.0_WP
      auxp%src_base_cs  = 0.0_WP
    end if

    ! Get C pointer from Fortran pointer
    cp = c_loc(auxp)

  end subroutine largo_BL_temporal_consistent_allocate


  subroutine largo_BL_temporal_consistent_deallocate(cp)

    type(largo_workspace_ptr), intent(out)  :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp

    call c_f_pointer(cp, auxp)

    ! Deallocate arrays for species
    if (allocated(auxp%Ts_cs))        deallocate(auxp%Ts_cs       )
    if (allocated(auxp%ffluc_cs))     deallocate(auxp%ffluc_cs    )
    if (allocated(auxp%fav_cs))       deallocate(auxp%fav_cs      )
    if (allocated(auxp%dfav_cs))      deallocate(auxp%dfav_cs     )
    if (allocated(auxp%field_cs))     deallocate(auxp%field_cs    )
    if (allocated(auxp%TsArms_cs))    deallocate(auxp%TsArms_cs   )
    if (allocated(auxp%Arms_cs))      deallocate(auxp%Arms_cs     )
    if (allocated(auxp%dArms_cs))     deallocate(auxp%dArms_cs    )
    if (allocated(auxp%ygArms_cs))    deallocate(auxp%ygArms_cs   )
    if (allocated(auxp%TsFull_cs))    deallocate(auxp%TsFull_cs   )
    if (allocated(auxp%gr_DA_cs))     deallocate(auxp%gr_DA_cs    )
    if (allocated(auxp%gr_DA_rms_cs)) deallocate(auxp%gr_DA_rms_cs)

    if (allocated(auxp%base_cs    ))  deallocate(auxp%base_cs     )
    if (allocated(auxp%ddy_base_cs )) deallocate(auxp%ddy_base_cs )
    if (allocated(auxp%ddt_base_cs )) deallocate(auxp%ddt_base_cs )
    if (allocated(auxp%src_base_cs )) deallocate(auxp%src_base_cs )

    ! Deallocate array of derived types
    deallocate(auxp)

    ! Nullify C pointer
    cp = c_null_ptr

  end subroutine largo_BL_temporal_consistent_deallocate


  subroutine largo_BL_temporal_consistent_init(cp, gr_delta, gr_DA, gr_DA_rms)

    real(WP), intent(in)                  :: gr_delta
    real(WP), dimension(*), intent(in)    :: gr_DA
    real(WP), dimension(*), intent(in)    :: gr_DA_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer   :: auxp

    ! get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! set growth rates
    auxp%gr_delta   = gr_delta

    ! Mean defect amplitude growth rates
    auxp%gr_DA_rho  = gr_DA(irho )
    auxp%gr_DA_U    = gr_DA(irhoU)
    auxp%gr_DA_V    = gr_DA(irhoV)
    auxp%gr_DA_W    = gr_DA(irhoW)
    auxp%gr_DA_E    = gr_DA(irhoE)

    do is=1, ns_
      auxp%gr_DA_cs(is) = gr_DA(5+is)
    end do

    ! rms amplitude growth rates
    auxp%gr_DA_rms_rho  = gr_DA_rms(irho )
    auxp%gr_DA_rms_U    = gr_DA_rms(irhoU)
    auxp%gr_DA_rms_V    = gr_DA_rms(irhoV)
    auxp%gr_DA_rms_W    = gr_DA_rms(irhoW)
    auxp%gr_DA_rms_E    = gr_DA_rms(irhoE)

    do is=1, ns_
      auxp%gr_DA_rms_cs(is) = gr_DA_rms(5+is)
    end do

  end subroutine largo_BL_temporal_consistent_init


  subroutine largo_BL_temporal_consistent_preStep_baseflow(cp, base, &
                              ddy_base, ddt_base, ddx_base, src_base)

    real(WP), dimension(*), intent(in)    :: base
    real(WP), dimension(*), intent(in)    :: ddy_base
    real(WP), dimension(*), intent(in)    :: ddt_base
    real(WP), dimension(*), intent(in)    :: ddx_base
    real(WP), dimension(*), intent(in)    :: src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer   :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Store baseflow information
    if (base(irho) > 0.0_WP) then
      auxp%base_rho = base(irho )
      auxp%base_U   = base(irhoU)/base(irho )
      auxp%base_V   = base(irhoV)/base(irho )
      auxp%base_W   = base(irhoW)/base(irho )
      auxp%base_E   = base(irhoE)/base(irho )

      auxp%ddy_base_rho = ddy_base(irho )
      auxp%ddy_base_U   = ddy_base(irhoU)/base(irho ) - auxp%base_U/base(irho ) * ddy_base(irho )
      auxp%ddy_base_V   = ddy_base(irhoV)/base(irho ) - auxp%base_V/base(irho ) * ddy_base(irho )
      auxp%ddy_base_W   = ddy_base(irhoW)/base(irho ) - auxp%base_W/base(irho ) * ddy_base(irho )
      auxp%ddy_base_E   = ddy_base(irhoE)/base(irho ) - auxp%base_E/base(irho ) * ddy_base(irho )

      auxp%ddt_base_rho = ddt_base(irho )
      auxp%ddt_base_U   = ddt_base(irhoU)/base(irho ) - auxp%base_U/base(irho ) * ddt_base(irho )
      auxp%ddt_base_V   = ddt_base(irhoV)/base(irho ) - auxp%base_V/base(irho ) * ddt_base(irho )
      auxp%ddt_base_W   = ddt_base(irhoW)/base(irho ) - auxp%base_W/base(irho ) * ddt_base(irho )
      auxp%ddt_base_E   = ddt_base(irhoE)/base(irho ) - auxp%base_E/base(irho ) * ddt_base(irho )
    end if


!!$     auxp%src_base_rho = src_base(irho )
!!$     auxp%src_base_U   = src_base(irhoU)
!!$     auxp%src_base_V   = src_base(irhoV)
!!$     auxp%src_base_W   = src_base(irhoW)
!!$     auxp%src_base_E   = src_base(irhoE)
!!$
    do is=1, ns_
      auxp%base_cs    (is) =     base(5+is)/base(irho )
      auxp%ddy_base_cs(is) = ddy_base(5+is)/base(irho ) - base(5+is)/base(irho )**2 * ddy_base(irho )
      auxp%ddt_base_cs(is) = ddt_base(5+is)/base(irho ) - base(5+is)/base(irho )**2 * ddt_base(irho )
!!$       auxp%src_base_cs(is) = src_base(5+is)
    end do

  end subroutine largo_BL_temporal_consistent_preStep_baseflow


  subroutine largo_BL_temporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: ddy_mean
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp

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
    auxp%Ts_rho =  - auxp%ddt_base_rho - auxp%gr_DA_rho * (auxp%mean_rho-auxp%base_rho) + y * auxp%gr_delta * (ddy_mean(irho )-auxp%ddy_base_rho)
    auxp%Ts_U   =  - auxp%ddt_base_U   - auxp%gr_DA_U   * (auxp%fav_U   -auxp%base_U  ) + y * auxp%gr_delta * (auxp%dfav_U    -auxp%ddy_base_U  )
    auxp%Ts_V   =  - auxp%ddt_base_V   - auxp%gr_DA_V   * (auxp%fav_V   -auxp%base_V  ) + y * auxp%gr_delta * (auxp%dfav_V    -auxp%ddy_base_V  )
    auxp%Ts_W   =  - auxp%ddt_base_W   - auxp%gr_DA_W   * (auxp%fav_W   -auxp%base_W  ) + y * auxp%gr_delta * (auxp%dfav_W    -auxp%ddy_base_W  )
    auxp%Ts_E   =  - auxp%ddt_base_E   - auxp%gr_DA_E   * (auxp%fav_E   -auxp%base_E  ) + y * auxp%gr_delta * (auxp%dfav_E    -auxp%ddy_base_E  )

    do is=1, ns_
      auxp%Ts_cs(is)  = - auxp%ddt_base_cs(is) - auxp%gr_DA_cs(is) * (auxp%fav_cs(is)-auxp%base_cs(is)) + y * auxp%gr_delta * (auxp%dfav_cs(is)-auxp%ddy_base_cs(is)) 
    end do

  end subroutine largo_BL_temporal_consistent_preStep_sEtaMean


  subroutine largo_BL_temporal_consistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: dmean_rqq
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in)          :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Compute \mean{ru"u"} for each u component and for energy
    auxp%rhoupup  = mean_rqq(irhoU) &
      &          - auxp%mean_rho * auxp%fav_u**2

    auxp%rhovpvp  = mean_rqq(irhoV) &
      &          - auxp%mean_rho * auxp%fav_v**2

    auxp%rhowpwp  = mean_rqq(irhoW) &
      &          - auxp%mean_rho * auxp%fav_w**2

    auxp%rhoEpEp  = mean_rqq(irhoE) &
      &          - auxp%mean_rho * auxp%fav_E**2

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

    auxp%ygArms_U = 0.0_WP
    if (auxp%Arms_U   > eps) auxp%ygArms_U   = y * auxp%gr_delta * auxp%dArms_U/auxp%Arms_U

    auxp%ygArms_V = 0.0_WP
    if (auxp%Arms_V   > eps) auxp%ygArms_V   = y * auxp%gr_delta * auxp%dArms_V/auxp%Arms_V

    auxp%ygArms_W = 0.0_WP
    if (auxp%Arms_W   > eps) auxp%ygArms_W   = y * auxp%gr_delta * auxp%dArms_W/auxp%Arms_W

    auxp%ygArms_E = 0.0_WP
    if (auxp%Arms_E   > eps) auxp%ygArms_E   = y * auxp%gr_delta * auxp%dArms_E/auxp%Arms_E

    do is=1, ns_
      auxp%ygArms_cs(is) = 0.0_WP
      if (auxp%Arms_cs(is) > eps)  auxp%ygArms_cs(is) = y * auxp%gr_delta * auxp%dArms_cs(is)/auxp%Arms_cs(is)
    end do

  end subroutine largo_BL_temporal_consistent_preStep_sEta_


  subroutine largo_BL_temporal_consistent_preStep_sEta_innerxz(cp, qflow)

    real(WP), dimension(*), intent(in) :: qflow
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp

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

    auxp%TsArms_rho = auxp%fluc_rho / auxp%mean_rho * auxp%Ts_rho
    auxp%TsArms_U   = auxp%ffluc_U  * (- auxp%gr_DA_rms_U   + auxp%ygArms_U   )
    auxp%TsArms_V   = auxp%ffluc_V  * (- auxp%gr_DA_rms_V   + auxp%ygArms_V   )
    auxp%TsArms_W   = auxp%ffluc_W  * (- auxp%gr_DA_rms_W   + auxp%ygArms_W   )
    auxp%TsArms_E   = auxp%ffluc_E  * (- auxp%gr_DA_rms_E   + auxp%ygArms_E   )

    do is=1, ns_
      auxp%field_cs (is)  = qflow(5+is)/qflow(irho)
      auxp%ffluc_cs (is)  = auxp%field_cs(is) - auxp%fav_cs   (is)
      auxp%TsArms_cs(is)  = auxp%ffluc_cs(is) * (- auxp%gr_DA_rms_cs(is) + auxp%ygArms_cs(is))
    end do

    ! Compute mean plus fluctuations slow time derivative
    auxp%TsFull_rho = auxp%Ts_rho + auxp%TsArms_rho
    auxp%TsFull_U   = auxp%Ts_U   + auxp%TsArms_U
    auxp%TsFull_V   = auxp%Ts_V   + auxp%TsArms_V
    auxp%TsFull_W   = auxp%Ts_W   + auxp%TsArms_W
    auxp%TsFull_E   = auxp%Ts_E   + auxp%TsArms_E

    do is=1, ns_
      auxp%TsFull_cs(is) = auxp%Ts_cs(is) + auxp%TsArms_cs(is) 
    end do

  end subroutine largo_BL_temporal_consistent_preStep_sEta_innerxz


  subroutine largo_BL_temporal_consistent_preStep_sEta_innery(cp, y, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: ddy_mean
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: dmean_rqq
    type(largo_workspace_ptr), intent(in)        :: cp

    call largo_BL_temporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_temporal_consistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

  end subroutine largo_BL_temporal_consistent_preStep_sEta_innery


  subroutine largo_BL_temporal_consistent_preStep_sEta(cp, y, qflow, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

    real(WP), intent(in)                  :: y
    real(WP), dimension(*), intent(in) :: qflow
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: ddy_mean
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: dmean_rqq
    type(largo_workspace_ptr), intent(in)        :: cp


    call largo_BL_temporal_consistent_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_temporal_consistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)
    call largo_BL_temporal_consistent_preStep_sEta_innerxz(cp, qflow)

  end subroutine largo_BL_temporal_consistent_preStep_sEta


#define DECLARE_SUBROUTINE(token)token (cp, A, B, src);\
  type(largo_workspace_ptr), intent(in)   :: cp;\
  real(WP)       , intent(in)             :: A, B;\
  real(WP)       , intent(inout)          :: src;\
  type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp;\
  call c_f_pointer(cp, auxp)\

  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_continuity_sEtaMean)
    src = A * src + B * auxp%Ts_rho
  end subroutine largo_BL_temporal_consistent_continuity_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_xMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%Ts_U &
      &                  + auxp%fav_U    * auxp%Ts_rho  )
  end subroutine largo_BL_temporal_consistent_xMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_yMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%Ts_V &
      &                  + auxp%fav_V    * auxp%Ts_rho  )
  end subroutine largo_BL_temporal_consistent_yMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_zMomentum_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%Ts_W &
      &                  + auxp%fav_W    * auxp%Ts_rho  )
  end subroutine largo_BL_temporal_consistent_zMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_energy_sEtaMean)
    src = A * src + B * (  auxp%mean_rho * auxp%Ts_E &
      &                  + auxp%fav_E    * auxp%Ts_rho  )
  end subroutine largo_BL_temporal_consistent_energy_sEtaMean


  subroutine largo_BL_temporal_consistent_ispecies_sEtaMean (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * (  auxp%mean_rho   * auxp%Ts_cs(is) &
      &                  + auxp%fav_cs(is) * auxp%Ts_rho  )
  end subroutine largo_BL_temporal_consistent_ispecies_sEtaMean


  subroutine largo_BL_temporal_consistent_species_sEtaMean (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_temporal_consistent_ispecies_sEtaMean (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_consistent_species_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_continuity_sEta_)
    src = A * src + B * auxp%TsFull_rho
  end subroutine largo_BL_temporal_consistent_continuity_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_xMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%TsFull_U &
      &                  + auxp%field_U   * auxp%TsFull_rho  )
  end subroutine largo_BL_temporal_consistent_xMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_yMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%TsFull_V &
      &                  + auxp%field_V   * auxp%TsFull_rho  )
  end subroutine largo_BL_temporal_consistent_yMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_zMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%TsFull_W &
      &                  + auxp%field_W   * auxp%TsFull_rho  )
  end subroutine largo_BL_temporal_consistent_zMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_consistent_energy_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%TsFull_E &
      &                  + auxp%field_E   * auxp%TsFull_rho  )
  end subroutine largo_BL_temporal_consistent_energy_sEta_


  subroutine largo_BL_temporal_consistent_ispecies_sEta_ (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_temporal_consistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * (  auxp%field_rho    * auxp%TsFull_cs(is) &
      &                  + auxp%field_cs(is) * auxp%TsFull_rho  )
  end subroutine largo_BL_temporal_consistent_ispecies_sEta_


  subroutine largo_BL_temporal_consistent_species_sEta_ (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)   :: cp
    real(WP)       , intent(in)             :: A, B
    real(WP), dimension(*), intent(inout)   :: srcvec ! "*" = ns_
    integer(c_int)                          :: is

    do is = 1, ns_
      call largo_BL_temporal_consistent_ispecies_sEta_ (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_consistent_species_sEta_


  subroutine largo_BL_temporal_consistent_sEta (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_temporal_consistent_continuity_sEta_ (cp, A, B, srcvec(irho ))
    call largo_BL_temporal_consistent_xMomentum_sEta_  (cp, A, B, srcvec(irhou))
    call largo_BL_temporal_consistent_yMomentum_sEta_  (cp, A, B, srcvec(irhov))
    call largo_BL_temporal_consistent_zMomentum_sEta_  (cp, A, B, srcvec(irhow))
    call largo_BL_temporal_consistent_energy_sEta_     (cp, A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_temporal_consistent_ispecies_sEta_ (cp, A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_temporal_consistent_sEta

end module largo_BL_temporal_consistent
