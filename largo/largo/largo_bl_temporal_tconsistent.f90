module largo_BL_temporal_tconsistent

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

  type :: largo_BL_temporal_tconsistent_workspace_type

    real(WP) :: gr_delta   = 1.0_WP

    real(WP) :: dts_rho = 0.0_WP
    real(WP) :: dts_U   = 0.0_WP
    real(WP) :: dts_V   = 0.0_WP
    real(WP) :: dts_W   = 0.0_WP
    real(WP) :: dts_E   = 0.0_WP

    real(WP) :: dts_rhoU = 0.0_WP
    real(WP) :: dts_rhoV = 0.0_WP
    real(WP) :: dts_rhoW = 0.0_WP
    real(WP) :: dts_rhoE = 0.0_WP

!!$     real(WP) :: ygrms_rho   = 0.0_WP
!!$     real(WP) :: ygrms_rhoU  = 0.0_WP
!!$     real(WP) :: ygrms_rhoV  = 0.0_WP
!!$     real(WP) :: ygrms_rhoW  = 0.0_WP
!!$     real(WP) :: ygrms_rhoE  = 0.0_WP
!!$
!!$     real(WP) :: fluc_rho   = 0.0_WP
!!$     real(WP) :: fluc_rhoU  = 0.0_WP
!!$     real(WP) :: fluc_rhoV  = 0.0_WP
!!$     real(WP) :: fluc_rhoW  = 0.0_WP
!!$     real(WP) :: fluc_rhoE  = 0.0_WP

    real(WP) :: dtsArms_rho = 0.0_WP
    real(WP) :: dtsArms_U   = 0.0_WP
    real(WP) :: dtsArms_V   = 0.0_WP
    real(WP) :: dtsArms_W   = 0.0_WP
    real(WP) :: dtsArms_E   = 0.0_WP

    real(WP) :: dtsFull_rho = 0.0_WP
    real(WP) :: dtsFull_U   = 0.0_WP
    real(WP) :: dtsFull_V   = 0.0_WP
    real(WP) :: dtsFull_W   = 0.0_WP
    real(WP) :: dtsFull_E   = 0.0_WP

    real(WP), allocatable, dimension(:) :: dts_rhos
    real(WP), allocatable, dimension(:) :: ygrms_rhos
    real(WP), allocatable, dimension(:) :: fluc_rhos
    real(WP), allocatable, dimension(:) :: dtsArms_rhos


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

    real(WP) :: ygArms_rho = 0.0_WP
    real(WP) :: ygArms_U   = 0.0_WP
    real(WP) :: ygArms_V   = 0.0_WP
    real(WP) :: ygArms_W   = 0.0_WP
    real(WP) :: ygArms_E   = 0.0_WP

    real(WP), allocatable, dimension(:) ::  field_rhos
    real(WP), allocatable, dimension(:) ::   mean_rhos
    real(WP), allocatable, dimension(:) ::   Arms_rhos
    real(WP), allocatable, dimension(:) ::  dArms_rhos
    real(WP), allocatable, dimension(:) :: ygArms_rhos

  end type largo_BL_temporal_tconsistent_workspace_type

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


  public  :: largo_BL_temporal_tconsistent_allocate
  public  :: largo_BL_temporal_tconsistent_deallocate
  public  :: largo_BL_temporal_tconsistent_init
  public  :: largo_BL_temporal_tconsistent_preStep_sEta
  public  :: largo_BL_temporal_tconsistent_preStep_sEta_innery
  public  :: largo_BL_temporal_tconsistent_preStep_sEta_innerxz
  public  :: largo_BL_temporal_tconsistent_preStep_sEtaMean
  public  :: largo_BL_temporal_tconsistent_continuity_sEtaMean
  public  :: largo_BL_temporal_tconsistent_xMomentum_sEtaMean
  public  :: largo_BL_temporal_tconsistent_yMomentum_sEtaMean
  public  :: largo_BL_temporal_tconsistent_zMomentum_sEtaMean
  public  :: largo_BL_temporal_tconsistent_energy_sEtaMean
  public  :: largo_BL_temporal_tconsistent_ispecies_sEtaMean
  public  :: largo_BL_temporal_tconsistent_species_sEtaMean
  public  :: largo_BL_temporal_tconsistent_continuity_sEta_
  public  :: largo_BL_temporal_tconsistent_xMomentum_sEta_
  public  :: largo_BL_temporal_tconsistent_yMomentum_sEta_
  public  :: largo_BL_temporal_tconsistent_zMomentum_sEta_
  public  :: largo_BL_temporal_tconsistent_energy_sEta_
  public  :: largo_BL_temporal_tconsistent_ispecies_sEta_
  public  :: largo_BL_temporal_tconsistent_species_sEta_
  public  :: largo_BL_temporal_tconsistent_sEta

  private :: largo_BL_temporal_tconsistent_preStep_sEta_

  public  :: largo_BL_temporal_tconsistent_preStep_baseflow

contains

  subroutine largo_BL_temporal_tconsistent_allocate(cp, neq, ns)

    type(largo_workspace_ptr), intent(out) :: cp
    ! neq=number of equations, might be needed later
    integer(c_int), intent(in)   :: neq
    integer(c_int), intent(in)   :: ns    ! number of species
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp

    ! Allocate derived type variable
    allocate(auxp)

    ! Initialize values
    neq_ = neq
    ns_  = ns

    ! Allocate arrays for species
    if (ns_ > 0) then
      allocate(auxp%dts_rhos    (1:ns_))
      allocate(auxp%ygrms_rhos  (1:ns_))
      allocate(auxp%fluc_rhos   (1:ns_))
      allocate(auxp%mean_rhos   (1:ns_))
      allocate(auxp%field_rhos  (1:ns_))
      allocate(auxp%dtsArms_rhos(1:ns_))
      allocate(auxp%Arms_rhos   (1:ns_))
      allocate(auxp%dArms_rhos  (1:ns_))
      allocate(auxp%ygArms_rhos (1:ns_))
    end if

    ! Get C pointer from Fortran pointer
    cp = c_loc(auxp)

  end subroutine largo_BL_temporal_tconsistent_allocate


  subroutine largo_BL_temporal_tconsistent_deallocate(cp)

    type(largo_workspace_ptr), intent(inout)  :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp

    call c_f_pointer(cp, auxp)

    ! Deallocate arrays for species
    if (allocated(auxp%dts_rhos))     deallocate(auxp%dts_rhos    )
    if (allocated(auxp%ygrms_rhos))   deallocate(auxp%ygrms_rhos  )
    if (allocated(auxp%fluc_rhos))    deallocate(auxp%fluc_rhos   )
    if (allocated(auxp%mean_rhos))    deallocate(auxp%mean_rhos   )
    if (allocated(auxp%field_rhos))   deallocate(auxp%field_rhos  )
    if (allocated(auxp%dtsArms_rhos)) deallocate(auxp%dtsArms_rhos)
    if (allocated(auxp%Arms_rhos))    deallocate(auxp%Arms_rhos   )
    if (allocated(auxp%dArms_rhos))   deallocate(auxp%dArms_rhos  )
    if (allocated(auxp%ygArms_rhos))  deallocate(auxp%ygArms_rhos )

    ! Deallocate array of derived types
    deallocate(auxp)

    ! Nullify C pointer
    cp = c_null_ptr

  end subroutine largo_BL_temporal_tconsistent_deallocate


  subroutine largo_BL_temporal_tconsistent_init(cp, gr_delta, gr_DA, gr_DA_rms)

    real(WP), intent(in)                  :: gr_delta
    real(WP), dimension(*), intent(in)    :: gr_DA
    real(WP), dimension(*), intent(in)    :: gr_DA_rms
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer   :: auxp

    ! get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! set growth rates
    auxp%gr_delta   = gr_delta

!!$     auxp%gr_DA_rho  = gr_DA(irho )
!!$     auxp%gr_DA_rhoU = gr_DA(irhoU)
!!$     auxp%gr_DA_rhoV = gr_DA(irhoV)
!!$     auxp%gr_DA_rhoW = gr_DA(irhoW)
!!$     auxp%gr_DA_rhoE = gr_DA(irhoE)
!!$
!!$     do is=1, ns_
!!$       auxp%gr_DA_rhos(is) = gr_DA(5+is)
!!$     end do

  end subroutine largo_BL_temporal_tconsistent_init


  subroutine largo_BL_temporal_tconsistent_preStep_baseflow(cp, base, &
                              ddy_base, ddt_base, ddx_base, src_base)

    real(WP), dimension(*), intent(in)    :: base
    real(WP), dimension(*), intent(in)    :: ddy_base
    real(WP), dimension(*), intent(in)    :: ddt_base
    real(WP), dimension(*), intent(in)    :: ddx_base
    real(WP), dimension(*), intent(in)    :: src_base
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer   :: auxp

!!$     ! Get Fortran pointer from C pointer
!!$     call c_f_pointer(cp, auxp)
!!$
!!$     ! Store baseflow information
!!$     auxp%base_rho  = base(irho )
!!$     auxp%base_rhoU = base(irhoU)
!!$     auxp%base_rhoV = base(irhoV)
!!$     auxp%base_rhoW = base(irhoW)
!!$     auxp%base_rhoE = base(irhoE)
!!$
!!$     auxp%ddy_base_rho  = ddy_base(irho )
!!$     auxp%ddy_base_rhoU = ddy_base(irhoU)
!!$     auxp%ddy_base_rhoV = ddy_base(irhoV)
!!$     auxp%ddy_base_rhoW = ddy_base(irhoW)
!!$     auxp%ddy_base_rhoE = ddy_base(irhoE)
!!$
!!$     auxp%ddt_base_rho  = ddt_base(irho )
!!$     auxp%ddt_base_rhoU = ddt_base(irhoU)
!!$     auxp%ddt_base_rhoV = ddt_base(irhoV)
!!$     auxp%ddt_base_rhoW = ddt_base(irhoW)
!!$     auxp%ddt_base_rhoE = ddt_base(irhoE)
!!$
!!$     auxp%src_base_rho  = src_base(irho )
!!$     auxp%src_base_rhoU = src_base(irhoU)
!!$     auxp%src_base_rhoV = src_base(irhoV)
!!$     auxp%src_base_rhoW = src_base(irhoW)
!!$     auxp%src_base_rhoE = src_base(irhoE)
!!$
!!$     do is=1, ns_
!!$       auxp%base_rhos    (is) =     base(5+is)
!!$       auxp%ddy_base_rhos(is) = ddy_base(5+is)
!!$       auxp%ddt_base_rhos(is) = ddt_base(5+is)
!!$       auxp%src_base_rhos(is) = src_base(5+is)
!!$     end do

  end subroutine largo_BL_temporal_tconsistent_preStep_baseflow


  subroutine largo_BL_temporal_tconsistent_preStep_sEtaMean(cp, y, mean, ddy_mean)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: ddy_mean
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    ! Compute/get mean and Favre averages
    auxp%mean_rho = mean(irho)
    auxp%fav_U    = mean(irhoU)/mean(irho)
    auxp%fav_V    = mean(irhoV)/mean(irho)
    auxp%fav_W    = mean(irhoW)/mean(irho)
    auxp%fav_E    = mean(irhoE)/mean(irho)
    do is=1, ns_
      auxp%mean_rhos(is) = mean(5+is)
    end do

    ! Compute/get y-derivative of Favre averages
    auxp%dmean_rho = ddy_mean(irho)
    auxp%dfav_U    = ddy_mean(irhoU)/mean(irho) &
      &             - mean(irhoU)/mean(irho)**2 * ddy_mean(irho)
    auxp%dfav_V    = ddy_mean(irhoV)/mean(irho) &
      &             - mean(irhoV)/mean(irho)**2 * ddy_mean(irho)
    auxp%dfav_W    = ddy_mean(irhoW)/mean(irho) &
      &             - mean(irhoW)/mean(irho)**2 * ddy_mean(irho)
    auxp%dfav_E    = ddy_mean(irhoE)/mean(irho) &
      &             - mean(irhoE)/mean(irho)**2 * ddy_mean(irho)

    ! These ones depend on y only
    auxp%dts_rho = y * auxp%gr_delta * ddy_mean(irho )
    auxp%dts_U   = y * auxp%gr_delta * auxp%dfav_U
    auxp%dts_V   = y * auxp%gr_delta * auxp%dfav_V
    auxp%dts_W   = y * auxp%gr_delta * auxp%dfav_W
    auxp%dts_E   = y * auxp%gr_delta * auxp%dfav_E

    ! These ones depend on y only
    auxp%dts_rhoU = y * auxp%gr_delta * ddy_mean(irhoU)
    auxp%dts_rhoV = y * auxp%gr_delta * ddy_mean(irhoV)
    auxp%dts_rhoW = y * auxp%gr_delta * ddy_mean(irhoW)
    auxp%dts_rhoE = y * auxp%gr_delta * ddy_mean(irhoE)

    do is=1, ns_
      auxp%dts_rhos(is)  = y * auxp%gr_delta * ddy_mean(5+is)
    end do

  end subroutine largo_BL_temporal_tconsistent_preStep_sEtaMean


  subroutine largo_BL_temporal_tconsistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: dmean_rqq
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in)          :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp

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
    auxp%Arms_rho = auxp%mean_rho

    ! Compute Arms_U=sqrt{2*tke}
    auxp%Arms_U   = sqrt(auxp%rhoupup + auxp%rhovpvp + auxp%rhowpwp)
    auxp%Arms_V   = auxp%Arms_U
    auxp%Arms_W   = auxp%Arms_U
    auxp%Arms_E   = sqrt(auxp%rhoEpEp)

    ! Assign dArms_rho
    auxp%dArms_rho = auxp%dmean_rho

    ! Compute d(Armsu)/dy=d(sqrt{2*tke})/dy
    auxp%dArms_U  = 0.0_WP
    auxp%dArms_V  = 0.0_WP
    auxp%dArms_W  = 0.0_WP

    if (auxp%Arms_U > eps) then
      auxp%dArms_U  = 0.5_WP / auxp%Arms_U &
        &          * (auxp%drhoupup + auxp%drhovpvp + auxp%drhowpwp)
      auxp%dArms_V  = auxp%dArms_U
      auxp%dArms_W  = auxp%dArms_U
    end if

    auxp%dArms_E  = 0.0_WP
    if (auxp%Arms_E > eps) auxp%dArms_E  = 0.5_WP / auxp%Arms_E * auxp%drhoEpEp

    do is=1, ns_
      auxp%Arms_rhos(is)  = auxp%mean_rho
      auxp%dArms_rhos(is) = auxp%dmean_rho
    end do

    ! These ones depend on y only
    auxp%ygArms_rho  = 0.0_WP
    if (auxp%Arms_rho > eps) auxp%ygArms_rho = y * auxp%gr_delta * auxp%dArms_rho/auxp%Arms_rho

    auxp%ygArms_U = 0.0_WP
    if (auxp%Arms_U   > eps) auxp%ygArms_U   = y * auxp%gr_delta * auxp%dArms_U/auxp%Arms_U

    auxp%ygArms_V = 0.0_WP
    if (auxp%Arms_V   > eps) auxp%ygArms_V   = y * auxp%gr_delta * auxp%dArms_V/auxp%Arms_V

    auxp%ygArms_W = 0.0_WP
    if (auxp%Arms_W   > eps) auxp%ygArms_W   = y * auxp%gr_delta * auxp%dArms_W/auxp%Arms_W

    auxp%ygArms_E = 0.0_WP
    if (auxp%Arms_E   > eps) auxp%ygArms_E   = y * auxp%gr_delta * auxp%dArms_E/auxp%Arms_E

    do is=1, ns_
      auxp%ygArms_rhos(is) = 0.0_WP
      if (auxp%Arms_rhos(is) > eps)  auxp%ygArms_rhos(is) = y * auxp%gr_delta * auxp%dArms_rhos(is)/auxp%Arms_rhos(is)
    end do

  end subroutine largo_BL_temporal_tconsistent_preStep_sEta_


  subroutine largo_BL_temporal_tconsistent_preStep_sEta_innerxz(cp, qflow)

    real(WP), dimension(*), intent(in) :: qflow
    integer(c_int) :: is
    type(largo_workspace_ptr), intent(in) :: cp
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp

    ! Get Fortran pointer from C pointer
    call c_f_pointer(cp, auxp)

    auxp%field_rho = qflow(irho)
    auxp%field_U   = qflow(irhoU)/qflow(irho)
    auxp%field_V   = qflow(irhoV)/qflow(irho)
    auxp%field_W   = qflow(irhoW)/qflow(irho)
    auxp%field_E   = qflow(irhoE)/qflow(irho)

    auxp%fluc_rho = auxp%field_rho - auxp%mean_rho
    auxp%ffluc_U  = auxp%field_U   - auxp%fav_U
    auxp%ffluc_V  = auxp%field_V   - auxp%fav_V
    auxp%ffluc_W  = auxp%field_W   - auxp%fav_W
    auxp%ffluc_E  = auxp%field_E   - auxp%fav_E

    auxp%dtsArms_rho = auxp%fluc_rho * auxp%ygArms_rho
    auxp%dtsArms_U   = auxp%ffluc_U  * auxp%ygArms_U
    auxp%dtsArms_V   = auxp%ffluc_V  * auxp%ygArms_V
    auxp%dtsArms_W   = auxp%ffluc_W  * auxp%ygArms_W
    auxp%dtsArms_E   = auxp%ffluc_E  * auxp%ygArms_E

    do is=1, ns_
      auxp%field_rhos(is)   = qflow(5+is)
      auxp%fluc_rhos(is)    = auxp%field_rhos(is) - auxp%mean_rhos(is)
      auxp%dtsArms_rhos(is) = auxp%fluc_rhos(is) * auxp%ygArms_rhos(is)
    end do

    ! Compute mean plus fluctuations slow time derivative
    auxp%dtsFull_rho = auxp%dts_rho + auxp%dtsArms_rho
    auxp%dtsFull_U   = auxp%dts_U   + auxp%dtsArms_U
    auxp%dtsFull_V   = auxp%dts_V   + auxp%dtsArms_V
    auxp%dtsFull_W   = auxp%dts_W   + auxp%dtsArms_W
    auxp%dtsFull_E   = auxp%dts_E   + auxp%dtsArms_E


!!$     auxp%fluc_rho  = qflow(irho ) - mean(irho )
!!$     auxp%fluc_rhoU = qflow(irhoU) - mean(irhoU)
!!$     auxp%fluc_rhoV = qflow(irhoV) - mean(irhoV)
!!$     auxp%fluc_rhoW = qflow(irhoW) - mean(irhoW)
!!$     auxp%fluc_rhoE = qflow(irhoE) - mean(irhoE)
!!$
!!$     auxp%dtsRms_rho  = auxp%fluc_rho  * auxp%ygrms_rho
!!$     auxp%dtsRms_rhoU = auxp%fluc_rhoU * auxp%ygrms_rhoU
!!$     auxp%dtsRms_rhoV = auxp%fluc_rhoV * auxp%ygrms_rhoV
!!$     auxp%dtsRms_rhoW = auxp%fluc_rhoW * auxp%ygrms_rhoW
!!$     auxp%dtsRms_rhoE = auxp%fluc_rhoE * auxp%ygrms_rhoE
!!$
!!$     do is=1, ns_
!!$       auxp%fluc_rhos(is)  = qflow(5+is) - mean(5+is)
!!$       auxp%dtsRms_rhos(is) = auxp%fluc_rhos(is) * auxp%ygrms_rhos(is)
!!$     end do

  end subroutine largo_BL_temporal_tconsistent_preStep_sEta_innerxz


  subroutine largo_BL_temporal_tconsistent_preStep_sEta_innery(cp, y, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

    real(WP), intent(in)               :: y
    real(WP), dimension(*), intent(in) :: mean
    real(WP), dimension(*), intent(in) :: rms
    real(WP), dimension(*), intent(in) :: mean_rqq
    real(WP), dimension(*), intent(in) :: ddy_mean
    real(WP), dimension(*), intent(in) :: ddy_rms
    real(WP), dimension(*), intent(in) :: dmean_rqq
    type(largo_workspace_ptr), intent(in)        :: cp

    call largo_BL_temporal_tconsistent_preStep_sEtaMean(cp, y, mean, ddy_mean)
    call largo_BL_temporal_tconsistent_preStep_sEta_(cp, y, rms, mean_rqq, ddy_rms, dmean_rqq)

  end subroutine largo_BL_temporal_tconsistent_preStep_sEta_innery


  subroutine largo_BL_temporal_tconsistent_preStep_sEta(cp, y, qflow, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)

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
    call largo_BL_temporal_tconsistent_preStep_sEta_innery(cp, y, mean, rms, mean_rqq, ddy_mean, ddy_rms, dmean_rqq)
    call largo_BL_temporal_tconsistent_preStep_sEta_innerxz(cp, qflow)

  end subroutine largo_BL_temporal_tconsistent_preStep_sEta


#define DECLARE_SUBROUTINE(token)token (cp, A, B, src);\
  type(largo_workspace_ptr), intent(in)   :: cp;\
  real(WP)       , intent(in)             :: A, B;\
  real(WP)       , intent(inout)          :: src;\
  type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp;\
  call c_f_pointer(cp, auxp)\

  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_continuity_sEtaMean)
    src = A * src + B * auxp%dts_rho
  end subroutine largo_BL_temporal_tconsistent_continuity_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_xMomentum_sEtaMean)
    src = A * src + B * auxp%dts_rhoU
  end subroutine largo_BL_temporal_tconsistent_xMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_yMomentum_sEtaMean)
    src = A * src + B * auxp%dts_rhoV
  end subroutine largo_BL_temporal_tconsistent_yMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_zMomentum_sEtaMean)
    src = A * src + B * auxp%dts_rhoW
  end subroutine largo_BL_temporal_tconsistent_zMomentum_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_energy_sEtaMean)
    src = A * src + B * auxp%dts_rhoE
  end subroutine largo_BL_temporal_tconsistent_energy_sEtaMean


  subroutine largo_BL_temporal_tconsistent_ispecies_sEtaMean (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    call c_f_pointer(cp, auxp)
    src = A * src + B * auxp%dts_rhos(is)
  end subroutine largo_BL_temporal_tconsistent_ispecies_sEtaMean


  subroutine largo_BL_temporal_tconsistent_species_sEtaMean (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = ns_
    integer(c_int)                            :: is

    do is = 1, ns_
      call largo_BL_temporal_tconsistent_ispecies_sEtaMean (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_tconsistent_species_sEtaMean


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_continuity_sEta_)
    src = A * src + B * auxp%dtsFull_rho
  end subroutine largo_BL_temporal_tconsistent_continuity_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_xMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%dtsFull_U &
      &                  + auxp%field_U   * auxp%dtsFull_rho  )
  end subroutine largo_BL_temporal_tconsistent_xMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_yMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%dtsFull_V &
      &                  + auxp%field_V   * auxp%dtsFull_rho  )
  end subroutine largo_BL_temporal_tconsistent_yMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_zMomentum_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%dtsFull_W &
      &                  + auxp%field_W   * auxp%dtsFull_rho  )
  end subroutine largo_BL_temporal_tconsistent_zMomentum_sEta_


  subroutine DECLARE_SUBROUTINE(largo_BL_temporal_tconsistent_energy_sEta_)
    src = A * src + B * (  auxp%field_rho * auxp%dtsFull_E &
      &                  + auxp%field_E   * auxp%dtsFull_rho  )
  end subroutine largo_BL_temporal_tconsistent_energy_sEta_


  subroutine largo_BL_temporal_tconsistent_ispecies_sEta_ (cp, A, B, src, is)
    type(largo_workspace_ptr), intent(in)    :: cp
    real(WP)       , intent(in)              :: A, B
    real(WP)       , intent(inout)           :: src
    type(largo_BL_temporal_tconsistent_workspace_type), pointer :: auxp
    integer(c_int), intent(in)               :: is

    ! FIXME: method not implemented
    call c_f_pointer(cp, auxp)
    src = A * src + B * ( auxp%dts_rhos(is) + auxp%dtsArms_rhos(is) )
  end subroutine largo_BL_temporal_tconsistent_ispecies_sEta_


  subroutine largo_BL_temporal_tconsistent_species_sEta_ (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)   :: cp
    real(WP)       , intent(in)             :: A, B
    real(WP), dimension(*), intent(inout)   :: srcvec ! "*" = ns_
    integer(c_int)                          :: is

    ! FIXME: method not implemented
    do is = 1, ns_
      call largo_BL_temporal_tconsistent_ispecies_sEta_ (cp, A, B, srcvec(is), is)
    end do
  end subroutine largo_BL_temporal_tconsistent_species_sEta_


  subroutine largo_BL_temporal_tconsistent_sEta (cp, A, B, srcvec)
    type(largo_workspace_ptr), intent(in)     :: cp
    real(WP)       , intent(in)               :: A, B
    real(WP), dimension(*), intent(inout)     :: srcvec ! "*" = 5+ns_
    integer(c_int)                            :: is

    call largo_BL_temporal_tconsistent_continuity_sEta_ (cp, A, B, srcvec(irho ))
    call largo_BL_temporal_tconsistent_xMomentum_sEta_  (cp, A, B, srcvec(irhou))
    call largo_BL_temporal_tconsistent_yMomentum_sEta_  (cp, A, B, srcvec(irhov))
    call largo_BL_temporal_tconsistent_zMomentum_sEta_  (cp, A, B, srcvec(irhow))
    call largo_BL_temporal_tconsistent_energy_sEta_     (cp, A, B, srcvec(irhoE))

    do is = 1, ns_
      call largo_BL_temporal_tconsistent_ispecies_sEta_ (cp, A, B, srcvec(5+is), is)
    end do

  end subroutine largo_BL_temporal_tconsistent_sEta

end module largo_BL_temporal_tconsistent
