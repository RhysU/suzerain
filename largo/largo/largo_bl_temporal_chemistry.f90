module largo_BL_temporal_chemistry

  implicit none ! is in effect throughout entire module

  private
  integer, parameter :: WP = 8

  type largo_BL_temporal_chemistry_type

    real(WP) :: dts_rho    = 0.0_WP
    real(WP) :: dts_rhoU   = 0.0_WP
    real(WP) :: dts_rhoV   = 0.0_WP
    real(WP) :: dts_rhoW   = 0.0_WP
    real(WP) :: dts_rhoE   = 0.0_WP

    real(WP) :: grms_rho   = 0.0_WP
    real(WP) :: grms_rhoU  = 0.0_WP
    real(WP) :: grms_rhoV  = 0.0_WP
    real(WP) :: grms_rhoW  = 0.0_WP
    real(WP) :: grms_rhoE  = 0.0_WP

    real(WP) :: dtsRms_rho   = 0.0_WP
    real(WP) :: dtsRms_rhoU  = 0.0_WP
    real(WP) :: dtsRms_rhoV  = 0.0_WP
    real(WP) :: dtsRms_rhoW  = 0.0_WP
    real(WP) :: dtsRms_rhoE  = 0.0_WP

  end type largo_BL_temporal_chemistry_type

  type largo_BL_temporal_species_type

    real(WP) :: dts_rho_s    = 0.0_WP
    real(WP) :: grms_rho_s   = 0.0_WP
    real(WP) :: dtsRms_rho_s   = 0.0_WP

  end type largo_BL_temporal_species_type


  integer, parameter :: irho  = 1
  integer, parameter :: irhoU = 2
  integer, parameter :: irhoV = 3
  integer, parameter :: irhoW = 4
  integer, parameter :: irhoE = 5
  real(WP), parameter :: eps = 1.0E-10_WP

  public :: largo_BL_temporal_chemistry_type
  public :: largo_BL_temporal_species_type
!!$   public :: largo_BL_temporal_chemistry_preStep_sEtaMean
  public :: largo_BL_temporal_chemistry_preStep_sEta
  public :: largo_BL_temporal_chemistry_species_preStep_sEta

  public :: largo_BL_temporal_chemistry_continuity_sEtaMean
  public :: largo_BL_temporal_chemistry_xMomentum_sEtaMean
  public :: largo_BL_temporal_chemistry_yMomentum_sEtaMean
  public :: largo_BL_temporal_chemistry_zMomentum_sEtaMean
  public :: largo_BL_temporal_chemistry_energy_sEtaMean
  public :: largo_BL_temporal_chemistry_species_sEtaMean
  public :: largo_BL_temporal_chemistry_continuity_sEtaRms
  public :: largo_BL_temporal_chemistry_xMomentum_sEtaRms
  public :: largo_BL_temporal_chemistry_yMomentum_sEtaRms
  public :: largo_BL_temporal_chemistry_zMomentum_sEtaRms
  public :: largo_BL_temporal_chemistry_energy_sEtaRms
  public :: largo_BL_temporal_chemistry_species_sEtaRms

contains

  subroutine largo_BL_temporal_chemistry_preStep_sEta(y, qflow, mean, rms, ddy_mean, ddy_rms, aux)

    real(WP), intent(in) :: y

    real(WP), dimension(1:5), intent(in)  :: qflow
    real(WP), dimension(1:5), intent(in)  :: mean
    real(WP), dimension(1:5), intent(in)  :: rms
    real(WP), dimension(1:5), intent(in)  :: ddy_mean
    real(WP), dimension(1:5), intent(in)  :: ddy_rms
    type(largo_BL_temporal_chemistry_type), intent(out) :: aux

    !integer :: is
!!$     real(WP) :: inverse_rho
    real(WP) :: fluc_rho   = 0.0_WP
    real(WP) :: fluc_rhoU  = 0.0_WP
    real(WP) :: fluc_rhoV  = 0.0_WP
    real(WP) :: fluc_rhoW  = 0.0_WP
    real(WP) :: fluc_rhoE  = 0.0_WP

    !real(WP) :: fluc_rho_s = 0.0_WP

    ! These ones depend on y only
    aux%dts_rho  = y * ddy_mean(irho )
    aux%dts_rhoU = y * ddy_mean(irhoU)
    aux%dts_rhoV = y * ddy_mean(irhoV)
    aux%dts_rhoW = y * ddy_mean(irhoW)
    aux%dts_rhoE = y * ddy_mean(irhoE)

    aux%grms_rho  = 0.0_WP
    if (rms(irho ) > eps)  aux%grms_rho  = ddy_rms(irho )/(rms(irho ))

    aux%grms_rhoU = 0.0_WP
    if (rms(irhoU) > eps)  aux%grms_rhoU = ddy_rms(irhoU)/(rms(irhoU))

    aux%grms_rhoV = 0.0_WP
    if (rms(irhoV) > eps)  aux%grms_rhoV = ddy_rms(irhoV)/(rms(irhoV))

    aux%grms_rhoW = 0.0_WP
    if (rms(irhoW) > eps)  aux%grms_rhoW = ddy_rms(irhoW)/(rms(irhoW))

    aux%grms_rhoE = 0.0_WP
    if (rms(irhoE) > eps)  aux%grms_rhoE = ddy_rms(irhoE)/(rms(irhoE))

    fluc_rho  = qflow(irho ) - mean(irho )
    fluc_rhoU = qflow(irhoU) - mean(irhoU)
    fluc_rhoV = qflow(irhoV) - mean(irhoV)
    fluc_rhoW = qflow(irhoW) - mean(irhoW)
    fluc_rhoE = qflow(irhoE) - mean(irhoE)

    aux%dtsRms_rho  = y * fluc_rho  * aux%grms_rho
    aux%dtsRms_rhoU = y * fluc_rhoU * aux%grms_rhoU
    aux%dtsRms_rhoV = y * fluc_rhoV * aux%grms_rhoV
    aux%dtsRms_rhoW = y * fluc_rhoW * aux%grms_rhoW
    aux%dtsRms_rhoE = y * fluc_rhoE * aux%grms_rhoE

  end subroutine largo_BL_temporal_chemistry_preStep_sEta


  subroutine largo_BL_temporal_chemistry_species_preStep_sEta(y, ns, qflow, mean, rms, ddy_mean, ddy_rms, aux_s)

    real(WP), intent(in) :: y
    integer, intent(in) :: ns

    real(WP), dimension(1:ns), intent(in)  :: qflow
    real(WP), dimension(1:ns), intent(in)  :: mean
    real(WP), dimension(1:ns), intent(in)  :: rms
    real(WP), dimension(1:ns), intent(in)  :: ddy_mean
    real(WP), dimension(1:ns), intent(in)  :: ddy_rms
    type(largo_BL_temporal_species_type), intent(out) :: aux_s(1:ns)

    integer :: is
    real(WP) :: fluc_rho_s = 0.0_WP

    do is=1, ns
      aux_s(is)%dts_rho_s  = y * ddy_mean(is)

      aux_s(is)%grms_rho_s = 0.0_WP
      if (rms(is) > eps)  aux_s(is)%grms_rho_s = ddy_rms(is)/(rms(is))

      fluc_rho_s = qflow(is) - mean(is)
      aux_s(is)%dtsRms_rho_s = y * fluc_rho_s * aux_s(is)%grms_rho_s
    end do

  end subroutine largo_BL_temporal_chemistry_species_preStep_sEta


  subroutine largo_BL_temporal_chemistry_sEta(aux, A, B, src)

    type(largo_BL_temporal_chemistry_type), intent(in) :: aux
    real(WP), intent(in) :: A, B
    real(WP), intent(inout) :: src

    !integer :: is

    entry largo_BL_temporal_chemistry_continuity_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dts_rho
    return

    entry largo_BL_temporal_chemistry_xMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dts_rhoU
    return

    entry largo_BL_temporal_chemistry_yMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dts_rhoV
    return

    entry largo_BL_temporal_chemistry_zMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dts_rhoW
    return

    entry largo_BL_temporal_chemistry_energy_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dts_rhoE
    return

    entry largo_BL_temporal_chemistry_continuity_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dtsRms_rho
    return

    entry largo_BL_temporal_chemistry_xMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dtsRms_rhoU
    return

    entry largo_BL_temporal_chemistry_yMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dtsRms_rhoV
    return

    entry largo_BL_temporal_chemistry_zMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dtsRms_rhoW
    return

    entry largo_BL_temporal_chemistry_energy_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dtsRms_rhoE
    return

  end subroutine largo_BL_temporal_chemistry_sEta


  subroutine largo_BL_temporal_chemistry_species_sEta(aux_s, A, B, src)

    type(largo_BL_temporal_species_type), intent(in) :: aux_s
    real(WP), intent(in) :: A, B
    real(WP), intent(inout) :: src

    !integer :: is

    entry largo_BL_temporal_chemistry_species_sEtaMean(aux_s, A, B, src)
    src = A * src + B * aux_s%dts_rho_s
    return

    entry largo_BL_temporal_chemistry_species_sEtaRms(aux_s, A, B, src)
    src = A * src + B * aux_s%dtsRms_rho_s
    return

  end subroutine largo_BL_temporal_chemistry_species_sEta


end module largo_BL_temporal_chemistry
