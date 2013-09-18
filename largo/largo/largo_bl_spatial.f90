module largo_BL_spatial

  implicit none ! is in effect throughout entire module

  private
  integer, parameter :: WP = 8

  type largo_BL_spatial_type
    real(WP) :: uu         = 0.0_WP
    real(WP) :: vv         = 0.0_WP
    real(WP) :: ww         = 0.0_WP
    real(WP) :: pp         = 0.0_WP
    real(WP) :: hh         = 0.0_WP

    real(WP) :: dxs_rho    = 0.0_WP
    real(WP) :: dxs_rhoU   = 0.0_WP
    real(WP) :: dxs_rhoV   = 0.0_WP
    real(WP) :: dxs_rhoW   = 0.0_WP
    real(WP) :: dxs_rhoE   = 0.0_WP
    real(WP) :: dxsMean_p  = 0.0_WP

    real(WP) :: grms_rho   = 0.0_WP
    real(WP) :: grms_rhoU  = 0.0_WP
    real(WP) :: grms_rhoV  = 0.0_WP
    real(WP) :: grms_rhoW  = 0.0_WP
    real(WP) :: grms_rhoE  = 0.0_WP
    real(WP) :: grms_p     = 0.0_WP

    real(WP) :: dxsRms_rho   = 0.0_WP
    real(WP) :: dxsRms_rhoU  = 0.0_WP
    real(WP) :: dxsRms_rhoV  = 0.0_WP
    real(WP) :: dxsRms_rhoW  = 0.0_WP
    real(WP) :: dxsRms_rhoE  = 0.0_WP
    real(WP) :: dxsRms_p     = 0.0_WP

  end type largo_BL_spatial_type

  integer, parameter :: irho  = 1
  integer, parameter :: irhoU = 2
  integer, parameter :: irhoV = 3
  integer, parameter :: irhoW = 4
  integer, parameter :: irhoE = 5
  integer, parameter :: ip    = 6
  real(WP), parameter :: eps = 1.0E-10_WP

  public :: largo_BL_spatial_type
!!$   public :: largo_BL_spatial_preStep_sEtaMean
  public :: largo_BL_spatial_preStep_sEta
  public :: largo_BL_spatial_continuity_sEtaMean
  public :: largo_BL_spatial_xMomentum_sEtaMean
  public :: largo_BL_spatial_yMomentum_sEtaMean
  public :: largo_BL_spatial_zMomentum_sEtaMean
  public :: largo_BL_spatial_energy_sEtaMean
  public :: largo_BL_spatial_continuity_sEtaRms
  public :: largo_BL_spatial_xMomentum_sEtaRms
  public :: largo_BL_spatial_yMomentum_sEtaRms
  public :: largo_BL_spatial_zMomentum_sEtaRms
  public :: largo_BL_spatial_energy_sEtaRms

contains

  subroutine largo_BL_spatial_preStep_sEta(y, qflow, mean, rms, ddy_mean, ddy_rms, aux)

    real(WP), intent(in) :: y
!!$     real(WP), intent(in) :: gamma
    real(WP), dimension(1:6), intent(in)  :: qflow
    real(WP), dimension(1:6), intent(in)  :: mean
    real(WP), dimension(1:6), intent(in)  :: rms
    real(WP), dimension(1:6), intent(in)  :: ddy_mean
    real(WP), dimension(1:6), intent(in)  :: ddy_rms
    type(largo_BL_spatial_type), intent(out) :: aux

    real(WP) :: inverse_rho
    real(WP) :: fluc_rho   = 0.0_WP
    real(WP) :: fluc_rhoU  = 0.0_WP
    real(WP) :: fluc_rhoV  = 0.0_WP
    real(WP) :: fluc_rhoW  = 0.0_WP
    real(WP) :: fluc_rhoE  = 0.0_WP
    real(WP) :: fluc_p     = 0.0_WP

    ! These ones depend on y only
    aux%dxs_rho   = y * ddy_mean(irho )
    aux%dxs_rhoU  = y * ddy_mean(irhoU)
    aux%dxs_rhoV  = y * ddy_mean(irhoV)
    aux%dxs_rhoW  = y * ddy_mean(irhoW)
    aux%dxs_rhoE  = y * ddy_mean(irhoE)
    aux%dxsMean_p = y * ddy_mean(ip)

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

    aux%grms_p    = 0.0_WP
    if (rms(ip)    > eps)  aux%grms_p    = ddy_rms(ip   )/(rms(ip   ))

    ! These ones depend on x, y, z
    inverse_rho = 1.0_WP/qflow(irho)
    aux%uu = qflow(irhoU) * inverse_rho
    aux%vv = qflow(irhoV) * inverse_rho
    aux%ww = qflow(irhoW) * inverse_rho
!!$     aux%pp = (gamma-1.0_WP) * (mean(irhoE) &
!!$       &  - 0.5_WP * qflow(irho)*(aux%uu**2+aux%vv**2+aux%ww**2))
    aux%pp = qflow(ip   )
    aux%hh = (qflow(irhoE) + aux%pp) * inverse_rho

    fluc_rho  = qflow(irho ) - mean(irho )
    fluc_rhoU = qflow(irhoU) - mean(irhoU)
    fluc_rhoV = qflow(irhoV) - mean(irhoV)
    fluc_rhoW = qflow(irhoW) - mean(irhoW)
    fluc_rhoE = qflow(irhoE) - mean(irhoE)
    fluc_p    = qflow(ip)    - mean(ip)

    aux%dxsRms_rho  = y * fluc_rho  * aux%grms_rho
    aux%dxsRms_rhoU = y * fluc_rhoU * aux%grms_rhoU
    aux%dxsRms_rhoV = y * fluc_rhoV * aux%grms_rhoV
    aux%dxsRms_rhoW = y * fluc_rhoW * aux%grms_rhoW
    aux%dxsRms_rhoE = y * fluc_rhoE * aux%grms_rhoE
    aux%dxsRms_p    = y * fluc_p    * aux%grms_p

!!$     aux%dxsMean_p = (gamma-1.0_WP) * (aux%dxs_rhoE &
!!$       &  - aux%uu * aux%dxs_rhoU &
!!$       &  - aux%vv * aux%dxs_rhoV &
!!$       &  - aux%ww * aux%dxs_rhoW &
!!$       &  + 0.5_WP * (aux%uu**2 + aux%vv**2 + aux%ww**2) * aux%dxs_rho)
!!$
!!$     aux%dxsRms_p  = (gamma-1.0_WP) * (aux%dxsRms_rhoE &
!!$       &  - aux%uu * aux%dxsRms_rhoU &
!!$       &  - aux%vv * aux%dxsRms_rhoV &
!!$       &  - aux%ww * aux%dxsRms_rhoW &
!!$       &  + 0.5_WP * (aux%uu**2 + aux%vv**2 + aux%ww**2) * aux%dxsRms_rho)

  end subroutine largo_BL_spatial_preStep_sEta


  subroutine largo_BL_spatial_sEta(aux, A, B, src)

    type(largo_BL_spatial_type), intent(in) :: aux
    real(WP), intent(in) :: A, B
    real(WP), intent(inout) :: src

    entry largo_BL_spatial_continuity_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dxs_rhoU
    return

    entry largo_BL_spatial_xMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (2.0_wp * aux%uu * aux%dxs_rhoU - aux%uu**2 * aux%dxs_rho + aux%dxsMean_p)
    return

    entry largo_BL_spatial_yMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%uu * aux%dxs_rhoV + aux%vv * aux%dxs_rhoU - aux%uu*aux%vv * aux%dxs_rho )
    return

    entry largo_BL_spatial_zMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%uu * aux%dxs_rhoW + aux%ww * aux%dxs_rhoU - aux%uu*aux%ww * aux%dxs_rho )
    return

    entry largo_BL_spatial_energy_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%uu * (aux%dxs_rhoE+aux%dxsMean_p) + aux%hh*(aux%dxs_rhoU - aux%uu*aux%dxs_rho))
    return

    entry largo_BL_spatial_continuity_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoU
    return

    entry largo_BL_spatial_xMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (2.0_wp * aux%uu * aux%dxsRms_rhoU - aux%uu**2 * aux%dxsRms_rho + aux%dxsRms_p)
    return

    entry largo_BL_spatial_yMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (aux%uu * aux%dxsRms_rhoV + aux%vv * aux%dxsRms_rhoU - aux%uu*aux%vv * aux%dxsRms_rho )
    return

    entry largo_BL_spatial_zMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (aux%uu * aux%dxsRms_rhoW + aux%ww * aux%dxsRms_rhoU - aux%uu*aux%ww * aux%dxsRms_rho )
    return

    entry largo_BL_spatial_energy_sEtaRms(aux, A, B, src)
    src = A * src + B * (aux%uu * (aux%dxsRms_rhoE+aux%dxsRms_p) + aux%hh*(aux%dxsRms_rhoU - aux%uu*aux%dxsRms_rho))
    return

  end subroutine largo_BL_spatial_sEta

end module largo_BL_spatial
