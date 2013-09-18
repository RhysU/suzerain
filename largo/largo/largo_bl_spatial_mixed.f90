module largo_BL_spatial_mixed

  implicit none ! is in effect throughout entire module

  private
  integer, parameter :: WP = 8

  type largo_BL_spatial_mixed_type
    real(WP) :: ufav       = 0.0_WP
    real(WP) :: vfav       = 0.0_WP
    real(WP) :: wfav       = 0.0_WP
    real(WP) :: pbar       = 0.0_WP
    real(WP) :: hfav       = 0.0_WP

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

    real(WP) :: dxsRms_rho   = 0.0_WP
    real(WP) :: dxsRms_rhoU  = 0.0_WP
    real(WP) :: dxsRms_rhoV  = 0.0_WP
    real(WP) :: dxsRms_rhoW  = 0.0_WP
    real(WP) :: dxsRms_rhoE  = 0.0_WP

    real(WP) :: dx_over_dt   = 0.0_WP
  end type largo_BL_spatial_mixed_type

  integer, parameter :: irho  = 1
  integer, parameter :: irhoU = 2
  integer, parameter :: irhoV = 3
  integer, parameter :: irhoW = 4
  integer, parameter :: irhoE = 5
  integer, parameter :: ip    = 6
  real(WP), parameter :: eps = 1.0E-10_WP

  public :: largo_BL_spatial_mixed_type
!!$   public :: largo_BL_spatial_mixed_preStep_sEtaMean
  public :: largo_BL_spatial_mixed_preStep_sEta
  public :: largo_BL_spatial_mixed_continuity_sEtaMean
  public :: largo_BL_spatial_mixed_xMomentum_sEtaMean
  public :: largo_BL_spatial_mixed_yMomentum_sEtaMean
  public :: largo_BL_spatial_mixed_zMomentum_sEtaMean
  public :: largo_BL_spatial_mixed_energy_sEtaMean
  public :: largo_BL_spatial_mixed_continuity_sEtaRms
  public :: largo_BL_spatial_mixed_xMomentum_sEtaRms
  public :: largo_BL_spatial_mixed_yMomentum_sEtaRms
  public :: largo_BL_spatial_mixed_zMomentum_sEtaRms
  public :: largo_BL_spatial_mixed_energy_sEtaRms

contains

  subroutine largo_BL_spatial_mixed_preStep_sEta(y, qflow, mean, rms, ddy_mean, ddy_rms, dx_over_dt, aux)

    real(WP), intent(in) :: y
!!$     real(WP), intent(in) :: gamma
    real(WP), dimension(1:5), intent(in)  :: qflow
    real(WP), dimension(1:6), intent(in)  :: mean
    real(WP), dimension(1:5), intent(in)  :: rms
    real(WP), dimension(1:6), intent(in)  :: ddy_mean
    real(WP), dimension(1:5), intent(in)  :: ddy_rms
    real(WP), intent(in) :: dx_over_dt
    type(largo_BL_spatial_mixed_type), intent(out) :: aux

    real(WP) :: inverse_rho
    real(WP) :: fluc_rho   = 0.0_WP
    real(WP) :: fluc_rhoU  = 0.0_WP
    real(WP) :: fluc_rhoV  = 0.0_WP
    real(WP) :: fluc_rhoW  = 0.0_WP
    real(WP) :: fluc_rhoE  = 0.0_WP

    ! These ones depend on y only
    aux%dxs_rho   = y * ddy_mean(irho )
    aux%dxs_rhoU  = y * ddy_mean(irhoU)
    aux%dxs_rhoV  = y * ddy_mean(irhoV)
    aux%dxs_rhoW  = y * ddy_mean(irhoW)
    aux%dxs_rhoE  = y * ddy_mean(irhoE)
    aux%dxsMean_p = y * ddy_mean(ip   )

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

    inverse_rho = 1.0_WP/mean(irho)
    aux%ufav = mean(irhoU) * inverse_rho
    aux%vfav = mean(irhoV) * inverse_rho
    aux%wfav = mean(irhoW) * inverse_rho
    aux%pbar = mean(ip   )
!!$     aux%pbar = (gamma-1.0_WP) * (mean(irhoE) &
!!$       &  - 0.5_WP * mean(irho)*(aux%ufav**2+aux%vfav**2+aux%wfav**2))
    aux%hfav = (mean(irhoE) + aux%pbar) * inverse_rho

    ! These ones depend on x, y, z
    fluc_rho  = qflow(irho ) - mean(irho )
    fluc_rhoU = qflow(irhoU) - mean(irhoU)
    fluc_rhoV = qflow(irhoV) - mean(irhoV)
    fluc_rhoW = qflow(irhoW) - mean(irhoW)
    fluc_rhoE = qflow(irhoE) - mean(irhoE)

    aux%dxsRms_rho  = y * dx_over_dt * fluc_rho  * aux%grms_rho
    aux%dxsRms_rhoU = y * dx_over_dt * fluc_rhoU * aux%grms_rhoU
    aux%dxsRms_rhoV = y * dx_over_dt * fluc_rhoV * aux%grms_rhoV
    aux%dxsRms_rhoW = y * dx_over_dt * fluc_rhoW * aux%grms_rhoW
    aux%dxsRms_rhoE = y * dx_over_dt * fluc_rhoE * aux%grms_rhoE

!!$     aux%dxsMean_p = (gamma-1.0_WP) * (aux%dxs_rhoE &
!!$       &  - aux%ufav * aux%dxs_rhoU &
!!$       &  - aux%vfav * aux%dxs_rhoV &
!!$       &  - aux%wfav * aux%dxs_rhoW &
!!$       &  + 0.5_WP * (aux%ufav**2 + aux%vfav**2 + aux%wfav**2) * aux%dxs_rho)

  end subroutine largo_BL_spatial_mixed_preStep_sEta


  subroutine largo_BL_spatial_mixed_sEta(aux, A, B, src)

    type(largo_BL_spatial_mixed_type), intent(in) :: aux
    real(WP), intent(in) :: A, B
    real(WP), intent(inout) :: src

    entry largo_BL_spatial_mixed_continuity_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dxs_rhoU
    return

    entry largo_BL_spatial_mixed_xMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (2.0_wp * aux%ufav * aux%dxs_rhoU - aux%ufav**2 * aux%dxs_rho + aux%dxsMean_p)
    return

    entry largo_BL_spatial_mixed_yMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%ufav * aux%dxs_rhoV + aux%vfav * aux%dxs_rhoU - aux%ufav*aux%vfav * aux%dxs_rho )
    return

    entry largo_BL_spatial_mixed_zMomentum_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%ufav * aux%dxs_rhoW + aux%wfav * aux%dxs_rhoU - aux%ufav*aux%wfav * aux%dxs_rho )
    return

    entry largo_BL_spatial_mixed_energy_sEtaMean(aux, A, B, src)
    src = A * src + B * (aux%ufav * (aux%dxs_rhoE+aux%dxsMean_p) + aux%hfav*(aux%dxs_rhoU - aux%ufav*aux%dxs_rho))
    return

    entry largo_BL_spatial_mixed_continuity_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rho
    return

    entry largo_BL_spatial_mixed_xMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoU
    return

    entry largo_BL_spatial_mixed_yMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoV
    return

    entry largo_BL_spatial_mixed_zMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoW
    return

    entry largo_BL_spatial_mixed_energy_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoE
    return

  end subroutine largo_BL_spatial_mixed_sEta

end module largo_BL_spatial_mixed
