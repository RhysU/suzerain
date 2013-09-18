module largo_BL_spatial_first_order

  implicit none ! is in effect throughout entire module

  private
  integer, parameter :: WP = 8

  type largo_BL_spatial_first_order_type
    real(WP) :: fluc_rho   = 0.0_WP
    real(WP) :: fluc_rhoU  = 0.0_WP
    real(WP) :: fluc_rhoV  = 0.0_WP
    real(WP) :: fluc_rhoW  = 0.0_WP
    real(WP) :: fluc_rhoE  = 0.0_WP
    real(WP) :: fluc_p     = 0.0_WP

    real(WP) :: meanRho    = 0.0_WP
    real(WP) :: ufav       = 0.0_WP
    real(WP) :: vfav       = 0.0_WP
    real(WP) :: wfav       = 0.0_WP
    real(WP) :: Efav       = 0.0_WP
    real(WP) :: pbar       = 0.0_WP
    real(WP) :: pinvrho    = 0.0_WP

    real(WP) :: dxs_rho    = 0.0_WP
    real(WP) :: dxs_rhoU   = 0.0_WP
    real(WP) :: dxs_rhoV   = 0.0_WP
    real(WP) :: dxs_rhoW   = 0.0_WP
    real(WP) :: dxs_rhoE   = 0.0_WP
    real(WP) :: dxs_p      = 0.0_WP

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

    real(WP) :: dxs_ufav   = 0.0_WP
    real(WP) :: dxs_vfav   = 0.0_WP
    real(WP) :: dxs_wfav   = 0.0_WP
    real(WP) :: dxs_Efav   = 0.0_WP
    real(WP) :: dxs_pinvrho  = 0.0_WP

  end type largo_BL_spatial_first_order_type

  integer, parameter :: irho  = 1
  integer, parameter :: irhoU = 2
  integer, parameter :: irhoV = 3
  integer, parameter :: irhoW = 4
  integer, parameter :: irhoE = 5
  integer, parameter :: ip    = 6
  real(WP), parameter :: eps = 1.0E-10_WP

  public :: largo_BL_spatial_first_order_type
  public :: largo_BL_spatial_first_order_preStep_sEta
  public :: largo_BL_spatial_first_order_continuity_sEtaMean
  public :: largo_BL_spatial_first_order_xMomentum_sEtaMean
  public :: largo_BL_spatial_first_order_yMomentum_sEtaMean
  public :: largo_BL_spatial_first_order_zMomentum_sEtaMean
  public :: largo_BL_spatial_first_order_energy_sEtaMean
  public :: largo_BL_spatial_first_order_continuity_sEtaRms
  public :: largo_BL_spatial_first_order_xMomentum_sEtaRms
  public :: largo_BL_spatial_first_order_yMomentum_sEtaRms
  public :: largo_BL_spatial_first_order_zMomentum_sEtaRms
  public :: largo_BL_spatial_first_order_energy_sEtaRms

contains

  subroutine largo_BL_spatial_first_order_preStep_sEta(y, qflow, mean, rms, ddy_mean, ddy_rms, aux)

    real(WP), intent(in) :: y
    real(WP), dimension(1:6), intent(in)  :: qflow
    real(WP), dimension(1:6), intent(in)  :: mean
    real(WP), dimension(1:6), intent(in)  :: rms
    real(WP), dimension(1:6), intent(in)  :: ddy_mean
    real(WP), dimension(1:6), intent(in)  :: ddy_rms
    type(largo_BL_spatial_first_order_type), intent(out) :: aux

    real(WP) :: inverse_rho

    ! These ones depend on y only
    ! Mean
    aux%dxs_rho  = y * ddy_mean(irho )
    aux%dxs_rhoU = y * ddy_mean(irhoU)
    aux%dxs_rhoV = y * ddy_mean(irhoV)
    aux%dxs_rhoW = y * ddy_mean(irhoW)
    aux%dxs_rhoE = y * ddy_mean(irhoE)
    aux%dxs_p    = y * ddy_mean(ip)

    ! Fluctuations
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

    inverse_rho = 1.0_WP/mean(irho)
    aux%meanrho = mean(irho)
    aux%ufav = mean(irhoU) * inverse_rho
    aux%vfav = mean(irhoV) * inverse_rho
    aux%wfav = mean(irhoW) * inverse_rho
    aux%pbar  = mean(ip)
    aux%Efav  = mean(irhoE) * inverse_rho

    !! pinvrho = pbar/meanrho
    aux%pinvrho = aux%pbar * inverse_rho

    ! Favre
    aux%dxs_ufav  = (aux%dxs_rhoU - aux%ufav    * aux%dxs_rho) * inverse_rho
    aux%dxs_vfav  = (aux%dxs_rhoV - aux%vfav    * aux%dxs_rho) * inverse_rho
    aux%dxs_wfav  = (aux%dxs_rhoW - aux%wfav    * aux%dxs_rho) * inverse_rho
    aux%dxs_Efav  = (aux%dxs_rhoE - aux%Efav    * aux%dxs_rho) * inverse_rho
    aux%dxs_pinvrho = (aux%dxs_p  - aux%pinvrho * aux%dxs_rho) * inverse_rho

    ! These ones depend on x, y, z
    aux%fluc_rho  = qflow(irho ) - mean(irho )
    aux%fluc_rhoU = qflow(irhoU) - mean(irhoU)
    aux%fluc_rhoV = qflow(irhoV) - mean(irhoV)
    aux%fluc_rhoW = qflow(irhoW) - mean(irhoW)
    aux%fluc_rhoE = qflow(irhoE) - mean(irhoE)
    aux%fluc_p    = qflow(ip)    - mean(ip)

    aux%dxsRms_rho  = y * aux%fluc_rho  * aux%grms_rho
    aux%dxsRms_rhoU = y * aux%fluc_rhoU * aux%grms_rhoU
    aux%dxsRms_rhoV = y * aux%fluc_rhoV * aux%grms_rhoV
    aux%dxsRms_rhoW = y * aux%fluc_rhoW * aux%grms_rhoW
    aux%dxsRms_rhoE = y * aux%fluc_rhoE * aux%grms_rhoE
    aux%dxsRms_p    = y * aux%fluc_p    * aux%grms_p

  end subroutine largo_BL_spatial_first_order_preStep_sEta


  subroutine largo_BL_spatial_first_order_sEta(aux, A, B, src)

    type(largo_BL_spatial_first_order_type), intent(in) :: aux
    real(WP), intent(in) :: A, B
    real(WP), intent(inout) :: src

    entry largo_BL_spatial_first_order_continuity_sEtaMean(aux, A, B, src)
    src = A * src + B * aux%dxs_rhoU
    return

    entry largo_BL_spatial_first_order_xMomentum_sEtaMean(aux, A, B, src)
!!$     src = A * src + B * (2.0_wp * aux%meanrho * aux%ufav * aux%dxs_ufav + aux%ufav**2 * aux%dxs_rho + aux%dxs_p)
    src = A * src + B * (2.0_wp * aux%ufav * aux%dxs_rhoU                  - aux%ufav**2 * aux%dxs_rho         + aux%dxs_p)
    return

    entry largo_BL_spatial_first_order_yMomentum_sEtaMean(aux, A, B, src)
!!$     src = A * src + B * (aux%meanrho * aux%vfav * aux%dxs_ufav + aux%meanrho * aux%ufav * aux%dxs_vfav + aux%ufav*aux%vfav * aux%dxs_rho )
    src = A * src + B * (aux%ufav * aux%dxs_rhoV + aux%vfav * aux%dxs_rhoU - aux%ufav*aux%vfav * aux%dxs_rho )
    return

    entry largo_BL_spatial_first_order_zMomentum_sEtaMean(aux, A, B, src)
!!$     src = A * src + B * (aux%meanrho * aux%wfav * aux%dxs_ufav + aux%meanrho * aux%ufav * aux%dxs_wfav + aux%ufav*aux%wfav * aux%dxs_rho )
    src = A * src + B * (aux%ufav * aux%dxs_rhoW + aux%wfav * aux%dxs_rhoU - aux%ufav*aux%wfav * aux%dxs_rho )
    return

    entry largo_BL_spatial_first_order_energy_sEtaMean(aux, A, B, src)
!!$     src = A * src + B * (aux%ufav * (aux%dxs_rhoE+aux%dxsMean_p) + aux%hfav*(aux%dxs_rhoU - aux%ufav*aux%dxs_rho))
    src = A * src + B * (aux%ufav * (aux%dxs_rhoE+aux%dxs_p) + (aux%Efav + aux%pinvrho) *(aux%dxs_rhoU - aux%ufav*aux%dxs_rho))
    return

    entry largo_BL_spatial_first_order_continuity_sEtaRms(aux, A, B, src)
    src = A * src + B * aux%dxsRms_rhoU
    return

    entry largo_BL_spatial_first_order_xMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (  2.0_WP * aux%ufav                 * aux%dxsRms_rhoU &
     &                   + 2.0_WP * aux%fluc_rhoU            * aux%dxs_ufav &
     &                   - 2.0_WP * aux%fluc_rho  * aux%ufav * aux%dxs_ufav &
     &                   -          aux%ufav**2              * aux%dxsRms_rho &
     &                   + aux%dxsRms_p)
    return

    entry largo_BL_spatial_first_order_yMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (  aux%ufav                 * aux%dxsRms_rhoV &
     &                   + aux%fluc_rhoV            * aux%dxs_ufav &
     &                   + aux%vfav                 * aux%dxsRms_rhoU &
     &                   + aux%fluc_rhoU            * aux%dxs_vfav &
     &                   - aux%fluc_rho  * aux%ufav * aux%dxs_vfav &
     &                   - aux%fluc_rho  * aux%vfav * aux%dxs_ufav &
     &                   - aux%ufav      * aux%vfav * aux%dxsRms_rho)
    return

    entry largo_BL_spatial_first_order_zMomentum_sEtaRms(aux, A, B, src)
    src = A * src + B * (  aux%ufav                 * aux%dxsRms_rhoW &
     &                   + aux%fluc_rhoW            * aux%dxs_ufav &
     &                   + aux%wfav                 * aux%dxsRms_rhoU &
     &                   + aux%fluc_rhoU            * aux%dxs_wfav &
     &                   - aux%fluc_rho  * aux%ufav * aux%dxs_wfav &
     &                   - aux%fluc_rho  * aux%wfav * aux%dxs_ufav &
     &                   - aux%ufav      * aux%wfav * aux%dxsRms_rho)
    return

    entry largo_BL_spatial_first_order_energy_sEtaRms(aux, A, B, src)
    ! Energy terms
    src = A * src + B * (  aux%ufav                 * aux%dxsRms_rhoE &
     &                   + aux%fluc_rhoE            * aux%dxs_ufav &
     &                   + aux%Efav                 * aux%dxsRms_rhoU &
     &                   + aux%fluc_rhoU            * aux%dxs_Efav &
     &                   - aux%fluc_rho  * aux%ufav * aux%dxs_Efav &
     &                   - aux%fluc_rho  * aux%Efav * aux%dxs_ufav &
     &                   - aux%ufav      * aux%Efav * aux%dxsRms_rho)
    ! Pressure terms
    src =     src + B * (  aux%ufav                    * aux%dxsRms_p &
     &                   + aux%fluc_p                  * aux%dxs_ufav &
     &                   + aux%pinvrho                 * aux%dxsRms_rhoU &
     &                   + aux%fluc_rhoU               * aux%dxs_pinvrho &
     &                   - aux%fluc_rho  * aux%ufav    * aux%dxs_pinvrho &
     &                   - aux%fluc_rho  * aux%pinvrho * aux%dxs_ufav &
     &                   - aux%ufav      * aux%pinvrho * aux%dxsRms_rho)
    return

  end subroutine largo_BL_spatial_first_order_sEta

end module largo_BL_spatial_first_order
