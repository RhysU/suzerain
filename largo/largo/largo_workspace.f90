module largo_workspace

! Use ISO_C_BINDING to expose C-friendly API
  use, intrinsic :: iso_c_binding, only:   WP => c_double    &
                                         , c_int             &
                                         , largo_ptr =>c_ptr &
                                         , largo_workspace_ptr =>c_ptr


  implicit none ! is in effect throughout entire module

  
  abstract interface
    subroutine allocate_rans(cp, fmodel)
      import
      character(*), intent(in)              :: fmodel
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine allocate_rans
  end interface


  abstract interface
    subroutine init(cp, grDelta, grDA, grDArms)
      import
      real(WP), intent(in)                  :: grDelta
      real(WP), dimension(*), intent(in)    :: grDA
      real(WP), dimension(*), intent(in)    :: grDArms
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine init
  end interface


  abstract interface
    subroutine get_ntvar_rans(cp, ntvar)
      import
      integer(c_int), intent(out)           :: ntvar
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine get_ntvar_rans
  end interface


  abstract interface
    subroutine init_rans(cp, grDAturb)
      import
      real(WP), dimension(*), intent(in)    :: grDAturb
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine init_rans
  end interface


  abstract interface
    subroutine finalize(cp)
      import
      type(largo_workspace_ptr), intent(inout) :: cp
    end subroutine finalize
  end interface


  abstract interface
    subroutine prestep_mean(cp, y, mean, ddy_mean)
      import
      real(WP), intent(in)                  :: y
      real(WP), dimension(*), intent(in)    :: mean
      real(WP), dimension(*), intent(in)    :: ddy_mean
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_mean
  end interface


  abstract interface
    subroutine prestep_innerxz(cp, qflow)
      import
      real(WP), dimension(*), intent(in)    :: qflow
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_innerxz
  end interface


  abstract interface
    subroutine prestep_innery(cp, y,     mean,     rms,     mean_rqq, &
                                     ddy_mean, ddy_rms, ddy_mean_rqq)
      import
      real(WP), intent(in)                  :: y
      real(WP), dimension(*), intent(in)    :: mean
      real(WP), dimension(*), intent(in)    :: rms
      real(WP), dimension(*), intent(in)    :: mean_rqq
      real(WP), dimension(*), intent(in)    :: ddy_mean
      real(WP), dimension(*), intent(in)    :: ddy_rms
      real(WP), dimension(*), intent(in)    :: ddy_mean_rqq
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_innery
  end interface


  abstract interface
    subroutine prestep(cp, y, qflow,     mean,     rms,     mean_rqq, &
                                     ddy_mean, ddy_rms, ddy_mean_rqq)
      import
      real(WP), intent(in)                  :: y
      real(WP), dimension(*), intent(in)    :: qflow
      real(WP), dimension(*), intent(in)    :: mean
      real(WP), dimension(*), intent(in)    :: rms
      real(WP), dimension(*), intent(in)    :: mean_rqq
      real(WP), dimension(*), intent(in)    :: ddy_mean
      real(WP), dimension(*), intent(in)    :: ddy_rms
      real(WP), dimension(*), intent(in)    :: ddy_mean_rqq
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep
  end interface


  abstract interface
    subroutine prestep_baseflow(cp,     base, ddy_base, &
                                    ddt_base, ddx_base, src_base)
      import
      real(WP), dimension(*), intent(in)    ::     base
      real(WP), dimension(*), intent(in)    :: ddy_base
      real(WP), dimension(*), intent(in)    :: ddt_base
      real(WP), dimension(*), intent(in)    :: ddx_base
      real(WP), dimension(*), intent(in)    :: src_base
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_baseflow
  end interface


  abstract interface
    subroutine prestep_setamean_rans(cp, y, mean, ddy_mean)
      import
      real(WP), intent(in)                  :: y
      real(WP), dimension(*), intent(in)    :: mean
      real(WP), dimension(*), intent(in)    :: ddy_mean
      type(largo_workspace_ptr), intent(in) :: cp
    end subroutine prestep_setamean_rans
  end interface


  abstract interface
    subroutine source(cp, A, B, src)
      import
      type(largo_workspace_ptr), intent(in)  :: cp
      real(WP)       , intent(in)            :: A, B
      real(WP)       , intent(inout)         :: src
    end subroutine source
  end interface


  abstract interface
    subroutine isource(cp, A, B, src, is)
      import
      type(largo_workspace_ptr), intent(in)  :: cp
      real(WP)       , intent(in)            :: A, B
      real(WP)       , intent(inout)         :: src
      integer(c_int) , intent(in)            :: is
    end subroutine isource
  end interface


  abstract interface
    subroutine sourcevec(cp, A, B, srcvec)
      import
      type(largo_workspace_ptr), intent(in)  :: cp
      real(WP), intent(in)                   :: A, B
      real(WP), dimension(*), intent(inout)  :: srcvec
    end subroutine sourcevec
  end interface


  ! largo type - for generic interface
  type :: largo_type

    ! Procedure pointers to largo methods
    procedure(init),             pointer, nopass :: largo_init              => NULL()
    procedure(finalize),         pointer, nopass :: largo_finalize          => NULL()
    procedure(prestep_mean),     pointer, nopass :: largo_prestep_mean      => NULL()
    procedure(prestep_innerxz),  pointer, nopass :: largo_prestep_innerxz   => NULL()
    procedure(prestep_innery),   pointer, nopass :: largo_prestep_innery    => NULL()
    procedure(prestep),          pointer, nopass :: largo_prestep           => NULL()
    procedure(source),           pointer, nopass :: largo_continuity_mean   => NULL()
    procedure(source),           pointer, nopass :: largo_xmomentum_mean    => NULL()
    procedure(source),           pointer, nopass :: largo_ymomentum_mean    => NULL()
    procedure(source),           pointer, nopass :: largo_zmomentum_mean    => NULL()
    procedure(source),           pointer, nopass :: largo_energy_mean       => NULL()
    procedure(sourcevec),        pointer, nopass :: largo_species_mean      => NULL()
    procedure(isource),          pointer, nopass :: largo_ispecies_mean     => NULL()
    procedure(source),           pointer, nopass :: largo_continuity_rms    => NULL()
    procedure(source),           pointer, nopass :: largo_xmomentum_rms     => NULL()
    procedure(source),           pointer, nopass :: largo_ymomentum_rms     => NULL()
    procedure(source),           pointer, nopass :: largo_zmomentum_rms     => NULL()
    procedure(source),           pointer, nopass :: largo_energy_rms        => NULL()
    procedure(sourcevec),        pointer, nopass :: largo_species_rms       => NULL()
    procedure(isource),          pointer, nopass :: largo_ispecies_rms      => NULL()

    ! Wrapper procedures
    procedure(source),           pointer, nopass :: largo_continuity        => NULL()
    procedure(source),           pointer, nopass :: largo_xmomentum         => NULL()
    procedure(source),           pointer, nopass :: largo_ymomentum         => NULL()
    procedure(source),           pointer, nopass :: largo_zmomentum         => NULL()
    procedure(source),           pointer, nopass :: largo_energy            => NULL()
    procedure(sourcevec),        pointer, nopass :: largo_species           => NULL()
    procedure(sourcevec),        pointer, nopass :: largo_all_sources_mean  => NULL()
    procedure(sourcevec),        pointer, nopass :: largo_all_sources       => NULL()

    ! Baseflow
    procedure(prestep_baseflow), pointer, nopass :: largo_prestep_baseflow    => NULL()
    procedure(prestep_baseflow), pointer, nopass :: largo_init_wall_baseflow  => NULL()

    ! RANS Turbulence variables
    procedure(get_ntvar_rans),        pointer, nopass :: largo_get_ntvar_rans        => NULL()
    procedure(init_rans),             pointer, nopass :: largo_init_rans             => NULL()
    procedure(prestep_setamean_rans), pointer, nopass :: largo_prestep_setamean_rans => NULL()
    procedure(sourcevec),             pointer, nopass :: largo_sources_mean_rans     => NULL()

    ! Pointer to largo field-related data (workspace)
    type(largo_workspace_ptr) :: cp

    ! Number of equations
    integer(c_int) :: neq = 0

    ! Number of species
    integer(c_int) :: ns  = 0

    ! Number of variables (includes pressure if relevant)
    integer(c_int) :: nvar = 0

    ! Number of turbulence variables
    integer(c_int) :: ntvar = 0

  end type largo_type

contains

end module largo_workspace
