!==============================================================================
! mhdunits : convert normalized MHD quantities back to physical units (common to all tools)
!
!   Optionally write a &units group in the namelist file of each tool (default: Jupiter):
!
!     &units
!       phys_units = .true.,   ! .false. keeps normalized values (same as the legacy tools)
!       planet = 'jupiter',    ! 'jupiter' / 'saturn' / 'earth' / 'custom'
!       ! giving only the 3 basic quantities (and ion mass) computes the rest automatically:
!       ns = 1.0e10, bs = 4.2e-4, xs = 7.14e7, ion_mass = 1.0,
!       ! derived quantities given here take precedence (to match the code constants exactly):
!       ! rhos = 1.67e-17, vs = 9.167e7, ts = 0.779, ps = 0.1403, temps = 8.76e7,
!       length_unit = 'RJ',    ! distance unit of plots (default: Jupiter RJ, Saturn RS, Earth RE)
!       temp_unit = 'eV',      ! 'eV', 'keV', 'K'
!     /
!
!   basic: ns [m^-3] number density, bs [T] equatorial surface field, xs [m] planetary radius, ion_mass [amu]
!   derived (unless given):
!     rhos = ns ion_mass mp,  vs = bs/sqrt(mu0 rhos),  ts = xs/vs,
!     ps = bs^2/mu0 (= rhos vs^2),  temps [eV] = ps/(ns kB)/11600   (T = p/rho * temps)
!   planet defaults (normalization constants of the Jupiter / Saturn / Earth MHD codes):
!     jupiter: ns=1e10, bs=4.2e-4, xs=7.14e7,  rhos=1.67e-17, vs=9.167e7, ts=0.779, ps=0.1403
!     saturn : ns=1e10, bs=2.0e-5, xs=6.004e7, rhos=1.67e-17, vs=4.359e6, ts=13.77, ps=3.185e-4
!     earth  : ns=1e10, bs=3.12e-5, xs=6.378e6, rhos=1.67e-17, vs=6.80e6, ts=0.937, ps=7.75e-4
!     temps = ps/(ns kB)/11600 (Jupiter: the given value 8.76e7 eV)
!   Giving basic quantities (ns, bs, xs, ion_mass) in the namelist recomputes the derived ones.
!
!   kinds of quantities and display units
!     'n'     number density rho*ns    [cm^-3]
!     'rhom'  mass density rho*rhos    [kg m^-3]
!     'v'     velocity                 [km/s]
!     'p'     pressure                 [nPa]
!     'b'     magnetic field           [nT]
!     'temp'  temperature p/rho        [eV] (temp_unit)
!     'pdyn'  dynamic pressure rho v^2 (rhos vs^2) [nPa]
!     'beta'  plasma beta              [-] (2 mu0 p / B^2 in physical units)
!     'j'     current density bs/(mu0 xs) [uA/m^2]
!     'e'     electric field vs bs     [mV/m]
!     'phi'   potential vs bs xs       [kV]
!     'omega' vorticity / angular velocity 1/ts [1/s]
!     'eflux' p sqrt(p/rho)            [mW/m^2]
!     'len'   length                   [length_unit] (not converted)
!     ''      dimensionless / unknown  (not converted)
!
!   All computations are done in double precision (real64); only results are stored in single.
!==============================================================================
module mhdunits
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
  public :: units_t, UC, units_read, ufac, ulabel, units_report, kind_of_name, ulen, default_header

  real(real64), parameter :: mu0 = 4.0e-7_real64*3.14159265358979323846_real64
  real(real64), parameter :: amu = 1.67e-27_real64        ! proton mass (same digits as the MHD code)
  real(real64), parameter :: kb  = 1.3806e-23_real64
  real(real64), parameter :: ev_k = 11600.0_real64        ! 1 eV = 11600 K (same as the MHD code)

  type :: units_t
    logical :: phys = .true.
    character(len=16) :: planet = 'jupiter'
    real(real64) :: ns = 1.0e10_real64, bs = 4.2e-4_real64, xs = 7.14e7_real64, &
                    ion_mass = 1.0_real64, &
                    rhos = 1.67e-17_real64, vs = 9.167e7_real64, ts = 0.779_real64, &
                    ps = 0.1403_real64, temps = 8.76e7_real64
    character(len=8) :: temp_unit = 'eV'
    character(len=8) :: length_unit = 'RJ'
    character(len=64) :: derived = ''        ! quantities computed automatically (for the log)
  end type units_t

  type(units_t), save :: UC                 ! settings shared by all tools (set by units_read)

contains

  !----------------------------------------------------------------------------
  ! read &units from the namelist file (defaults if absent). unit: an opened unit number
  !----------------------------------------------------------------------------
  subroutine units_read(unit)
    integer, intent(in) :: unit
    logical :: phys_units, base_changed
    character(len=16) :: planet
    real(real64) :: ns, bs, xs, ion_mass, rhos, vs, ts, ps, temps
    character(len=8) :: temp_unit, length_unit
    integer :: ios
    character(len=256) :: msg
    namelist /units/ phys_units, planet, ns, bs, xs, ion_mass, rhos, vs, ts, ps, temps, &
                     temp_unit, length_unit
    phys_units = .true. ; planet = 'jupiter' ; temp_unit = 'eV' ; length_unit = ''
    ns = -1 ; bs = -1 ; xs = -1 ; ion_mass = -1
    rhos = -1 ; vs = -1 ; ts = -1 ; ps = -1 ; temps = -1
    rewind(unit)
    read(unit, nml=units, iostat=ios, iomsg=msg)
    if (ios > 0) then
      write(*,'(a)') ' ERROR: namelist &units : '//trim(msg)
      error stop 1
    end if
    rewind(unit)
    UC = units_t()
    UC%phys = phys_units
    UC%planet = lower(planet)
    ! ---- defaults of each planet
    select case (trim(UC%planet))
    case ('jupiter', '')
      UC%planet = 'jupiter' ; UC%length_unit = 'RJ'      ! type defaults (Jupiter code constants)
    case ('saturn')                                     ! Saturn code constants
      UC%ns = 1.0e10_real64 ; UC%bs = 2.0e-5_real64 ; UC%xs = 6.004e7_real64
      UC%rhos = 1.67e-17_real64 ; UC%vs = 4.359e6_real64 ; UC%ts = 13.77_real64
      UC%ps = 3.185e-4_real64 ; UC%temps = UC%ps/(UC%ns*kb)/ev_k
      UC%length_unit = 'RS'
    case ('earth')                                      ! Earth code constants
      UC%ns = 1.0e10_real64 ; UC%bs = 3.12e-5_real64 ; UC%xs = 6.378e6_real64
      UC%rhos = 1.67e-17_real64 ; UC%vs = 6.80e6_real64 ; UC%ts = 0.937_real64
      UC%ps = 7.75e-4_real64 ; UC%temps = UC%ps/(UC%ns*kb)/ev_k
      UC%length_unit = 'RE'
    case default
      if (any([ns, bs, xs] <= 0.0_real64)) then
        write(*,'(a)') ' ERROR: &units : planet = '//trim(planet)// &
                       ' has no preset; give at least ns, bs, xs (and ion_mass, length_unit)'
        error stop 1
      end if
      UC%length_unit = 'R'
    end select
    ! ---- basic quantities given -> recompute derived ones (derived values given are applied afterwards)
    base_changed = .false.
    if (ns       > 0.0_real64) then ; UC%ns = ns             ; base_changed = .true. ; end if
    if (bs       > 0.0_real64) then ; UC%bs = bs             ; base_changed = .true. ; end if
    if (xs       > 0.0_real64) then ; UC%xs = xs             ; base_changed = .true. ; end if
    if (ion_mass > 0.0_real64) then ; UC%ion_mass = ion_mass ; base_changed = .true. ; end if
    if (base_changed) call derive()
    UC%derived = ''
    if (base_changed) then
      if (rhos  <= 0.0_real64) UC%derived = trim(UC%derived)//' rhos'
      if (vs    <= 0.0_real64) UC%derived = trim(UC%derived)//' vs'
      if (ts    <= 0.0_real64) UC%derived = trim(UC%derived)//' ts'
      if (ps    <= 0.0_real64) UC%derived = trim(UC%derived)//' ps'
      if (temps <= 0.0_real64) UC%derived = trim(UC%derived)//' temps'
    end if
    if (rhos  > 0.0_real64) UC%rhos  = rhos
    if (vs    > 0.0_real64) UC%vs    = vs
    if (ts    > 0.0_real64) UC%ts    = ts
    if (ps    > 0.0_real64) UC%ps    = ps
    if (temps > 0.0_real64) UC%temps = temps
    if (len_trim(length_unit) > 0) UC%length_unit = length_unit
    UC%temp_unit = temp_unit
    if (all(UC%temp_unit /= ['eV ', 'keV', 'K  '])) then
      write(*,'(a)') ' ERROR: &units : temp_unit must be eV, keV or K'
      error stop 1
    end if
  contains
    subroutine derive()
      UC%rhos  = UC%ns*UC%ion_mass*amu
      UC%vs    = UC%bs/sqrt(mu0*UC%rhos)
      UC%ts    = UC%xs/UC%vs
      UC%ps    = UC%bs**2/mu0
      UC%temps = UC%ps/(UC%ns*kb)/ev_k
    end subroutine derive
  end subroutine units_read

  ! plot header (when header in the namelist is empty)
  function default_header() result(s)
    character(len=128) :: s
    select case (trim(UC%planet))
    case ('jupiter') ; s = 'High resolution Jovian Magnetosphere'
    case ('saturn')  ; s = 'High resolution Kronian Magnetosphere'
    case ('earth')   ; s = 'High resolution Earth Magnetosphere'
    case default     ; s = 'MHD simulation of the Magnetosphere'
    end select
  end function default_header

  ! distance unit (for axis labels)
  function ulen() result(s)
    character(len=8) :: s
    s = UC%length_unit
  end function ulen

  !----------------------------------------------------------------------------
  ! factor normalized -> display unit. 1 if phys_units=.false. (except 2 of beta = 2p/B^2)
  !----------------------------------------------------------------------------
  pure real(real64) function ufac(kind) result(f)
    character(len=*), intent(in) :: kind
    f = 1.0_real64
    if (.not. UC%phys) then
      if (trim(kind) == 'beta') f = 2.0_real64
      return
    end if
    select case (trim(kind))
    case ('n')     ; f = UC%ns*1.0e-6_real64                       ! cm^-3
    case ('rhom')  ; f = UC%rhos                                   ! kg m^-3
    case ('v')     ; f = UC%vs*1.0e-3_real64                       ! km/s
    case ('p')     ; f = UC%ps*1.0e9_real64                        ! nPa
    case ('b')     ; f = UC%bs*1.0e9_real64                        ! nT
    case ('temp')
      select case (trim(UC%temp_unit))
      case ('keV') ; f = UC%temps*1.0e-3_real64
      case ('K')   ; f = UC%temps*11600.0_real64
      case default ; f = UC%temps
      end select
    case ('pdyn')  ; f = UC%rhos*UC%vs**2*1.0e9_real64              ! nPa
    case ('beta')  ; f = 2.0_real64*mu0*UC%ps/UC%bs**2              ! 2 mu0 p / B^2
    case ('j')     ; f = UC%bs/(mu0*UC%xs)*1.0e6_real64             ! uA/m^2
    case ('e')     ; f = UC%vs*UC%bs*1.0e3_real64                   ! mV/m
    case ('phi')   ; f = UC%vs*UC%bs*UC%xs*1.0e-3_real64             ! kV
    case ('omega') ; f = 1.0_real64/UC%ts                          ! 1/s
    case ('eflux') ; f = UC%ps*sqrt(UC%ps/UC%rhos)*1.0e3_real64      ! mW/m^2
    case default   ; f = 1.0_real64
    end select
  end function ufac

  ! display unit string ('norm.' if phys_units=.false.)
  function ulabel(kind) result(s)
    character(len=*), intent(in) :: kind
    character(len=12) :: s
    s = ''
    if (len_trim(kind) == 0 .or. trim(kind) == 'beta' .or. trim(kind) == 'len') then
      if (trim(kind) == 'len') s = UC%length_unit
      return
    end if
    if (.not. UC%phys) then
      s = 'norm.' ; return
    end if
    select case (trim(kind))
    case ('n')     ; s = 'cm^-3'
    case ('rhom')  ; s = 'kg/m^3'
    case ('v')     ; s = 'km/s'
    case ('p', 'pdyn') ; s = 'nPa'
    case ('b')     ; s = 'nT'
    case ('temp')  ; s = UC%temp_unit
    case ('j')     ; s = 'uA/m^2'
    case ('e')     ; s = 'mV/m'
    case ('phi')   ; s = 'kV'
    case ('omega') ; s = '1/s'
    case ('eflux') ; s = 'mW/m^2'
    end select
  end function ulabel

  ! variable / quantity name -> kind ('' if unknown)
  function kind_of_name(name) result(k)
    character(len=*), intent(in) :: name
    character(len=8) :: k
    select case (lower(trim(name)))
    case ('rho', 'n', 'dens', 'density') ; k = 'n'
    case ('vx', 'vy', 'vz', 'absv', 'v') ; k = 'v'
    case ('p', 'pres', 'pressure')       ; k = 'p'
    case ('bx', 'by', 'bz', 'absb', 'absbz', 'b') ; k = 'b'
    case ('temp', 't')                   ; k = 'temp'
    case ('pdyn')                        ; k = 'pdyn'
    case ('beta')                        ; k = 'beta'
    case ('par', 'perp', 'wx', 'wy', 'wz', 'wmag', 'paravort', 'perpvort') ; k = 'omega'
    case default                         ; k = ''
    end select
  end function kind_of_name

  subroutine units_report()
    if (.not. UC%phys) then
      write(*,'(a)') ' units: normalized (phys_units = .false.)'
      return
    end if
    write(*,'(a,a,a,a)') ' units: physical, planet = ', trim(UC%planet), '   length unit = ', &
         trim(UC%length_unit)
    write(*,'(a,es11.4,a,es11.4,a,es11.4,a,f6.2,a,es11.4)') '   ns=', UC%ns, ' m^-3  bs=', UC%bs, &
         ' T  xs=', UC%xs, ' m  ion_mass=', UC%ion_mass, '  rhos=', UC%rhos
    write(*,'(a,es11.4,a,es11.4,a,es11.4,a,es11.4,a)') '   vs=', UC%vs, ' m/s  ts=', UC%ts, &
         ' s  ps=', UC%ps, ' Pa  temps=', UC%temps, ' eV'
    if (len_trim(UC%derived) > 0) write(*,'(a)') '   derived from ns, bs, xs, ion_mass:'//trim(UC%derived)
    write(*,'(a)') '   display units: n [cm^-3], v [km/s], p [nPa], B [nT], T ['// &
         trim(UC%temp_unit)//'], J [uA/m^2], phi [kV], omega [1/s]'
    ! consistency of the normalization (ps should equal rhos vs^2 = bs^2/mu0)
    if (abs(UC%ps/(UC%rhos*UC%vs**2) - 1.0_real64) > 0.05_real64) &
      write(*,'(a,es11.4,a,es11.4,a)') '   NOTE: ps =', UC%ps, ' Pa differs from rhos*vs^2 =', &
         UC%rhos*UC%vs**2, ' Pa (check the pressure normalisation; beta uses 2 mu0 p/B^2)'
  end subroutine units_report

  pure function lower(s) result(o)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i
    o = s
    do i = 1, len(s)
      if (iachar(s(i:i)) >= 65 .and. iachar(s(i:i)) <= 90) o(i:i) = achar(iachar(s(i:i))+32)
    end do
  end function lower

end module mhdunits
