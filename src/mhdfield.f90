!==============================================================================
! mhdfield : partial reading, interpolation and field-line tracing of 3-D MHD data (m.data format)
!
!   * grid_t  : description of the data grid (grid position of the planet center, half domain or not)
!   * block_t : 3-D arrays of variables read over any range with any thinning (stride)
!               full x range with stride=1 is read plane by plane, otherwise row by row
!               reading only the rows needed (about 1/s^2 of the data for stride=s)
!   * sample  : trilinear interpolation; z<0 mirrored with parity for half-domain data
!   * trace   : field-line tracing (2nd-order Runge-Kutta) with point thinning
!   * evaluation of quantity names (rho, vx, ..., temp, absb, beta, log10:<qty> etc.)
!
!   All positions are simulation coordinates centered on the planet [planetary radii] (x: anti-sunward, z: north)
!==============================================================================
module mhdfield
  use mhd_common
  use mhdunits
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: grid_t, block_t, grid_init, load_block, sample, qty_t, qty_parse, &
            qty_eval, qty_vars, qty_needs, qty_loaded, qty_unit, trace_line, nan32, VAR_RHO, VAR_VX, VAR_VY,   &
            VAR_VZ, VAR_P, VAR_BX, VAR_BY, VAR_BZ, PARITY8, field_t, fsample, &
            END_PLANET, END_OUT, END_MAXLEN, END_NULL, END_LOOP

  integer, parameter :: VAR_RHO = 1, VAR_VX = 2, VAR_VY = 3, VAR_VZ = 4, VAR_P = 5, &
                        VAR_BX = 6, VAR_BY = 7, VAR_BZ = 8
  integer, parameter :: PARITY8(8) = [1, 1, 1, -1, 1, -1, -1, 1]
  character(len=3), parameter :: VARNAME8(8) = ['rho', 'vx ', 'vy ', 'vz ', 'p  ', 'bx ', 'by ', 'bz ']
  integer, parameter :: END_PLANET = 1, END_OUT = 2, END_MAXLEN = 3, END_NULL = 4, END_LOOP = 5

  type :: grid_t
    integer :: nx2 = 0, ny2 = 0, nz2 = 0, nvar = 8
    real    :: h = 1.0, c(3) = 0.0          ! grid spacing, grid position of the planet center (i0,j0,k0)
    logical :: half = .true.                ! z>=0 only (symmetric about the equator)
    real    :: xlo(3) = 0.0, xhi(3) = 0.0   ! data domain (includes z<0 for half-domain data)
  end type grid_t

  type :: block_t
    integer :: n(3) = 0, nv = 0
    real    :: org(3) = 0.0, dx = 1.0
    integer, allocatable :: var(:)          ! slot -> variable number in the file
    real(real32), allocatable :: a(:,:,:,:)
  end type block_t

  ! a pair of a fine block (near the planet) and a coarse block (whole domain)
  type :: field_t
    type(grid_t)  :: g
    type(block_t) :: fine, coarse
    logical :: has_fine = .false.
  end type field_t

  ! quantity (result of parsing the name)
  type :: qty_t
    character(len=32) :: name = ''
    integer :: base = 0        ! 1..nvar: a variable itself, <0: derived quantity
    character(len=8) :: kind = ''   ! kind of unit (mhdunits)
    real(real64) :: fac = 1.0_real64 ! normalized -> display unit
    integer :: op(4) = 0       ! operators applied afterwards (1 log10, 2 abs, 3 neg)
    integer :: nop = 0
  end type qty_t
  integer, parameter :: Q_TEMP = -1, Q_ABSB = -2, Q_ABSBZ = -3, Q_ABSV = -4, &
                        Q_BETA = -5, Q_PDYN = -6

contains

  pure function nan32() result(x)
    real(real32) :: x
    x = ieee_value(x, ieee_quiet_nan)
  end function nan32

  !----------------------------------------------------------------------------
  subroutine grid_init(g, nx2, ny2, nz2, nvar, h, i0, j0, k0, nxp, half)
    type(grid_t), intent(out) :: g
    integer, intent(in) :: nx2, ny2, nz2, nvar, nxp
    real, intent(in) :: h, i0, j0, k0
    logical, intent(in) :: half
    g%nx2 = nx2 ; g%ny2 = ny2 ; g%nz2 = nz2 ; g%nvar = nvar ; g%h = h ; g%half = half
    g%c(1) = i0 ; g%c(2) = j0 ; g%c(3) = k0
    if (i0 < 0.0) g%c(1) = real(nx2+1)/2.0 - real(nxp)
    if (j0 < 0.0) g%c(2) = real(ny2+1)/2.0
    if (k0 < 0.0) g%c(3) = merge(1.5, real(nz2+1)/2.0, half)
    g%xlo = [(1.0-g%c(1))*h, (1.0-g%c(2))*h, (1.0-g%c(3))*h]
    g%xhi = [(real(nx2)-g%c(1))*h, (real(ny2)-g%c(2))*h, (real(nz2)-g%c(3))*h]
    if (half) g%xlo(3) = -g%xhi(3)
  end subroutine grid_init

  !----------------------------------------------------------------------------
  ! read variables vars(:) over [lo,hi] with stride s
  !----------------------------------------------------------------------------
  subroutine load_block(fin, g, vars, s, lo, hi, b)
    type(seqfile), intent(in) :: fin
    type(grid_t), intent(in) :: g
    integer, intent(in) :: vars(:), s
    real, intent(in) :: lo(3), hi(3)
    type(block_t), intent(inout) :: b
    integer :: ia(3), ib(3), nn(3), d, m, kk, jj, k, j, nxr
    real(real32), allocatable :: row(:), pl(:,:)
    real :: l, u
    nn = [g%nx2, g%ny2, g%nz2]
    do d = 1, 3
      l = max(lo(d), -1.0e8) ; u = min(hi(d), 1.0e8)
      ia(d) = max(1,     floor(g%c(d) + l/g%h))
      ib(d) = min(nn(d), ceiling(g%c(d) + u/g%h))
    end do
    if (g%half) ia(3) = 1                        ! k=1 (one plane at z<0) is needed for mirroring
    if (any(ib - ia < s)) call fatal('load_block: region too small')
    b%n  = (ib - ia)/s + 1
    b%nv = size(vars)
    b%dx = g%h*real(s)
    b%org = [(real(ia(1))-g%c(1))*g%h, (real(ia(2))-g%c(2))*g%h, (real(ia(3))-g%c(3))*g%h]
    if (allocated(b%a)) deallocate(b%a)
    if (allocated(b%var)) deallocate(b%var)
    allocate(b%a(b%n(1), b%n(2), b%n(3), b%nv), b%var(b%nv))
    b%var = vars
    nxr = ib(1) - ia(1) + 1
    if (s == 1 .and. ia(1) == 1 .and. ib(1) == g%nx2) then
      ! full x range: read the needed j range at once
      allocate(pl(g%nx2, b%n(2)))
      do m = 1, b%nv
        do kk = 1, b%n(3)
          k = ia(3) + kk - 1
          call seq_read_part(fin, rec(vars(m),k), int(ia(2)-1,int64)*g%nx2, pl, &
                             int(g%nx2,int64)*b%n(2))
          b%a(:,:,kk,m) = pl
        end do
      end do
    else
      allocate(row(nxr))
      do m = 1, b%nv
        do kk = 1, b%n(3)
          k = ia(3) + (kk-1)*s
          do jj = 1, b%n(2)
            j = ia(2) + (jj-1)*s
            call seq_read_part(fin, rec(vars(m),k), int(j-1,int64)*g%nx2 + ia(1) - 1, &
                               row, int(nxr,int64))
            b%a(:,jj,kk,m) = row(1:nxr:s)
          end do
        end do
      end do
    end if
  contains
    pure integer(int64) function rec(mv, kz)
      integer, intent(in) :: mv, kz
      rec = int(mv-1,int64)*g%nz2 + kz
    end function rec
  end subroutine load_block

  !----------------------------------------------------------------------------
  ! trilinear interpolation. slots(:): variable slots in the block, par(:): sign for z mirroring
  !----------------------------------------------------------------------------
  logical function sample(b, half, pos, slots, par, val)
    type(block_t), intent(in) :: b
    logical, intent(in) :: half
    real, intent(in) :: pos(3)
    integer, intent(in) :: slots(:), par(:)
    real, intent(out) :: val(:)
    real :: p(3), f(3), w(3), w1(3)
    integer :: i(3), m, s
    logical :: mir
    p = pos ; mir = .false.
    if (half .and. p(3) < 0.0) then
      p(3) = -p(3) ; mir = .true.
    end if
    f = (p - b%org)/b%dx + 1.0
    i = floor(f)
    sample = all(i >= 1) .and. all(i <= b%n - 1)
    if (.not. sample) return
    w = f - real(i) ; w1 = 1.0 - w
    do m = 1, size(slots)
      s = slots(m)
      val(m) = w1(3)*(w1(2)*(w1(1)*b%a(i(1),i(2),i(3),s) + w(1)*b%a(i(1)+1,i(2),i(3),s))      &
                    +  w(2)*(w1(1)*b%a(i(1),i(2)+1,i(3),s) + w(1)*b%a(i(1)+1,i(2)+1,i(3),s)))  &
             +  w(3)*(w1(2)*(w1(1)*b%a(i(1),i(2),i(3)+1,s) + w(1)*b%a(i(1)+1,i(2),i(3)+1,s))  &
                    +  w(2)*(w1(1)*b%a(i(1),i(2)+1,i(3)+1,s) + w(1)*b%a(i(1)+1,i(2)+1,i(3)+1,s)))
      if (mir) val(m) = val(m)*real(par(m))
    end do
  end function sample

  ! look in the fine block, then the coarse one. vars(:) = variable numbers in the file
  logical function fsample(fm, pos, vars, val)
    type(field_t), intent(in) :: fm
    real, intent(in) :: pos(3)
    integer, intent(in) :: vars(:)
    real, intent(out) :: val(:)
    integer :: sl(size(vars)), pr(size(vars)), m
    fsample = .false.
    if (any(pos < fm%g%xlo) .or. any(pos > fm%g%xhi)) return
    do m = 1, size(vars)
      pr(m) = 1
      if (vars(m) <= 8) pr(m) = PARITY8(vars(m))
    end do
    if (fm%has_fine) then
      if (slots_of(fm%fine, vars, sl)) then
        fsample = sample(fm%fine, fm%g%half, pos, sl, pr, val)
        if (fsample) return
      end if
    end if
    if (slots_of(fm%coarse, vars, sl)) fsample = sample(fm%coarse, fm%g%half, pos, sl, pr, val)
  end function fsample

  logical function slots_of(b, vars, sl)
    type(block_t), intent(in) :: b
    integer, intent(in) :: vars(:)
    integer, intent(out) :: sl(:)
    integer :: m, k
    slots_of = .true.
    do m = 1, size(vars)
      sl(m) = 0
      do k = 1, b%nv
        if (b%var(k) == vars(m)) sl(m) = k
      end do
      if (sl(m) == 0) slots_of = .false.
    end do
  end function slots_of

  !============================================================================
  ! physical quantities
  !============================================================================
  subroutine qty_parse(name, q)
    character(len=*), intent(in) :: name
    type(qty_t), intent(out) :: q
    character(len=32) :: s
    integer :: ic
    q%name = name
    s = adjustl(name)
    do
      ic = index(s, ':')
      if (ic == 0) exit
      q%nop = q%nop + 1
      select case (s(1:ic-1))
      case ('log10') ; q%op(q%nop) = 1
      case ('abs')   ; q%op(q%nop) = 2
      case ('neg')   ; q%op(q%nop) = 3
      case default   ; call fatal('unknown operator in '//trim(name))
      end select
      s = s(ic+1:)
    end do
    select case (trim(s))
    case ('rho') ; q%base = VAR_RHO
    case ('vx')  ; q%base = VAR_VX
    case ('vy')  ; q%base = VAR_VY
    case ('vz')  ; q%base = VAR_VZ
    case ('p')   ; q%base = VAR_P
    case ('bx')  ; q%base = VAR_BX
    case ('by')  ; q%base = VAR_BY
    case ('bz')  ; q%base = VAR_BZ
    case ('temp'); q%base = Q_TEMP
    case ('absb'); q%base = Q_ABSB
    case ('absbz'); q%base = Q_ABSBZ
    case ('absv'); q%base = Q_ABSV
    case ('beta'); q%base = Q_BETA
    case ('pdyn'); q%base = Q_PDYN
    case default
      if (s(1:3) == 'var') then
        read(s(4:), *) q%base
      else
        call fatal('unknown quantity '//trim(name))
      end if
    end select
    if (q%base > 0 .and. q%base <= 8) then
      q%kind = kind_of_name(VARNAME8(q%base))
    else if (q%base < 0) then
      q%kind = kind_of_name(s)
    end if
    q%fac = ufac(q%kind)
  end subroutine qty_parse

  ! display unit ('log10 ' is prefixed if an operator such as log10: is applied)
  function qty_unit(q) result(s)
    type(qty_t), intent(in) :: q
    character(len=24) :: s
    s = ulabel(q%kind)
    if (q%nop > 0 .and. len_trim(s) > 0) then
      if (any(q%op(1:q%nop) == 1)) s = 'log10 '//trim(s)
    end if
  end function qty_unit

  ! which of the file variables 1..8 this quantity needs (returned as a function)
  pure function qty_needs(q) result(nd)
    type(qty_t), intent(in) :: q
    logical :: nd(8)
    nd = .false.
    select case (q%base)
    case (Q_TEMP)
      nd(VAR_RHO) = .true. ; nd(VAR_P) = .true.
    case (Q_ABSB)
      nd(VAR_BX) = .true. ; nd(VAR_BY) = .true. ; nd(VAR_BZ) = .true.
    case (Q_ABSBZ)
      nd(VAR_BZ) = .true.
    case (Q_ABSV)
      nd(VAR_VX) = .true. ; nd(VAR_VY) = .true. ; nd(VAR_VZ) = .true.
    case (Q_BETA)
      nd(VAR_P) = .true. ; nd(VAR_BX) = .true. ; nd(VAR_BY) = .true. ; nd(VAR_BZ) = .true.
    case (Q_PDYN)
      nd(VAR_RHO) = .true. ; nd(VAR_VX) = .true. ; nd(VAR_VY) = .true. ; nd(VAR_VZ) = .true.
    case (1:8)
      nd(q%base) = .true.
    end select
  end function qty_needs

  ! for compatibility: add the variables needed to need(:)
  subroutine qty_vars(q, need)
    type(qty_t), intent(in) :: q
    logical, intent(inout) :: need(:)
    logical :: nd(8)
    integer :: m
    nd = qty_needs(q)
    do m = 1, 8
      if (nd(m)) need(m) = .true.
    end do
    if (q%base > 8 .and. q%base <= size(need)) need(q%base) = .true.
  end subroutine qty_vars

  ! are all variables needed to compute q read into block b?
  logical function qty_loaded(q, b)
    type(qty_t), intent(in) :: q
    type(block_t), intent(in) :: b
    logical :: nd(8)
    integer :: m
    nd = qty_needs(q)
    qty_loaded = .true.
    do m = 1, 8
      if (nd(m)) qty_loaded = qty_loaded .and. any(b%var == m)
    end do
    if (q%base > 8) qty_loaded = qty_loaded .and. any(b%var == q%base)
  end function qty_loaded

  ! v(1:nvar): all variables at one point (unused ones ignored). Computed in double precision, converted to display units
  pure real function qty_eval(q, v) result(xo)
    type(qty_t), intent(in) :: q
    real, intent(in) :: v(:)
    real(real64) :: x, d(8)
    integer :: k
    d = 0.0_real64
    d(1:min(8,size(v))) = real(v(1:min(8,size(v))), real64)
    select case (q%base)
    case (Q_TEMP) ; x = d(VAR_P)/max(d(VAR_RHO), 1.0e-30_real64)
    case (Q_ABSB) ; x = sqrt(d(VAR_BX)**2 + d(VAR_BY)**2 + d(VAR_BZ)**2)
    case (Q_ABSBZ); x = abs(d(VAR_BZ))
    case (Q_ABSV) ; x = sqrt(d(VAR_VX)**2 + d(VAR_VY)**2 + d(VAR_VZ)**2)
    case (Q_BETA) ; x = d(VAR_P)/max(d(VAR_BX)**2 + d(VAR_BY)**2 + d(VAR_BZ)**2, 1.0e-200_real64)
    case (Q_PDYN) ; x = d(VAR_RHO)*(d(VAR_VX)**2 + d(VAR_VY)**2 + d(VAR_VZ)**2)
    case default  ; x = real(v(q%base), real64)
    end select
    x = x*q%fac
    do k = q%nop, 1, -1
      select case (q%op(k))
      case (1) ; x = log10(max(x, 1.0e-30_real64))
      case (2) ; x = abs(x)
      case (3) ; x = -x
      end select
    end do
    xo = real(x)
  end function qty_eval

  !============================================================================
  ! field-line tracing from p0 in direction dir (+1: along B, -1: opposite)
  !   returns thinned points pts(3, npt). iend: reason of termination
  !   END_PLANET when inside r_stop. seg: maximum spacing of the points kept
  !============================================================================
  subroutine trace_line(fm, p0, dir, ds, maxlen, r_stop, seg, pts, npt, iend)
    type(field_t), intent(in) :: fm
    real, intent(in) :: p0(3), dir, ds, maxlen, r_stop, seg
    real, allocatable, intent(inout) :: pts(:,:)
    integer, intent(out) :: npt, iend
    real :: p(3), pm(3), b(3), bh(3), bh0(3), s, slast, turn, bb
    integer :: n, nmax
    real, allocatable :: tmp(:,:)
    if (.not. allocated(pts)) allocate(pts(3, 1024))
    npt = 1 ; pts(:,1) = p0
    p = p0 ; s = 0.0 ; slast = 0.0 ; turn = 0.0
    nmax = int(maxlen/ds)
    iend = END_MAXLEN
    if (.not. fsample(fm, p, [VAR_BX, VAR_BY, VAR_BZ], b)) then
      iend = END_OUT ; return
    end if
    bb = sqrt(sum(b**2))
    if (bb <= 0.0) then
      iend = END_NULL ; return
    end if
    bh0 = dir*b/bb
    do n = 1, nmax
      bh = dir*b/sqrt(sum(b**2))
      pm = p + 0.5*ds*bh
      if (.not. fsample(fm, pm, [VAR_BX, VAR_BY, VAR_BZ], b)) then
        iend = END_OUT ; exit
      end if
      bb = sqrt(sum(b**2))
      if (bb <= 0.0) then
        iend = END_NULL ; exit
      end if
      bh = dir*b/bb
      p = p + ds*bh
      s = s + ds
      turn = turn + sqrt(sum((bh - bh0)**2)) ; bh0 = bh
      if (s - slast >= seg .or. turn > 0.15) then
        call push(p) ; slast = s ; turn = 0.0
      end if
      if (sum(p**2) < r_stop**2) then
        iend = END_PLANET ; exit
      end if
      ! back at the starting point (closed field-line loop)
      if (s > 20.0*ds .and. sum((p - p0)**2) < (1.5*ds)**2) then
        call push(p0) ; iend = END_LOOP ; exit
      end if
      if (.not. fsample(fm, p, [VAR_BX, VAR_BY, VAR_BZ], b)) then
        iend = END_OUT ; exit
      end if
    end do
    if (any(pts(:,npt) /= p)) call push(p)
  contains
    subroutine push(q)
      real, intent(in) :: q(3)
      if (npt == size(pts,2)) then
        allocate(tmp(3, 2*npt)) ; tmp(:,1:npt) = pts(:,1:npt)
        call move_alloc(tmp, pts)
      end if
      npt = npt + 1 ; pts(:,npt) = q
    end subroutine push
  end subroutine trace_line

end module mhdfield
