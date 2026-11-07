!==============================================================================
! slice_plot : PostScript plots of slices of combined 3-D data (m.data / .vort)
!
!   Replaces the legacy fig_h/xyju02nt171m.f + sub3ahpar2.f (F90, no recompilation needed)
!
!   * quantities, slices (xy/xz/yz and position), colour ranges, vectors and layout are
!     all given in a namelist (default slice.nml)
!   * only the planes and variables needed for the slices are read from the file
!       xy slice : 2 planes x needed variables (legacy: all 8 variables x full 3-D = 24 GB x several in memory)
!       xz slice : only 2 rows of each z plane
!       yz slice : each z plane read and 2 columns extracted (yz_read='strided' reads 2 words each)
!     memory use is a few 2-D planes (a few hundred MB or less on the standard grid)
!   * slice positions are distances from the planet center (linear interpolation between grid points),
!     identical to the legacy "z=0 plane = mean of k=1,2" and "x=0 plane = mean of i=ixc,ixc+1"
!   * z<0 is drawn by mirroring with the equatorial symmetry (parity) (as before, mirror_z)
!   * quantities: variable names (rho,vx,...,bz or names given in varname) and
!       temp (=p/rho), absb (|B|), absbz (|Bz|), absv (|V|), beta (2p/B^2),
!       pdyn (rho V^2), log10:<qty>, abs:<qty>, neg:<qty>
!   * one file = one page (ps_per_file=.true. gives a separate PS per file)
!   * dump_planes=.true. also writes the plotted 2-D planes as raw binary (for Python etc.)
!
! Usage:  ./slice_plot [namelist file]      (default: slice.nml)
!         for PNG/PDF conversion see scripts/ps2png.sh
!==============================================================================
program slice_plot
  use mhd_common
  use psplot
  use mhdunits
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  integer, parameter :: maxf = 1000, maxp = 16, maxv = 32

  !---- &slice : files, grid, page --------------------------------------------
  character(len=512) :: datafile(maxf) = ''      ! explicit list of files (takes precedence)
  character(len=512) :: diri = './'              ! or a name pattern
  character(len=64)  :: prefix = 'ju02nt', suffix = 'm.data'
  integer :: ntap = 171, nfiles = 1, nstep = 1, ndigit = 3
  character(len=128) :: title(maxf) = ''         ! title of each page (takes precedence)
  character(len=128) :: title_head = ''          ! automatic title: title_head // t // t_unit
  real               :: t0 = 0.0, dt = 0.0       !   t = t0 + (i-1)*dt
  character(len=16)  :: t_unit = ' h', t_fmt = '(f0.2)'
  character(len=128) :: header = ''           ! empty: chosen from the planet name
  character(len=512) :: psfile = 'slice.ps'
  logical :: ps_per_file = .false.
  ! grid
  integer :: nx2 = 3002, ny2 = 2002, nz2 = 1002, nvar = 8, nxp = 500
  real    :: h = 0.15                            ! grid spacing [planetary radii]
  real    :: i0 = -1.0, j0 = -1.0, k0 = -1.0     ! grid position of the planet center (negative: automatic)
  character(len=16) :: varname(maxv) = ''
  integer :: parity(maxv) = 1                    ! sign for z mirroring (+1 symmetric, -1 antisymmetric)
  character(len=8) :: varunit(maxv) = ''         ! kind of unit (n,v,p,b,omega,...; empty: guessed from the name)
  logical :: mirror_z = .true.
  character(len=8) :: yz_read = 'full'           ! 'full' or 'strided'
  ! page / layout [cm]
  integer :: idci = 2, ity = 2, isz = 1          ! cm, landscape, A4 (as before)
  real    :: x0 = 2.2, y0 = 1.8, panel_h = 6.3, dxl = 2.2, dyl = 2.6, h1 = 0.30
  character(len=8) :: label_coord = 'jso'        ! 'jso': X=-x,Y=-y,Z=z  / 'sim'
  real    :: tick = 50.0                         ! tick interval [planetary radii] (0: none)
  logical :: colorbar = .true.
  logical :: dump_planes = .false.
  namelist /slice/ datafile, diri, prefix, suffix, ntap, nfiles, nstep, ndigit, &
                   title, title_head, t0, dt, t_unit, t_fmt, header, psfile,   &
                   ps_per_file, nx2, ny2, nz2, nvar, nxp, h, i0, j0, k0,       &
                   varname, varunit, parity, mirror_z, yz_read, idci, ity, isz, x0, y0, &
                   panel_h, dxl, dyl, h1, label_coord, tick, colorbar,         &
                   dump_planes

  !---- &panels ---------------------------------------------------------------
  integer :: npanel = 1
  character(len=4)  :: plane(maxp) = 'xy'        ! 'xy', 'xz', 'yz'
  real              :: pos(maxp) = 0.0           ! slice position (sim coordinates, from the planet center)
  character(len=32) :: qty(maxp) = 'rho'
  real              :: vmin(maxp) = 0.0, vmax(maxp) = 0.0   ! equal: automatic
  integer           :: palette(maxp) = 3         ! legacy ico (1..8)
  integer           :: nsmooth(maxp) = 2         ! number of legacy diff2 passes (mm)
  integer           :: stride(maxp) = 1          ! image thinning (smaller PS)
  character(len=4)  :: vec(maxp) = 'none'        ! 'v', 'b', 'none'
  real              :: vecmax(maxp) = 0.0        ! vector normalization (0: automatic)
  integer           :: nvec_h(maxp) = 0, nvec_v(maxp) = 0   ! number of seed points (0: automatic)
  integer           :: vec_steps(maxp) = 41, vec_pen(maxp) = 21
  integer           :: row(maxp) = 1, col(maxp) = 1
  character(len=32) :: label(maxp) = ''
  integer           :: hflip(maxp) = -1, vflip(maxp) = -1   ! -1: default (yz mirrored horizontally)
  real              :: h_lo(maxp) = -1.0e30, h_hi(maxp) = 1.0e30  ! plotted range
  real              :: v_lo(maxp) = -1.0e30, v_hi(maxp) = 1.0e30  ! (sim coordinates)
  namelist /panels/ npanel, plane, pos, qty, vmin, vmax, palette, nsmooth,     &
                    stride, vec, vecmax, nvec_h, nvec_v, vec_steps, vec_pen,   &
                    row, col, label, hflip, vflip, h_lo, h_hi, v_lo, v_hi

  !---- extracted planes ------------------------------------------------------
  type :: pdata
    character(len=4) :: plane = ''
    real :: pos = 0.0
    integer :: nh = 0, nv = 0
    real, allocatable :: ch(:), cv(:)            ! sim coordinates
    integer :: slot(maxv) = 0                    ! variable number -> 3rd index of d
    real(real32), allocatable :: d(:,:,:)
  end type pdata
  type(pdata), allocatable :: pd(:)
  integer :: pmap(maxp), nuniq

  type(seqfile) :: fin
  integer :: u, ios, ifile, nf, p, ino, nrow, ncol, c, r
  integer(int64) :: nplane
  real :: wcol(maxp), xcol(maxp), yrow(maxp), ytop
  real(real64) :: ta, tr
  logical :: need(maxv, maxp)
  character(len=:), allocatable :: fni, nmlfile, psname
  character(len=512) :: arg
  character(len=256) :: msg
  character(len=128) :: ttl

  !---- defaults & namelist ---------------------------------------------------
  varname(1:8) = ['rho', 'vx ', 'vy ', 'vz ', 'p  ', 'bx ', 'by ', 'bz ']
  parity(1:8)  = [1, 1, 1, -1, 1, -1, -1, 1]     ! same as the boundary condition fbb
  nmlfile = 'slice.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=u, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(u, nml=slice, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &slice : '//trim(msg))
  rewind(u)
  read(u, nml=panels, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &panels : '//trim(msg))
  call units_read(u)
  close(u)
  if (len_trim(header) == 0) header = default_header()
  if (npanel < 1 .or. npanel > maxp) call fatal('npanel out of range')
  if (i0 < 0.0) i0 = real(nx2/2 - nxp) + 0.5
  if (j0 < 0.0) j0 = real(ny2/2) + 0.5
  if (k0 < 0.0) k0 = 1.5
  do c = 1, maxv
    varname(c) = lower(varname(c))
    if (len_trim(varname(c)) == 0) write(varname(c),'(a,i0)') 'var', c
    if (len_trim(varunit(c)) == 0) varunit(c) = kind_of_name(varname(c))
  end do
  nplane = int(nx2,int64)*ny2

  nf = count(len_trim(datafile) > 0)
  if (nf == 0) nf = nfiles
  if (nf < 1) call fatal('no input file: nfiles must be >= 1 (or give datafile(1), ...)')

  !---- panels: needed variables, unique planes, layout -----------------------
  need = .false.
  do p = 1, npanel
    plane(p) = lower(plane(p)) ; qty(p) = lower(qty(p)) ; vec(p) = lower(vec(p))
    if (all(plane(p) /= ['xy', 'xz', 'yz'])) call fatal('plane must be xy, xz or yz: '//plane(p))
    call mark_needs(qty(p), need(:,p))
    if (vec(p) == 'v') call mark_needs('vx', need(:,p))
    if (vec(p) == 'v') call mark_needs('vy', need(:,p))
    if (vec(p) == 'v') call mark_needs('vz', need(:,p))
    if (vec(p) == 'b') call mark_needs('bx', need(:,p))
    if (vec(p) == 'b') call mark_needs('by', need(:,p))
    if (vec(p) == 'b') call mark_needs('bz', need(:,p))
    if (hflip(p) < 0) hflip(p) = merge(1, 0, plane(p) == 'yz')
    if (vflip(p) < 0) vflip(p) = 0
    if (len_trim(label(p)) == 0) then
      label(p) = qty(p)
      if (len_trim(qty_unit(qty(p))) > 0) label(p) = trim(qty(p))//' ['//trim(qty_unit(qty(p)))//']'
    end if
  end do
  ! unique (plane,pos)
  allocate(pd(npanel))
  nuniq = 0
  do p = 1, npanel
    pmap(p) = 0
    do c = 1, nuniq
      if (pd(c)%plane == plane(p) .and. pd(c)%pos == pos(p)) pmap(p) = c
    end do
    if (pmap(p) == 0) then
      nuniq = nuniq + 1
      pd(nuniq)%plane = plane(p)
      pd(nuniq)%pos   = pos(p)
      pmap(p) = nuniq
    end if
  end do

  write(*,'(a)') '----- slice_plot -----'
  write(*,'(a,3i6,a,i0)') ' grid: ', nx2, ny2, nz2, '  nvar=', nvar
  write(*,'(a,3f9.2)') ' planet center (i0,j0,k0): ', i0, j0, k0
  call units_report()
  do p = 1, npanel
    write(*,'(a,i2,4a,f8.2,6a)') ' panel', p, ': ', plane(p), ' at (sim) ', &
         'xyz'(normal_axis(plane(p)):normal_axis(plane(p))) // ' =', pos(p), '  qty=', trim(qty(p)), &
         '  vec=', trim(vec(p)), '  unit: ', trim(qty_unit(qty(p)))
  end do

  !---- loop over files (1 page per file) -------------------------------------
  ino = 0
  do ifile = 1, nf
    if (len_trim(datafile(ifile)) > 0) then
      fni = trim(datafile(ifile))
    else
      fni = trim(make_fname(diri, prefix, ntap + (ifile-1)*nstep, ndigit, suffix))
    end if
    if (ps_per_file .or. ifile == 1) then
      if (ifile > 1) call ps_close()
      if (ps_per_file) then
        psname = trim(strip_ext(psfile))//'_'//trim(zpad(ifile,3))//'.ps'
      else
        psname = trim(psfile)
      end if
      call ps_open(psname)
      ino = 0
      if (psname(1:1) == '/') then
        write(*,'(a)') ' PostScript: '//psname
      else
        call get_environment_variable('PWD', arg)
        write(*,'(a)') ' PostScript: '//trim(arg)//'/'//psname
      end if
    end if
    write(*,'(a)') ' input: '//fni
    ta = wtime()

    !-- read planes
    call seq_open_read(fin, fni, 4_int64*nplane)
    call seq_check_size(fin, int(nvar,int64)*nz2)
    do c = 1, nuniq
      call extract(pd(c), c)
    end do
    call seq_close(fin)
    tr = wtime() - ta

    !-- layout
    nrow = maxval(row(1:npanel)) ; ncol = maxval(col(1:npanel))
    wcol = 0.0
    do p = 1, npanel
      wcol(col(p)) = max(wcol(col(p)), panel_h*aspect(p))
    end do
    xcol(1) = x0
    do c = 2, ncol
      xcol(c) = xcol(c-1) + wcol(c-1) + dxl
    end do
    do r = 1, nrow
      yrow(r) = y0 + real(nrow - r)*(panel_h + dyl)
    end do
    ytop = y0 + real(nrow)*panel_h + real(nrow-1)*dyl

    !-- page
    ino = ino + 1
    call ps_page(idci, ity, isz, ino)
    call ps_pen(2)
    if (len_trim(title(ifile)) > 0) then
      ttl = title(ifile)
    else if (len_trim(title_head) > 0) then
      write(ttl, t_fmt) t0 + real(ifile-1)*dt
      ttl = adjustl(ttl)
      if (ttl(1:1) == '.') ttl = '0'//ttl(1:len(ttl)-1)
      if (ttl(1:2) == '-.') ttl = '-0'//ttl(2:)
      ttl = trim(title_head)//' '//trim(ttl)//trim(t_unit)
    else
      ttl = basename(fni)
    end if
    call ps_text(x0, ytop + 0.95 + h1*2.2, h1*1.5, trim(header))
    call ps_text(x0, ytop + 0.95, h1*1.4, trim(ttl))
    do p = 1, npanel
      call draw_panel(p, xcol(col(p)), yrow(row(p)))
    end do
    call ps_endpage()
    write(*,'(a,f8.2,a,f8.2,a)') '  read', tr, ' s,  total', wtime()-ta, ' s'
  end do
  call ps_close()
  write(*,'(a)') ' done.'

contains

  !============================================================================
  ! helpers
  !============================================================================
  pure function lower(s) result(o)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i
    o = s
    do i = 1, len(s)
      if (iachar(s(i:i)) >= 65 .and. iachar(s(i:i)) <= 90) o(i:i) = achar(iachar(s(i:i))+32)
    end do
  end function lower

  function basename(s) result(o)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: o
    integer :: i
    i = index(s, '/', back=.true.)
    o = trim(s(i+1:))
  end function basename

  function strip_ext(s) result(o)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: o
    integer :: i
    i = index(s, '.', back=.true.)
    if (i > index(s, '/', back=.true.) .and. i > 1) then
      o = s(1:i-1)
    else
      o = trim(s)
    end if
  end function strip_ext

  integer function var_index(name)
    character(len=*), intent(in) :: name
    integer :: m
    var_index = 0
    do m = 1, nvar
      if (trim(varname(m)) == trim(name)) then
        var_index = m ; return
      end if
    end do
  end function var_index

  integer function req(name)          ! var index, error if absent
    character(len=*), intent(in) :: name
    req = var_index(name)
    if (req == 0) call fatal('variable "'//trim(name)//'" not found in varname list')
  end function req

  ! axes: 1=x, 2=y, 3=z
  pure integer function normal_axis(pl)
    character(len=*), intent(in) :: pl
    select case (pl)
    case ('xy') ; normal_axis = 3
    case ('xz') ; normal_axis = 2
    case default; normal_axis = 1
    end select
  end function normal_axis

  pure integer function h_axis(pl)
    character(len=*), intent(in) :: pl
    h_axis = merge(2, 1, pl == 'yz')
  end function h_axis

  pure integer function v_axis(pl)
    character(len=*), intent(in) :: pl
    v_axis = merge(2, 3, pl == 'xy')
  end function v_axis

  function axis_name(ax) result(s)
    integer, intent(in) :: ax
    character(len=1) :: s
    if (label_coord == 'jso') then
      s = 'XYZ'(ax:ax)
    else
      s = 'xyz'(ax:ax)
    end if
  end function axis_name

  ! sim -> displayed coordinate sign (jso: X=-x, Y=-y, Z=z)
  pure real function dispsign(ax)
    integer, intent(in) :: ax
    dispsign = 1.0
    if (label_coord == 'jso' .and. ax <= 2) dispsign = -1.0
  end function dispsign

  !----------------------------------------------------------------------------
  ! which variables does quantity q need?
  !----------------------------------------------------------------------------
  recursive subroutine mark_needs(q, nd)
    character(len=*), intent(in) :: q
    logical, intent(inout) :: nd(maxv)
    character(len=32) :: s
    s = q
    if (index(s, ':') > 0) then
      call mark_needs(s(index(s,':')+1:), nd) ; return
    end if
    select case (trim(s))
    case ('temp')
      nd(req('rho')) = .true. ; nd(req('p')) = .true.
    case ('absb', 'beta')
      nd(req('bx')) = .true. ; nd(req('by')) = .true. ; nd(req('bz')) = .true.
      if (trim(s) == 'beta') nd(req('p')) = .true.
    case ('absbz')
      nd(req('bz')) = .true.
    case ('absv', 'pdyn')
      nd(req('vx')) = .true. ; nd(req('vy')) = .true. ; nd(req('vz')) = .true.
      if (trim(s) == 'pdyn') nd(req('rho')) = .true.
    case default
      nd(req(s)) = .true.
    end select
  end subroutine mark_needs

  !----------------------------------------------------------------------------
  ! kind of unit of quantity q (decided after removing operators such as log10:)
  !----------------------------------------------------------------------------
  recursive function qty_kind(q) result(k)
    character(len=*), intent(in) :: q
    character(len=8) :: k
    integer :: ic
    ic = index(q, ':')
    if (ic > 0) then
      k = qty_kind(q(ic+1:)) ; return
    end if
    select case (trim(q))
    case ('temp')           ; k = 'temp'
    case ('absb', 'absbz')  ; k = 'b'
    case ('absv')           ; k = 'v'
    case ('beta')           ; k = 'beta'
    case ('pdyn')           ; k = 'pdyn'
    case default            ; k = varunit(req(q))
    end select
  end function qty_kind

  function qty_unit(q) result(s)
    character(len=*), intent(in) :: q
    character(len=24) :: s
    s = ulabel(qty_kind(q))
    if (index(q, 'log10:') > 0 .and. len_trim(s) > 0) s = 'log10 '//trim(s)
  end function qty_unit

  !----------------------------------------------------------------------------
  ! evaluate quantity q on plane data P -> out(nh,nv)
  !   derived quantities are computed in double precision and converted to display units (mhdunits)
  !----------------------------------------------------------------------------
  recursive subroutine eval_qty(q, P, out)
    character(len=*), intent(in) :: q
    type(pdata), intent(in) :: P
    real(real32), intent(out) :: out(:,:)
    character(len=32) :: s, op
    real(real64), allocatable :: d(:,:)
    real(real64) :: f
    integer :: ic
    s = q
    ic = index(s, ':')
    if (ic > 0) then
      op = s(1:ic-1)
      call eval_qty(s(ic+1:), P, out)
      select case (trim(op))
      case ('log10') ; out = log10(max(out, tiny(1.0_real32)))
      case ('abs')   ; out = abs(out)
      case ('neg')   ; out = -out
      case default   ; call fatal('unknown operator '//trim(op))
      end select
      return
    end if
    f = ufac(qty_kind(s))
    select case (trim(s))
    case ('temp')
      d = dv(P,'p') / max(dv(P,'rho'), 1.0e-30_real64)
    case ('absb')
      d = sqrt(dv(P,'bx')**2 + dv(P,'by')**2 + dv(P,'bz')**2)
    case ('beta')
      d = dv(P,'p') / max(dv(P,'bx')**2 + dv(P,'by')**2 + dv(P,'bz')**2, 1.0e-200_real64)
    case ('absbz')
      d = abs(dv(P,'bz'))
    case ('absv')
      d = sqrt(dv(P,'vx')**2 + dv(P,'vy')**2 + dv(P,'vz')**2)
    case ('pdyn')
      d = dv(P,'rho')*(dv(P,'vx')**2 + dv(P,'vy')**2 + dv(P,'vz')**2)
    case default
      d = dv(P,s)
    end select
    out = real(d*f, real32)
  end subroutine eval_qty

  function dv(P, nm) result(a)
    type(pdata), intent(in) :: P
    character(len=*), intent(in) :: nm
    real(real64) :: a(P%nh, P%nv)
    a = real(P%d(:,:,P%slot(req(nm))), real64)
  end function dv

  !============================================================================
  ! extraction of a plane from the 3-D file
  !============================================================================
  subroutine extract(P, ip)
    type(pdata), intent(inout) :: P
    integer, intent(in) :: ip
    logical :: nd(maxv)
    integer :: m, ns, p2, i1, k, j
    real :: w, sgn, cpos
    real(real32), allocatable :: a(:), b(:), raw(:,:,:), buf(:,:)
    logical :: neg

    nd = .false.
    do p2 = 1, npanel
      if (pmap(p2) == ip) nd = nd .or. need(:,p2)
    end do
    ns = 0 ; P%slot = 0
    do m = 1, nvar
      if (nd(m)) then
        ns = ns + 1 ; P%slot(m) = ns
      end if
    end do
    if (allocated(P%d))  deallocate(P%d)
    if (allocated(P%ch)) deallocate(P%ch)
    if (allocated(P%cv)) deallocate(P%cv)

    select case (P%plane)
    !------------------------------------------------------------------ xy
    case ('xy')
      cpos = P%pos ; neg = .false.
      if (mirror_z .and. cpos < 0.0) then
        cpos = -cpos ; neg = .true.
      end if
      call weights(k0 + cpos/h, nz2, i1, w)
      P%nh = nx2 ; P%nv = ny2
      allocate(P%d(nx2, ny2, ns), a(nplane), b(nplane))
      do m = 1, nvar
        if (P%slot(m) == 0) cycle
        call seq_read_rec(fin, rec(m,i1), a, nplane)
        if (w > 0.0) then
          call seq_read_rec(fin, rec(m,i1+1), b, nplane)
          a = (1.0-w)*a + w*b
        end if
        sgn = 1.0 ; if (neg) sgn = real(parity(m))
        P%d(:,:,P%slot(m)) = sgn*reshape(a, [nx2, ny2])
      end do
      allocate(P%ch(nx2), P%cv(ny2))
      P%ch = [((real(j)-i0)*h, j=1,nx2)]
      P%cv = [((real(j)-j0)*h, j=1,ny2)]
    !------------------------------------------------------------------ xz
    case ('xz')
      call weights(j0 + P%pos/h, ny2, i1, w)
      allocate(raw(nx2, nz2, ns), buf(nx2, 2))
      do m = 1, nvar
        if (P%slot(m) == 0) cycle
        do k = 1, nz2
          call seq_read_part(fin, rec(m,k), int(i1-1,int64)*nx2, buf, 2_int64*nx2)
          raw(:,k,P%slot(m)) = (1.0-w)*buf(:,1) + w*buf(:,2)
        end do
      end do
      call finish_vertical(P, raw, nx2, [((real(j)-i0)*h, j=1,nx2)])
    !------------------------------------------------------------------ yz
    case ('yz')
      call weights(i0 + P%pos/h, nx2, i1, w)
      allocate(raw(ny2, nz2, ns))
      if (yz_read == 'strided') then
        allocate(a(2))
        do m = 1, nvar
          if (P%slot(m) == 0) cycle
          do k = 1, nz2
            do j = 1, ny2
              call seq_read_part(fin, rec(m,k), int(j-1,int64)*nx2 + i1 - 1, a, 2_int64)
              raw(j,k,P%slot(m)) = (1.0-w)*a(1) + w*a(2)
            end do
          end do
        end do
      else
        allocate(buf(nx2, ny2))
        do m = 1, nvar
          if (P%slot(m) == 0) cycle
          do k = 1, nz2
            call seq_read_rec(fin, rec(m,k), buf, nplane)
            raw(:,k,P%slot(m)) = (1.0-w)*buf(i1,:) + w*buf(i1+1,:)
          end do
        end do
      end if
      call finish_vertical(P, raw, ny2, [((real(j)-j0)*h, j=1,ny2)])
    end select

  end subroutine extract

  ! vertical axis = z : optional mirroring to z<0 using parity
  subroutine finish_vertical(P, raw, nh, chv)
    type(pdata), intent(inout) :: P
    real(real32), intent(in) :: raw(:,:,:)
    integer, intent(in) :: nh
    real, intent(in) :: chv(nh)
    integer :: m, r, ns, nv
    ns = size(raw, 3)
    P%nh = nh
    allocate(P%ch(nh)) ; P%ch = chv
    if (mirror_z) then
      ! rows 1..nz2-1 : k = nz2..2 (mirrored), rows nz2..2nz2-2 : k = 2..nz2
      ! (k=1 is the ghost plane outside z=0 and is not used; same layout as the legacy code)
      nv = 2*(nz2-1)
      allocate(P%d(nh, nv, ns), P%cv(nv))
      do m = 1, nvar
        if (P%slot(m) == 0) cycle
        do r = 1, nz2-1
          P%d(:, r, P%slot(m)) = real(parity(m))*raw(:, nz2+1-r, P%slot(m))
        end do
        do r = nz2, nv
          P%d(:, r, P%slot(m)) = raw(:, r-nz2+2, P%slot(m))
        end do
      end do
      P%cv = [((real(r - nz2 + 2) - k0)*h, r=1,nv)]
    else
      nv = nz2
      allocate(P%d(nh, nv, ns), P%cv(nv))
      P%d = raw
      P%cv = [((real(r) - k0)*h, r=1,nv)]
    end if
    P%nv = nv
  end subroutine finish_vertical

  ! fractional index -> lower index i1 and weight w (clamped)
  subroutine weights(fidx, n, i1, w)
    real, intent(in) :: fidx
    integer, intent(in) :: n
    integer, intent(out) :: i1
    real, intent(out) :: w
    i1 = floor(fidx)
    w  = fidx - real(i1)
    if (i1 < 1) then
      i1 = 1 ; w = 0.0
    else if (i1 >= n) then
      i1 = n-1 ; w = 1.0
    end if
    if (w < 1.0e-6) w = 0.0
    if (w > 1.0-1.0e-6 .and. i1 < n-1) then
      i1 = i1 + 1 ; w = 0.0
    end if
  end subroutine weights

  pure integer(int64) function rec(m, k)
    integer, intent(in) :: m, k
    rec = int(m-1,int64)*nz2 + k
  end function rec

  !============================================================================
  ! panel geometry (after zoom)
  !============================================================================
  subroutine panel_range(p, ih0, ih1, iv0, iv1)
    integer, intent(in) :: p
    integer, intent(out) :: ih0, ih1, iv0, iv1
    integer :: i
    associate(P0 => pd(pmap(p)))
      ih0 = 1 ; ih1 = P0%nh ; iv0 = 1 ; iv1 = P0%nv
      do i = 1, P0%nh
        if (P0%ch(i) >= h_lo(p)) then ; ih0 = i ; exit ; end if
      end do
      do i = P0%nh, 1, -1
        if (P0%ch(i) <= h_hi(p)) then ; ih1 = i ; exit ; end if
      end do
      do i = 1, P0%nv
        if (P0%cv(i) >= v_lo(p)) then ; iv0 = i ; exit ; end if
      end do
      do i = P0%nv, 1, -1
        if (P0%cv(i) <= v_hi(p)) then ; iv1 = i ; exit ; end if
      end do
    end associate
    if (ih1 - ih0 < 2 .or. iv1 - iv0 < 2) call fatal('panel range too small')
  end subroutine panel_range

  real function aspect(p)
    integer, intent(in) :: p
    integer :: a, b, c, d
    call panel_range(p, a, b, c, d)
    aspect = real(b-a)/real(d-c)
  end function aspect

  !============================================================================
  ! draw one panel with lower-left corner (xb,yb)
  !============================================================================
  subroutine draw_panel(p, xb, yb)
    integer, intent(in) :: p
    real, intent(in) :: xb, yb
    real(real32), allocatable :: full(:,:), q(:,:), vh(:,:), vv(:,:), img(:,:)
    integer :: ih0, ih1, iv0, iv1, nh, nv, ah, av, nsx, nsy, uu
    real :: xl, yl, lo, hi, umax, xpl, ypl, chs(2), cvs(2), vfac
    character(len=64) :: s1, s2
    character(len=:), allocatable :: dname

    associate(PL => pd(pmap(p)))
    call panel_range(p, ih0, ih1, iv0, iv1)
    nh = ih1 - ih0 + 1 ; nv = iv1 - iv0 + 1
    ah = h_axis(plane(p)) ; av = v_axis(plane(p))
    yl = panel_h
    xl = panel_h*real(nh-1)/real(nv-1)

    !-- scalar
    allocate(full(PL%nh, PL%nv))
    call eval_qty(qty(p), PL, full)
    q = full(ih0:ih1, iv0:iv1)
    deallocate(full)
    if (hflip(p) == 1) q = q(nh:1:-1, :)
    if (vflip(p) == 1) q = q(:, nv:1:-1)

    if (dump_planes) then
      dname = trim(strip_ext(psname))//'_f'//trim(zpad(ifile,3))//'_p'//trim(zpad(p,2))
      open(newunit=uu, file=dname//'.bin', access='stream', form='unformatted', status='replace')
      write(uu) q
      close(uu)
      open(newunit=uu, file=dname//'.txt', status='replace')
      write(uu,'(a)')    '# slice_plot plane dump (float32, native endian, Fortran order q(nh,nv))'
      write(uu,'(a,a)')  'file   = ', fni
      write(uu,'(a,a)')  'plane  = ', plane(p)
      write(uu,'(a,g0)') 'pos    = ', pos(p)
      write(uu,'(a,a)')  'qty    = ', trim(qty(p))
      write(uu,'(a,a)')  'unit   = ', trim(qty_unit(qty(p)))
      write(uu,'(a,i0,1x,i0)') 'nh nv  = ', nh, nv
      write(uu,'(a,2g16.7)') 'h axis (sim, '//trim(ulen())//') first,last = ', &
           merge(PL%ch(ih1), PL%ch(ih0), hflip(p)==1), merge(PL%ch(ih0), PL%ch(ih1), hflip(p)==1)
      write(uu,'(a,2g16.7)') 'v axis (sim, '//trim(ulen())//') first,last = ', &
           merge(PL%cv(iv1), PL%cv(iv0), vflip(p)==1), merge(PL%cv(iv0), PL%cv(iv1), vflip(p)==1)
      write(uu,'(a)')    'smoothing = none'
      close(uu)
    end if

    lo = vmin(p) ; hi = vmax(p)
    if (lo == hi) then
      lo = minval(q, mask=.not. ieee_is_nan(q)) ; hi = maxval(q, mask=.not. ieee_is_nan(q))
      if (lo == hi) hi = lo + 1.0
    end if

    call smooth2d(nsmooth(p), nh, nv, q)
    if (stride(p) > 1) then
      img = q(1:nh:stride(p), 1:nv:stride(p))
    else
      call move_alloc(q, img)
    end if
    call ps_image(size(img,1), size(img,2), xb, yb, xl, yl, img, lo, hi, palette(p), 64)
    deallocate(img)

    !-- vectors
    if (vec(p) == 'v' .or. vec(p) == 'b') then
      s1 = merge('vx', 'bx', vec(p) == 'v')
      s2 = merge('vy', 'by', vec(p) == 'v')
      if (ah == 2) s1 = merge('vy', 'by', vec(p) == 'v')
      if (av == 3) s2 = merge('vz', 'bz', vec(p) == 'v')
      vfac = real(ufac(merge('v', 'b', vec(p) == 'v')))
      vh = vfac*PL%d(ih0:ih1, iv0:iv1, PL%slot(req(s1)))
      vv = vfac*PL%d(ih0:ih1, iv0:iv1, PL%slot(req(s2)))
      if (hflip(p) == 1) then
        vh = -vh(nh:1:-1, :) ; vv = vv(nh:1:-1, :)
      end if
      if (vflip(p) == 1) then
        vh = vh(:, nv:1:-1) ; vv = -vv(:, nv:1:-1)
      end if
      umax = vecmax(p)
      if (umax <= 0.0) umax = maxval(sqrt(vh**2 + vv**2))
      nsy = nvec_v(p) ; nsx = nvec_h(p)
      if (nsy <= 0) nsy = 20
      if (nsx <= 0) nsx = max(1, nint(real(nsy)*xl/yl))
      call ps_vectors(nh, nv, vh, vv, nsx, nsy, xb, yb, xl, yl, umax, 0.01, &
                      vec_steps(p), 1, vec_pen(p))
      deallocate(vh, vv)
    end if

    !-- frame, planet lines
    call ps_pen(2) ; call ps_rgb(0.0, 0.0, 0.0)
    call ps_box(xb, yb, xl, yl) ; call ps_stroke()
    chs = [PL%ch(ih0), PL%ch(ih1)] ; cvs = [PL%cv(iv0), PL%cv(iv1)]
    if (hflip(p) == 1) chs = chs(2:1:-1)
    if (vflip(p) == 1) cvs = cvs(2:1:-1)
    call ps_pen(1) ; call ps_rgb(0.0, 0.0, 1.0)
    xpl = xb + xl*(0.0-chs(1))/(chs(2)-chs(1))
    ypl = yb + yl*(0.0-cvs(1))/(cvs(2)-cvs(1))
    if (xpl > xb .and. xpl < xb+xl) call ps_line(xpl, yb, xpl, yb+yl)
    if (ypl > yb .and. ypl < yb+yl) call ps_line(xb, ypl, xb+xl, ypl)
    call ps_stroke() ; call ps_rgb(0.0, 0.0, 0.0)

    !-- ticks & labels (displayed coordinates)
    call ps_pen(1)
    if (tick > 0.0) then
      call axis_ticks(xb, yb, xl, yl, chs, ah, .true.)
      call axis_ticks(xb, yb, xl, yl, cvs, av, .false.)
    end if
    call ps_text(xb+xl*0.5, yb-0.95, h1, axis_name(ah)//' ['//trim(ulen())//']', align='c')
    call ps_text(xb-0.75, yb+yl*0.5, h1, axis_name(av)//' ['//trim(ulen())//']', ang=90.0, align='c')
    ! panel title
    call ps_text(xb, yb+yl+0.15, h1*1.2, trim(label(p)))
    write(s1,'(f0.1)') dispsign(normal_axis(plane(p)))*pos(p) + 0.0
    if (abs(pos(p)) < 0.05) s1 = '0.0'
    if (s1(1:1) == '.') s1 = '0'//trim(s1)
    if (s1(1:2) == '-.') s1 = '-0'//trim(s1(2:))
    call ps_text(xb+xl, yb+yl+0.15, h1, axis_name(normal_axis(plane(p)))//' = '// &
                 trim(s1)//' '//trim(ulen()), align='r')
    if (vec(p) == 'v' .or. vec(p) == 'b') then
      write(s2,'(es9.2)') umax
      call ps_text(xb+0.25+0.25*len_trim(label(p)), yb+yl+0.15, h1*0.75, &
                   '(vec '//trim(vec(p))//': '//trim(adjustl(s2))//' '// &
                   trim(ulabel(merge('v', 'b', vec(p) == 'v')))//')')
    end if
    !-- colour bar
    if (colorbar) then
      call ps_colorbar(xb+xl*0.25, yb-1.6, xl*0.5, 0.22, palette(p))
      write(s1,'(es10.3)') lo ; write(s2,'(es10.3)') hi
      call ps_text(xb+xl*0.25-0.1, yb-1.58, h1*0.8, trim(adjustl(s1)), align='r')
      call ps_text(xb+xl*0.75+0.1, yb-1.58, h1*0.8, trim(adjustl(s2)))
      if (len_trim(qty_unit(qty(p))) > 0) &
        call ps_text(xb+xl*0.5, yb-1.95, h1*0.75, '['//trim(qty_unit(qty(p)))//']', align='c')
    end if
    end associate
  end subroutine draw_panel

  !----------------------------------------------------------------------------
  ! ticks every `tick` RJ along one panel edge; c(1:2) = sim coord at both ends
  !----------------------------------------------------------------------------
  subroutine axis_ticks(xb, yb, xl, yl, c, ax, horizontal)
    real, intent(in) :: xb, yb, xl, yl, c(2)
    integer, intent(in) :: ax
    logical, intent(in) :: horizontal
    real :: cmin, cmax, t, f, xp, yp, tl
    integer :: n0, n1, n
    character(len=16) :: s
    tl = 0.12
    cmin = min(c(1), c(2)) ; cmax = max(c(1), c(2))
    n0 = ceiling(cmin/tick) ; n1 = floor(cmax/tick)
    do n = n0, n1
      t = real(n)*tick
      f = (t - c(1))/(c(2) - c(1))
      write(s,'(i0)') nint(dispsign(ax)*t)
      if (horizontal) then
        xp = xb + xl*f
        call ps_line(xp, yb, xp, yb+tl)
        call ps_line(xp, yb+yl, xp, yb+yl-tl)
        call ps_text(xp, yb-0.45, h1*0.85, trim(s), align='c')
      else
        yp = yb + yl*f
        call ps_line(xb, yp, xb+tl, yp)
        call ps_line(xb+xl, yp, xb+xl-tl, yp)
        call ps_text(xb-0.1, yp-0.1, h1*0.85, trim(s), align='r')
      end if
    end do
    call ps_stroke()
  end subroutine axis_ticks

end program slice_plot
