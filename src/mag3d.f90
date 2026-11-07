!==============================================================================
! mag3d : 3-D visualization of magnetospheric MHD data (field lines, slices, isosurfaces, vectors)
!
!   Successor of the VRML tools of the late Prof. Tatsuki Ogino.
!   One scene per time step, written as
!     * binary glTF 2.0 (.glb)    : default. Opens in viewer/viewer.html (browser), Blender,
!                                    ParaView and the Windows 3D viewer
!     * X3D (.x3d + .png)         : the official successor of VRML. Opens in viewer/x3d_view.html
!                                    (X_ITE) and X3D software
!   Contents:
!     - field lines: from a polar foot-point grid (latitude x MLT) / an equatorial grid / a 3-D grid / any points.
!               coloured by topology (closed green, open blue, solar wind magenta, detached grey) or by a quantity.
!               lines near the open-closed boundary can be selected (legacy zsub33d)
!     - slices: colour images on xy/xz/yz planes (optionally translucent); "stacked slices" of equally
!             spaced parallel planes in up to 3 directions (as in the legacy code), and "selectable slices" shown one at a time
!     - isosurfaces: marching tetrahedra with shared vertices and gradient normals
!     - vector arrows (velocity, B), frame, axes, labels, colour bars, planet with an image (Saturn's rings)
!   Large data: the global field is read with stride thinning and the region near the planet at full
!   resolution, only the variables needed. Field-line tracing and isosurfaces are OpenMP parallel.
!
! Usage:  ./mag3d [namelist file]      (default: mag3d.nml)
!==============================================================================
program mag3d_main
  use mhd_common
  use mhdfield
  use mhdunits
  use scene3d
  use psplot, only: palette_rgb
  use isosurf, only: iso_extract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  integer, parameter :: maxf = 1000, mp = 16, mi = 8, mk = 3
  type :: lines_store
    integer :: n = 0, topo = 2          ! 0 closed, 1 open, 2 detached/undetermined, 3 solar wind
    real, allocatable :: p(:,:)
  end type lines_store
  real, parameter :: pi = 3.14159265
  !---- files ------------------------------------------------------------------
  character(len=512) :: datafile(maxf) = '', diri = './'
  character(len=64)  :: prefix = 'ju02nt', suffix = 'm.data'
  integer :: ntap = 171, nfiles = 1, nstep = 1, ndigit = 3
  character(len=128) :: title(maxf) = '', title_head = ''
  real :: t0 = 0.0, dt = 0.0
  character(len=16) :: t_unit = ' h', t_fmt = '(f0.2)'
  character(len=512) :: outbase = 'mag3d'
  logical :: out_glb = .true., out_x3d = .false.
  !---- grid / data -------------------------------------------------------------
  integer :: nx2 = 3002, ny2 = 2002, nz2 = 1002, nvar = 8, nxp = 500
  real    :: h = 0.15, i0 = -1.0, j0 = -1.0, k0 = -1.0
  logical :: half_domain = .true.
  integer :: stride = 4                       ! thinning of the global field
  real    :: box_lo(3) = -1.0e30, box_hi(3) = 1.0e30   ! region to read [sim coordinates]
  real    :: r_fine = 15.0                    ! full resolution inside this radius (0: off)
  real    :: b0 = -1.0, r_in = 7.0
  !---- display ------------------------------------------------------------------
  character(len=8) :: coord = 'jso'           ! 'jso': X=-x, Y=-y, Z=z / 'sim'
  real    :: bg(3) = [0.0, 0.0, 0.0]
  logical :: draw_frame = .true., draw_planet = .true., draw_inner = .false.
  real    :: label_size = 4.0
  ! planet image: 'auto' (../textures/<planet>.jpg|png relative to the executable), a file name, 'none' (plain colour)
  character(len=512) :: planet_texture = 'auto'
  real    :: planet_lon0 = 0.0                ! direction of image longitude 0 [deg] (0: toward the Sun)
  logical :: draw_rings = .true.              ! Saturn's rings (when planet='saturn')
  real    :: ring_r1 = 1.239, ring_r2 = 2.27  ! ring radii [planetary radii] (C ring inner to A ring outer edge)
  logical :: colorbar_3d = .false.            ! also draw colour bars inside the 3-D scene (for X3D, Blender)
  logical :: plane_labels = .true.               ! label "quantity [unit] range" above each slice
  !---- field lines ----------------------------------------------------------------
  logical :: fl_polar = .true.
  character(len=8) :: fl_hemi = 'both'        ! 'north' / 'south' / 'both'
  real    :: fl_lat_min = 70.0, fl_lat_max = 88.0
  integer :: fl_nlat = 7, fl_nmlt = 24
  real    :: fl_r_start = 10.0                ! radius to which foot points are mapped along the dipole
  logical :: fl_dipole_part = .true.          ! draw the dipole part from the surface to r_start
  character(len=8) :: fl_select = 'all'       ! 'all' / 'ocb' (only near the open-closed boundary)
  logical :: fl_equator = .false.
  real    :: eq_xmin = -30.0, eq_xmax = 60.0, eq_ymin = -30.0, eq_ymax = 30.0
  integer :: eq_nx = 10, eq_ny = 7
  integer :: fl_npt = 0
  ! from a 3-D grid of points (field lines not connected to the planet: solar wind, loops, etc.)
  logical :: fl_box = .false.
  real    :: fl_box_lo(3) = [-60.0, -60.0, -30.0], fl_box_hi(3) = [150.0, 60.0, 30.0]
  integer :: fl_box_n(3) = [8, 6, 3]
  ! which lines started from the equator / points / 3-D grid are kept:
  !   'all' / 'unconnected' (neither end on the planet) / 'connected' / 'closed' / 'open' /
  !   'solarwind' (both ends leave the domain) / 'detached' (loops, undetermined)
  character(len=16) :: fl_keep = 'all'
  real    :: fl_px(100) = 0.0, fl_py(100) = 0.0, fl_pz(100) = 0.0
  character(len=32) :: fl_color = 'topology'  ! 'topology' or a quantity name
  real    :: fl_vmin = 0.0, fl_vmax = 0.0
  integer :: fl_palette = 3
  real    :: fl_ds = 0.0, fl_maxlen = 2000.0, fl_seg = 1.0
  real    :: col_closed(3) = [0.0, 1.0, 0.0], col_open(3) = [0.0, 0.4, 1.0], &
             col_sw(3) = [1.0, 0.0, 1.0], col_undet(3) = [0.6, 0.6, 0.6]
  !---- planes -----------------------------------------------------------------------
  character(len=4)  :: pl_plane(mp) = ''      ! 'xy','xz','yz'
  real              :: pl_pos(mp) = 0.0
  character(len=32) :: pl_qty(mp) = ''        ! required for every index given in pl_plane
  real              :: pl_vmin(mp) = 0.0, pl_vmax(mp) = 0.0
  integer           :: pl_palette(mp) = 3
  real              :: pl_alpha(mp) = 0.85
  character(len=8)  :: pl_alpha_mode(mp) = 'const'   ! 'const' / 'value' (legacy pixect)
  real              :: pl_tra(mp) = 2.0              ! 'value': alpha = min(u*tra, 0.5*tra, 1)
  real              :: pl_h_lo(mp) = -1.0e30, pl_h_hi(mp) = 1.0e30
  real              :: pl_v_lo(mp) = -1.0e30, pl_v_hi(mp) = 1.0e30
  ! stacked slices (n parallel planes equally spaced in lo..hi)
  !   up to 3 sets (e.g. stk_plane = 'yz', 'xz', 'xy'). For one set, stk_plane = 'yz' as before
  character(len=4)  :: stk_plane(mk) = ''
  integer           :: stk_n(mk) = 11
  real              :: stk_lo(mk) = 0.0, stk_hi(mk) = 0.0
  character(len=32) :: stk_qty(mk) = ''
  real              :: stk_vmin(mk) = 0.0, stk_vmax(mk) = 0.0, stk_tra(mk) = 2.0
  real              :: stk_alpha(mk) = 0.4     ! maximum opacity (small, since overlapping planes add up)
  integer           :: stk_palette(mk) = 3
  ! selectable slices: sel_n parallel planes; the viewer slider shows one of them
  integer, parameter :: ms = 3
  character(len=4)  :: sel_plane(ms) = ''
  integer           :: sel_n(ms) = 21
  real              :: sel_lo(ms) = 0.0, sel_hi(ms) = 0.0
  character(len=32) :: sel_qty(ms) = ''
  real              :: sel_vmin(ms) = 0.0, sel_vmax(ms) = 0.0, sel_alpha(ms) = 0.9
  integer           :: sel_palette(ms) = 3
  !---- isosurfaces ----------------------------------------------------------------------
  character(len=32) :: iso_qty(mi) = ''
  real              :: iso_level(mi) = 0.0
  real              :: iso_color(3,mi) = 0.8, iso_alpha(mi) = 0.5
  integer           :: iso_stride = 1         ! extra thinning of the global field
  integer           :: iso_max_tri = 2000000  ! maximum triangles per isosurface (thinned automatically above this)
  integer           :: iso_smooth = 1         ! smoothing passes before extraction (1-2-1 filter, fewer small fragments)
  integer           :: ist = 1                ! thinning actually used (internal)
  !---- arrows -------------------------------------------------------------------------
  character(len=4)  :: arr_plane = ''
  real              :: arr_pos = 0.0
  character(len=4)  :: arr_vec = 'v'
  integer           :: arr_n = 30
  real              :: arr_scale = 0.0, arr_vmax = 0.0
  real              :: arr_h_lo = -1.0e30, arr_h_hi = 1.0e30, arr_v_lo = -1.0e30, arr_v_hi = 1.0e30

  namelist /mag3d/ datafile, diri, prefix, suffix, ntap, nfiles, nstep, ndigit, title,     &
       title_head, t0, dt, t_unit, t_fmt, outbase, out_glb, out_x3d, nx2, ny2, nz2, nvar,  &
       nxp, h, i0, j0, k0, half_domain, stride, box_lo, box_hi, r_fine, b0, r_in, coord,   &
       bg, draw_frame, draw_planet, draw_inner, label_size, plane_labels, planet_texture,     &
       planet_lon0, draw_rings, ring_r1, ring_r2, colorbar_3d, fl_box, fl_box_lo, fl_box_hi,   &
       fl_box_n, fl_keep, sel_plane, sel_n, sel_lo, sel_hi, sel_qty, sel_vmin, sel_vmax,       &
       sel_alpha, sel_palette, fl_polar, fl_hemi, fl_lat_min, &
       fl_lat_max, fl_nlat, fl_nmlt, fl_r_start, fl_dipole_part, fl_select, fl_equator,     &
       eq_xmin, eq_xmax, eq_ymin, eq_ymax, eq_nx, eq_ny, fl_npt, fl_px, fl_py, fl_pz,       &
       fl_color, fl_vmin, fl_vmax, fl_palette, fl_ds, fl_maxlen, fl_seg, col_closed,        &
       col_open, col_sw, col_undet, pl_plane, pl_pos, pl_qty, pl_vmin, pl_vmax, pl_palette, &
       pl_alpha, pl_alpha_mode, pl_tra, pl_h_lo, pl_h_hi, pl_v_lo, pl_v_hi, stk_plane,      &
       stk_n, stk_lo, stk_hi, stk_qty, stk_vmin, stk_vmax, stk_tra, stk_palette, stk_alpha, iso_qty,   &
       iso_level, iso_color, iso_alpha, iso_stride, iso_max_tri, iso_smooth, arr_plane, arr_pos, arr_vec, arr_n,     &
       arr_scale, arr_vmax, arr_h_lo, arr_h_hi, arr_v_lo, arr_v_hi

  type(field_t) :: fm
  type(scene_t) :: sc
  type(seqfile) :: fin
  logical :: need(32)
  integer, allocatable :: vlist(:)
  integer :: u, ios, ifile, nf, k
  real :: dsgn
  real(real64) :: ta, tb
  character(len=:), allocatable :: fni, nmlfile, ob
  character(len=512) :: arg
  character(len=256) :: msg
  character(len=128) :: ttl
  type(qty_t) :: qfl, qtmp

  nmlfile = 'mag3d.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=u, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(u, nml=mag3d, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &mag3d : '//trim(msg))
  call units_read(u)
  close(u)
  call grid_init(fm%g, nx2, ny2, nz2, nvar, h, i0, j0, k0, nxp, half_domain)
  dsgn = merge(-1.0, 1.0, coord == 'jso')
  if (fl_ds <= 0.0) fl_ds = 0.5*h*real(stride)

  ! ---- variables needed (collected here directly, without internal procedures)
  need = .false.
  need(VAR_BX) = .true. ; need(VAR_BY) = .true. ; need(VAR_BZ) = .true.
  if (fl_color /= 'topology') then
    call qty_parse(fl_color, qfl)
    need(1:8) = need(1:8) .or. qty_needs(qfl)
  end if
  do k = 1, mp
    if (len_trim(pl_plane(k)) > 0 .and. len_trim(pl_qty(k)) == 0) &
      call fatal('pl_qty('//trim(zpad(k,1))//') is not given for pl_plane('//trim(zpad(k,1))// &
                 ') (check the index numbers in the namelist)')
    if (len_trim(pl_plane(k)) == 0 .and. len_trim(pl_qty(k)) > 0) &
      write(*,'(a)') ' WARNING: pl_qty('//trim(zpad(k,1))//') is given but pl_plane('//trim(zpad(k,1))// &
                     ') is not: this plane is not drawn'
    if (len_trim(pl_plane(k)) > 0) then
      call qty_parse(pl_qty(k), qtmp) ; need(1:8) = need(1:8) .or. qty_needs(qtmp)
    end if
  end do
  do k = 1, mk
    if (len_trim(stk_plane(k)) > 0) then
      if (len_trim(stk_qty(k)) == 0) stk_qty(k) = stk_qty(1)       ! default: same quantity as set 1
      if (len_trim(stk_qty(k)) == 0) stk_qty(k) = 'temp'
      call qty_parse(stk_qty(k), qtmp) ; need(1:8) = need(1:8) .or. qty_needs(qtmp)
    end if
  end do
  do k = 1, ms
    if (len_trim(sel_plane(k)) > 0) then
      if (len_trim(sel_qty(k)) == 0) call fatal('sel_qty('//trim(zpad(k,1))//') is not given')
      call qty_parse(sel_qty(k), qtmp) ; need(1:8) = need(1:8) .or. qty_needs(qtmp)
    end if
  end do
  do k = 1, mi
    if (len_trim(iso_qty(k)) > 0) then
      call qty_parse(iso_qty(k), qtmp) ; need(1:8) = need(1:8) .or. qty_needs(qtmp)
    end if
  end do
  if (len_trim(arr_plane) > 0) then
    if (arr_vec == 'v') then
      need(VAR_VX) = .true. ; need(VAR_VY) = .true. ; need(VAR_VZ) = .true.
    end if
  end if
  vlist = pack([(k, k=1,32)], need)

  nf = count(len_trim(datafile) > 0)
  if (nf == 0) nf = nfiles
  if (nf < 1) call fatal('no input file: nfiles must be >= 1 (or give datafile(1), ...)')
  write(*,'(a)') '----- mag3d -----'
  write(*,'(a,3i6,a,l1)') ' grid: ', nx2, ny2, nz2, '  half_domain=', half_domain
  write(*,'(a,3f9.2)')    ' planet center (i0,j0,k0): ', fm%g%c
  call units_report()
  write(*,'(a,i0,a,f6.2,a,*(1x,i0))') ' stride ', stride, '  r_fine ', r_fine, '  variables:', vlist

  do ifile = 1, nf
    if (len_trim(datafile(ifile)) > 0) then
      fni = trim(datafile(ifile))
    else
      fni = trim(make_fname(diri, prefix, ntap + (ifile-1)*nstep, ndigit, suffix))
    end if
    if (len_trim(title(ifile)) > 0) then
      ttl = title(ifile)
    else if (len_trim(title_head) > 0) then
      write(ttl, t_fmt) t0 + real(ifile-1)*dt
      ttl = adjustl(ttl)
      if (ttl(1:1) == '.') ttl = '0'//ttl(1:len(ttl)-1)
      if (ttl(1:2) == '-.') ttl = '-0'//ttl(2:)
      ttl = trim(title_head)//' '//trim(ttl)//trim(t_unit)
    else
      ttl = fni(index(fni,'/',back=.true.)+1:)
    end if
    write(*,'(a)') ' input: '//fni
    ta = wtime()
    call seq_open_read(fin, fni, 4_int64*nx2*ny2)
    call seq_check_size(fin, int(nvar,int64)*nz2)
    call load_block(fin, fm%g, vlist, stride, box_lo, box_hi, fm%coarse)
    fm%has_fine = (r_fine > 0.0 .and. stride > 1)
    if (fm%has_fine) call load_block(fin, fm%g, vlist, 1, [-r_fine,-r_fine,-r_fine], &
                                     [r_fine,r_fine,r_fine], fm%fine)
    call seq_close(fin)
    write(*,'(a,f8.1,a,3i6)') '  read ', wtime()-ta, ' s   coarse grid', fm%coarse%n

    call scene_reset(sc)
    sc%title = ttl ; sc%bg = bg
    tb = wtime()
    if (draw_frame)  call make_frame()
    if (draw_planet) call make_planet()
    if (draw_inner)  call add_sphere(sc, 'inner_boundary', [0.,0.,0.], r_in, [0.7,0.7,0.7,0.25], 24, 48)
    call make_fieldlines()
    write(*,'(a,f8.1,a)') '  field lines ', wtime()-tb, ' s' ; tb = wtime()
    call make_planes()
    call make_isosurfaces()
    write(*,'(a,f8.1,a)') '  isosurfaces ', wtime()-tb, ' s'
    call make_arrows()
    if (colorbar_3d) call make_colorbars_3d()

    if (nf > 1) then
      ob = trim(outbase)//'_'//trim(zpad(ifile,3))
    else
      ob = trim(outbase)
    end if
    if (out_glb) then
      call write_glb(sc, ob//'.glb') ; write(*,'(a)') '  wrote '//ob//'.glb'
    end if
    if (out_x3d) then
      call write_x3d(sc, ob//'.x3d') ; write(*,'(a)') '  wrote '//ob//'.x3d'
    end if
    write(*,'(a,f8.1,a)') '  total ', wtime()-ta, ' s'
  end do
  write(*,'(a)') ' done.'

contains


  ! stop if variables needed for a quantity were not read (prevents drawing zeros)
  subroutine check_loaded(q)
    type(qty_t), intent(in) :: q
    character(len=64) :: s
    if (qty_loaded(q, fm%coarse)) return
    write(s,'(*(1x,i0))') fm%coarse%var
    call fatal('internal error: variables for "'//trim(q%name)//'" were not read (read:'//trim(s)// &
               '). Please report the compiler and its version.')
  end subroutine check_loaded

  ! sim -> display
  pure function dp(p) result(q)
    real, intent(in) :: p(3)
    real :: q(3)
    q = [dsgn*p(1), dsgn*p(2), p(3)]
  end function dp

  ! all variables at one point (only those needed)
  logical function vals_at(p, v)
    real, intent(in) :: p(3)
    real, intent(out) :: v(32)
    real :: w(size(vlist))
    v = 0.0
    vals_at = fsample(fm, p, vlist, w)
    if (vals_at) v(vlist) = w
  end function vals_at

  !============================================================================
  ! frame, axes, labels
  !============================================================================
  subroutine make_frame()
    real :: lo(3), hi(3), c(3,8), e(3,2)
    integer :: k, ed(2,12)
    real, parameter :: fc(4) = [0.7, 0.7, 0.7, 1.0]
    lo = fm%coarse%org
    hi = fm%coarse%org + fm%coarse%dx*real(fm%coarse%n - 1)
    if (half_domain) lo(3) = -hi(3)
    lo = max(lo, fm%g%xlo) ; hi = min(hi, fm%g%xhi)
    do k = 1, 8
      c(:,k) = [merge(hi(1),lo(1),btest(k-1,0)), merge(hi(2),lo(2),btest(k-1,1)), &
                merge(hi(3),lo(3),btest(k-1,2))]
    end do
    ed = reshape([1,2, 3,4, 5,6, 7,8, 1,3, 2,4, 5,7, 6,8, 1,5, 2,6, 3,7, 4,8], [2,12])
    do k = 1, 12
      e(:,1) = dp(c(:,ed(1,k))) ; e(:,2) = dp(c(:,ed(2,k)))
      call add_polyline(sc, 'frame', 2, e, reshape(fc, [4,1]))
    end do
    ! axes through the planet
    do k = 1, 3
      e = 0.0 ; e(k,1) = lo(k) ; e(k,2) = hi(k)
      e(:,1) = dp(e(:,1)) ; e(:,2) = dp(e(:,2))
      call add_polyline(sc, 'axes', 2, e, reshape([0.5,0.5,0.5,1.0], [4,1]))
    end do
    ! labels (display coordinates)
    call axis_labels(1, lo(1), hi(1)) ; call axis_labels(2, lo(2), hi(2)) ; call axis_labels(3, lo(3), hi(3))
    call add_label(sc, trim(sc%title), dp([0.5*(lo(1)+hi(1)), 0.0, 0.0]) + [0.0, 0.0, hi(3)+3.0*label_size], &
                   label_size*1.2, [1.0,1.0,1.0])
  end subroutine make_frame

  subroutine axis_labels(ax, a, b)
    integer, intent(in) :: ax
    real, intent(in) :: a, b
    real :: p(3)
    character(len=16) :: s
    character(len=1) :: nm
    nm = merge('XYZ'(ax:ax), 'xyz'(ax:ax), coord == 'jso')
    p = 0.0 ; p(ax) = a
    write(s,'(i0)') nint(merge(dsgn, 1.0, ax <= 2)*a)
    call add_label(sc, trim(s), dp(p), label_size, [0.8,0.8,0.8])
    p = 0.0 ; p(ax) = b
    write(s,'(i0)') nint(merge(dsgn, 1.0, ax <= 2)*b)
    call add_label(sc, trim(s), dp(p), label_size, [0.8,0.8,0.8])
    ! axis name at the positive end
    p = 0.0
    if (merge(dsgn, 1.0, ax <= 2) > 0.0) then
      p(ax) = b
    else
      p(ax) = a
    end if
    p = dp(p) ; p(3) = p(3) + 1.5*label_size
    call add_label(sc, nm//' ['//trim(ulen())//']', p, label_size, [1.0,1.0,0.6])
  end subroutine axis_labels

  !============================================================================
  ! planet (textured sphere) and Saturn's rings
  !============================================================================
  subroutine make_planet()
    integer :: it
    character(len=512) :: path
    it = 0
    path = planet_texture
    if (trim(path) == 'auto') path = find_texture(trim(UC%planet))
    if (len_trim(path) > 0 .and. trim(path) /= 'none') then
      it = add_image_file(sc, trim(path))
      if (it == 0) write(*,'(a)') '  WARNING: cannot read planet image '//trim(path)//' (drawn in a single colour)'
    end if
    if (it > 0) then
      call add_textured_sphere(sc, 'planet', [0.,0.,0.], 1.0, 48, 96, it, planet_lon0 + merge(0.0, 180.0, dsgn < 0.0))
    else
      call add_sphere(sc, 'planet', [0.,0.,0.], 1.0, [0.25,0.45,1.0,1.0], 24, 48)
    end if
    if (draw_rings .and. trim(UC%planet) == 'saturn') call make_rings()
  end subroutine make_planet

  ! look for ../textures/<name>.jpg (.png) relative to the executable
  function find_texture(name) result(path)
    character(len=*), intent(in) :: name
    character(len=512) :: path, exe, dir
    character(len=4), parameter :: ext(2) = ['.jpg', '.png']
    logical :: ex
    integer :: k, i
    path = ''
    call get_command_argument(0, exe)
    i = index(exe, '/', back=.true.)
    dir = '.'
    if (i > 1) dir = exe(1:i-1)
    do k = 1, 2
      path = trim(dir)//'/../textures/'//name//trim(ext(k))
      inquire(file=path, exist=ex)
      if (ex) return
      path = 'textures/'//name//trim(ext(k))
      inquire(file=path, exist=ex)
      if (ex) return
    end do
    path = ''
  end function find_texture

  ! Saturn's rings: an annulus in the equatorial plane. Radial brightness/opacity mimic the main rings (C, B, Cassini division, A)
  subroutine make_rings()
    integer, parameter :: nr = 64, nph = 128
    real, allocatable :: p(:,:), nn(:,:), col(:,:)
    integer, allocatable :: tri(:,:)
    integer :: i, j, k, im
    real :: r, ph, a, br
    allocate(p(3,(nr+1)*(nph+1)), nn(3,(nr+1)*(nph+1)), col(4,(nr+1)*(nph+1)), tri(3,2*nr*nph))
    k = 0
    do i = 0, nr
      r = ring_r1 + (ring_r2 - ring_r1)*real(i)/real(nr)
      call ring_profile(r, br, a)
      do j = 0, nph
        ph = 2.0*pi*real(j)/real(nph)
        k = k + 1
        p(:,k) = [r*cos(ph), r*sin(ph), 0.0] ; nn(:,k) = [0.0, 0.0, 1.0]
        col(:,k) = [0.93*br, 0.86*br, 0.72*br, a]
      end do
    end do
    k = 0
    do i = 0, nr-1
      do j = 0, nph-1
        k = k + 1 ; tri(:,k) = [i*(nph+1)+j+1, (i+1)*(nph+1)+j+1, (i+1)*(nph+1)+j+2]
        k = k + 1 ; tri(:,k) = [i*(nph+1)+j+1, (i+1)*(nph+1)+j+2, i*(nph+1)+j+2]
      end do
    end do
    im = add_mesh_begin(sc, 'planet_rings', [1.0, 1.0, 1.0, 1.0], .false.)
    call add_triangles(sc, im, size(p,2), p, nn, col, size(tri,2), tri)
  end subroutine make_rings

  pure subroutine ring_profile(r, br, a)
    real, intent(in) :: r
    real, intent(out) :: br, a
    if (r < 1.527) then            ! C ring
      br = 0.55 ; a = 0.25
    else if (r < 1.95) then        ! B ring
      br = 1.0 ; a = 0.85
    else if (r < 2.025) then       ! Cassini division
      br = 0.4 ; a = 0.08
    else if (r < 2.214 .or. r > 2.226) then   ! A ring (except the Encke gap)
      br = 0.85 ; a = 0.6
    else
      br = 0.4 ; a = 0.1
    end if
  end subroutine ring_profile

  !============================================================================
  ! colour bar: 33 colours of a palette
  !============================================================================
  function pal33(ico) result(c)
    integer, intent(in) :: ico
    real :: c(3,33)
    integer :: k
    do k = 1, 33
      call palette_rgb(real(k-1)/32.0, ico, c(1,k), c(2,k), c(3,k))
    end do
  end function pal33

  ! colour bars inside the 3-D scene (placed horizontally below the frame)
  subroutine make_colorbars_3d()
    real :: lo(3), hi(3), w, hb, x0, z0, img(4,128,4), cn(3,4)
    integer :: k, i
    character(len=96) :: s
    if (sc%nleg == 0) return
    lo = max(fm%coarse%org, fm%g%xlo)
    hi = min(fm%coarse%org + fm%coarse%dx*real(fm%coarse%n - 1), fm%g%xhi)
    if (half_domain) lo(3) = -hi(3)
    w = 0.4*(hi(1) - lo(1)) ; hb = 1.2*label_size
    do k = 1, sc%nleg
      do i = 1, 128
        img(1:3,i,:) = spread(sc%leg(k)%rgb(:, 1 + nint(32.0*real(i-1)/127.0)), 2, 4)
        img(4,i,:) = 1.0
      end do
      x0 = 0.5*(lo(1) + hi(1)) - 0.5*w
      z0 = lo(3) - real(k)*3.5*label_size
      cn(:,1) = dp([x0, lo(2), z0])     ; cn(:,2) = dp([x0+w, lo(2), z0])
      cn(:,3) = dp([x0+w, lo(2), z0+hb]) ; cn(:,4) = dp([x0, lo(2), z0+hb])
      call add_image_quad(sc, 'colorbar'//trim(zpad(k,2)), cn, 128, 4, img)
      write(s,'(es9.2)') sc%leg(k)%vmin
      call add_label(sc, trim(adjustl(s)), dp([x0, lo(2), z0 - 0.8*label_size]), 0.8*label_size, [1.,1.,1.])
      write(s,'(es9.2)') sc%leg(k)%vmax
      call add_label(sc, trim(adjustl(s)), dp([x0+w, lo(2), z0 - 0.8*label_size]), 0.8*label_size, [1.,1.,1.])
      call add_label(sc, trim(sc%leg(k)%label), dp([x0+0.5*w, lo(2), z0 + hb + 0.9*label_size]), &
                     0.8*label_size, [1.,1.,0.6])
    end do
  end subroutine make_colorbars_3d

  !============================================================================
  ! field lines
  !============================================================================
  subroutine make_fieldlines()
    integer :: nh, ih, ilat, imlt, nstart, n, iend1, iend2, topo, np1, np2, k, j
    real, allocatable :: sp(:,:), pa(:,:), pb(:,:), dip(:,:)
    integer, allocatable :: tp(:,:,:), kind(:)
    type(lines_store), allocatable :: res(:)
    real :: lat, th0, ph, ell, rm, thm, hs, b(3), dir
    logical :: ok
    ! ---- start points
    nstart = 0
    nh = 0
    if (fl_polar) nh = merge(2, 1, fl_hemi == 'both')
    n = nh*fl_nlat*fl_nmlt + eq_nx*eq_ny + fl_npt + 1
    if (fl_box) n = n + product(max(fl_box_n, 1))
    allocate(sp(3, n), kind(n))
    do ih = 1, nh
      hs = 1.0
      if (fl_hemi == 'south' .or. (fl_hemi == 'both' .and. ih == 2)) hs = -1.0
      do ilat = 1, fl_nlat
        lat = fl_lat_min + (fl_lat_max - fl_lat_min)*real(ilat-1)/real(max(fl_nlat-1,1))
        th0 = (90.0 - lat)*pi/180.0
        ell = 1.0/sin(th0)**2
        rm = min(fl_r_start, ell)
        thm = asin(min(sqrt(rm/ell), 1.0))
        do imlt = 1, fl_nmlt
          ph = 2.0*pi*real(imlt-1)/real(fl_nmlt) - pi       ! angle from noon (from MLT 0)
          nstart = nstart + 1
          sp(:,nstart) = [-rm*sin(thm)*cos(ph), -rm*sin(thm)*sin(ph), hs*rm*cos(thm)]
          kind(nstart) = ih
        end do
      end do
    end do
    if (fl_equator) then
      do j = 1, eq_ny
        do k = 1, eq_nx
          nstart = nstart + 1
          sp(:,nstart) = [eq_xmin + (eq_xmax-eq_xmin)*real(k-1)/real(max(eq_nx-1,1)), &
                          eq_ymin + (eq_ymax-eq_ymin)*real(j-1)/real(max(eq_ny-1,1)), 0.0]
          kind(nstart) = 0
        end do
      end do
    end if
    do k = 1, fl_npt
      nstart = nstart + 1 ; sp(:,nstart) = [fl_px(k), fl_py(k), fl_pz(k)] ; kind(nstart) = 0
    end do
    if (fl_box) then
      do ilat = 1, max(fl_box_n(3), 1)
        do j = 1, max(fl_box_n(2), 1)
          do k = 1, max(fl_box_n(1), 1)
            nstart = nstart + 1
            sp(:,nstart) = [lerp(1, k), lerp(2, j), lerp(3, ilat)]
            kind(nstart) = 0
          end do
        end do
      end do
    end if
    if (nstart == 0) return
    allocate(res(nstart))
    ! ---- trace (parallel)
    !$omp parallel do schedule(dynamic,1) private(k,pa,pb,np1,np2,iend1,iend2,b,dir,ok,topo,dip)
    do k = 1, nstart
      res(k)%n = 0
      ok = fsample(fm, sp(:,k), [VAR_BX, VAR_BY, VAR_BZ], b)
      if (.not. ok) cycle
      if (sum(b**2) <= 0.0) cycle
      if (kind(k) > 0) then
        ! from the polar region: one line outward
        dir = sign(1.0, sum(b*sp(:,k)))
        call trace_line(fm, sp(:,k), dir, fl_ds, fl_maxlen, r_in, fl_seg, pa, np1, iend1)
        topo = merge(0, merge(1, 2, iend1 == END_OUT), iend1 == END_PLANET)
        if (fl_dipole_part) then
          call dipole_segment(sp(:,k), dip)
          call store(res(k), reshape([dip, pa(:,1:np1)], [3, size(dip,2)+np1]), topo)
        else
          call store(res(k), pa(:,1:np1), topo)
        end if
      else
        ! from the equator / points: both directions
        call trace_line(fm, sp(:,k),  1.0, fl_ds, fl_maxlen, r_in, fl_seg, pa, np1, iend1)
        call trace_line(fm, sp(:,k), -1.0, fl_ds, fl_maxlen, r_in, fl_seg, pb, np2, iend2)
        if (iend1 == END_PLANET .and. iend2 == END_PLANET) then
          topo = 0
        else if (iend1 == END_OUT .and. iend2 == END_OUT) then
          topo = 3
        else if ((iend1 == END_PLANET .and. iend2 == END_OUT) .or. &
                 (iend1 == END_OUT .and. iend2 == END_PLANET)) then
          topo = 1
        else
          topo = 2
        end if
        if (keep_line(topo)) then
          if (iend1 == END_LOOP) then
            call store(res(k), pa(:,1:np1), topo)        ! closed loops in one direction only
          else
            call store(res(k), reshape([pb(:,np2:1:-1), pa(:,2:np1)], [3, np2+np1-1]), topo)
          end if
        end if
      end if
    end do
    !$omp end parallel do
    ! ---- selection near the open-closed boundary (polar start only)
    if (fl_select == 'ocb' .and. nh > 0) then
      allocate(tp(fl_nlat, fl_nmlt, nh))
      n = 0
      do ih = 1, nh
        do ilat = 1, fl_nlat
          do imlt = 1, fl_nmlt
            n = n + 1 ; tp(ilat,imlt,ih) = res(n)%topo
            if (res(n)%n == 0) tp(ilat,imlt,ih) = -1
          end do
        end do
      end do
      n = 0
      do ih = 1, nh
        do ilat = 1, fl_nlat
          do imlt = 1, fl_nmlt
            n = n + 1
            if (.not. near_ocb(tp(:,:,ih), ilat, imlt)) res(n)%n = 0
          end do
        end do
      end do
    end if
    ! ---- add to the scene
    do k = 1, nstart
      if (res(k)%n < 2) cycle
      call add_line(res(k))
    end do
    write(*,'(a,i0,a,4(1x,i0))') '  field lines: ', count(res(:)%n >= 2), &
         '  closed/open/detached/solar-wind:', (count(res(:)%n >= 2 .and. res(:)%topo == j), j=0,3)
    if (fl_color /= 'topology') call add_legend(sc, 'fieldlines', 'field lines: '//trim(fl_color)// &
         ' ['//trim(qty_unit(qfl))//']', fl_vmin, fl_vmax, pal33(fl_palette))
  end subroutine make_fieldlines

  ! coordinates of the 3-D grid seeds
  real function lerp(ax, i)
    integer, intent(in) :: ax, i
    if (fl_box_n(ax) <= 1) then
      lerp = 0.5*(fl_box_lo(ax) + fl_box_hi(ax))
    else
      lerp = fl_box_lo(ax) + (fl_box_hi(ax) - fl_box_lo(ax))*real(i-1)/real(fl_box_n(ax)-1)
    end if
  end function lerp

  ! selection by fl_keep (topo: 0 closed, 1 open, 2 loop/undetermined, 3 solar wind)
  logical function keep_line(topo)
    integer, intent(in) :: topo
    select case (trim(fl_keep))
    case ('unconnected') ; keep_line = (topo == 2 .or. topo == 3)
    case ('connected')   ; keep_line = (topo == 0 .or. topo == 1)
    case ('closed')      ; keep_line = (topo == 0)
    case ('open')        ; keep_line = (topo == 1)
    case ('solarwind')   ; keep_line = (topo == 3)
    case ('detached')    ; keep_line = (topo == 2)
    case default         ; keep_line = .true.
    end select
  end function keep_line

  logical function near_ocb(t, i, j)
    integer, intent(in) :: t(:,:), i, j
    integer :: di, dj, ii, jj, nm
    near_ocb = .false.
    if (t(i,j) < 0 .or. t(i,j) > 1) return
    nm = size(t,2)
    do di = -1, 1
      do dj = -1, 1
        ii = i + di ; jj = mod(j + dj - 1 + nm, nm) + 1
        if (ii < 1 .or. ii > size(t,1)) cycle
        if (t(ii,jj) >= 0 .and. t(ii,jj) <= 1 .and. t(ii,jj) /= t(i,j)) near_ocb = .true.
      end do
    end do
  end function near_ocb

  ! dipole field line from the surface (r=1) to p (same longitude, r = L sin^2)
  subroutine dipole_segment(p, d)
    real, intent(in) :: p(3)
    real, allocatable, intent(out) :: d(:,:)
    real :: r, rho, th, ell, th1, th2, t, phi, rr, hs
    integer :: n, k
    r = sqrt(sum(p**2)) ; rho = sqrt(p(1)**2 + p(2)**2)
    hs = sign(1.0, p(3))
    th = atan2(rho, abs(p(3)))                     ! angle from the pole
    ell = r/max(sin(th)**2, 1.0e-6)
    th1 = asin(min(1.0/sqrt(ell), 1.0)) ; th2 = th
    phi = atan2(p(2), p(1))
    n = max(2, nint(r/max(fl_seg,0.2)) + 2)
    allocate(d(3, n-1))
    do k = 1, n-1                                   ! the last point (= p) belongs to the traced part
      t = th1 + (th2 - th1)*real(k-1)/real(n-1)
      rr = ell*sin(t)**2
      d(:,k) = [rr*sin(t)*cos(phi), rr*sin(t)*sin(phi), hs*rr*cos(t)]
    end do
  end subroutine dipole_segment

  subroutine store(r, p, topo)
    type(lines_store), intent(inout) :: r
    real, intent(in) :: p(:,:)
    integer, intent(in) :: topo
    r%n = size(p,2) ; r%topo = topo
    r%p = p
  end subroutine store

  subroutine add_line(r)
    type(lines_store), intent(in) :: r
    real, allocatable :: q(:,:), c(:,:)
    real :: v(32), x, cl(3), lo, hi
    integer :: k
    character(len=24) :: grp
    allocate(q(3, r%n), c(4, r%n))
    do k = 1, r%n
      q(:,k) = dp(r%p(:,k))
    end do
    select case (r%topo)
    case (0) ; cl = col_closed ; grp = 'fieldlines_closed'
    case (1) ; cl = col_open   ; grp = 'fieldlines_open'
    case (3) ; cl = col_sw     ; grp = 'fieldlines_solarwind'
    case default ; cl = col_undet ; grp = 'fieldlines_detached'
    end select
    if (fl_color == 'topology') then
      c(1:3,:) = spread(cl, 2, r%n) ; c(4,:) = 1.0
    else
      lo = fl_vmin ; hi = fl_vmax
      if (lo == hi) then
        lo = 0.0 ; hi = 1.0
      end if
      do k = 1, r%n
        x = 0.0
        if (vals_at(r%p(:,k), v)) x = (qty_eval(qfl, v) - lo)/(hi - lo)
        call palette_rgb(min(max(x,0.0),1.0), fl_palette, c(1,k), c(2,k), c(3,k))
        c(4,k) = 1.0
      end do
    end if
    call add_polyline(sc, trim(grp), r%n, q, c)
  end subroutine add_line

  !============================================================================
  ! colour planes (pixel images)
  !============================================================================
  subroutine make_planes()
    integer :: k, s
    real :: z
    character(len=64) :: nm
    character(len=160) :: js
    do k = 1, mp
      if (len_trim(pl_plane(k)) == 0) cycle
      call one_plane(pl_plane(k), pl_pos(k), pl_qty(k), pl_vmin(k), pl_vmax(k), pl_palette(k), &
                     merge(1.0, pl_alpha(k), pl_alpha_mode(k) == 'value'), pl_alpha_mode(k), &
                     pl_tra(k), pl_h_lo(k), pl_h_hi(k), &
                     pl_v_lo(k), pl_v_hi(k), 'plane'//trim(zpad(k,2))//'_'//trim(pl_qty(k)), &
                     'plane'//trim(zpad(k,2)))
    end do
    ! ---- stacked slices (up to 3 sets, common colour range within a set)
    do s = 1, mk
      if (len_trim(stk_plane(s)) == 0) cycle
      if (stk_vmin(s) == stk_vmax(s)) call qty_range(stk_qty(s), stk_vmin(s), stk_vmax(s))
      do k = 1, stk_n(s)
        z = stk_lo(s) + (stk_hi(s) - stk_lo(s))*real(k-1)/real(max(stk_n(s)-1,1))
        nm = 'stack'//trim(zpad(s,1))//'_'//trim(stk_plane(s))//'_'//trim(zpad(k,2))
        call one_plane(stk_plane(s), z, stk_qty(s), stk_vmin(s), stk_vmax(s), stk_palette(s), stk_alpha(s), 'value', &
                       stk_tra(s), -1.0e30, 1.0e30, -1.0e30, 1.0e30, nm, 'stack'//trim(zpad(s,1))//'_')
      end do
    end do
    ! ---- selectable slices (one is chosen in the viewer)
    do s = 1, ms
      if (len_trim(sel_plane(s)) == 0) cycle
      if (sel_vmin(s) == sel_vmax(s)) call qty_range(sel_qty(s), sel_vmin(s), sel_vmax(s))
      do k = 1, sel_n(s)
        z = sel_lo(s) + (sel_hi(s) - sel_lo(s))*real(k-1)/real(max(sel_n(s)-1,1))
        nm = 'sel'//trim(zpad(s,1))//'_'//trim(sel_plane(s))//'_'//trim(zpad(k,2))
        call one_plane(sel_plane(s), z, sel_qty(s), sel_vmin(s), sel_vmax(s), sel_palette(s), &
                       sel_alpha(s), 'const', 1.0, -1.0e30, 1.0e30, -1.0e30, 1.0e30, nm, &
                       'sel'//trim(zpad(s,1))//'_')
        write(js,'(a,i0,a,i0,a,i0,a,a,a,g0.5,a)') '{"group":"sel', s, '","index":', k-1, ',"n":', &
             sel_n(s), ',"axis":"', trim(sel_axis(sel_plane(s))), '","pos":', &
             merge(dsgn, 1.0, sel_plane(s) /= 'xy')*z + 0.0, '}'
        call set_extras(sc, nm, trim(js))
      end do
    end do
  end subroutine make_planes

  ! name of the axis normal to the slice (display coordinates)
  function sel_axis(pl) result(a)
    character(len=*), intent(in) :: pl
    character(len=1) :: a
    select case (pl)
    case ('xy') ; a = merge('Z', 'z', coord == 'jso')
    case ('xz') ; a = merge('Y', 'y', coord == 'jso')
    case default; a = merge('X', 'x', coord == 'jso')
    end select
  end function sel_axis

  ! range of a quantity on the grid points of the global field
  subroutine qty_range(qn, a, b)
    character(len=*), intent(in) :: qn
    real, intent(out) :: a, b
    type(qty_t) :: q
    real(real32), allocatable :: f(:,:,:)
    integer :: sv
    call qty_parse(qn, q)
    sv = ist ; ist = 1
    call node_field(q, .false., f)
    ist = sv
    a = minval(f, mask=.not. ieee_is_nan(f)) ; b = maxval(f, mask=.not. ieee_is_nan(f))
    if (b <= a) b = a + 1.0
  end subroutine qty_range


  subroutine one_plane(pl, pos, qn, vmin, vmax, ico, alpha, amode, tra, hlo, hhi, vlo, vhi, name, legt)
    character(len=*), intent(in) :: pl, qn, amode, name, legt
    real, intent(in) :: pos, vmin, vmax, alpha, tra, hlo, hhi, vlo, vhi
    integer, intent(in) :: ico
    type(qty_t) :: q
    integer :: ah, av, an, nw, nh, i, j
    real :: lo(3), hi(3), a0, a1, b0r, b1, p(3), v(32), x, lo_v, hi_v, cn(3,4), dxp, dmin, dmax
    integer :: nvalid
    real, allocatable :: val(:,:), img(:,:,:)
    character(len=128) :: lab
    call qty_parse(qn, q)
    call check_loaded(q)
    select case (pl)
    case ('xy') ; ah = 1 ; av = 2 ; an = 3
    case ('xz') ; ah = 1 ; av = 3 ; an = 2
    case default; ah = 2 ; av = 3 ; an = 1
    end select
    lo = max(fm%coarse%org, fm%g%xlo)
    hi = min(fm%coarse%org + fm%coarse%dx*real(fm%coarse%n - 1), fm%g%xhi)
    if (half_domain) lo(3) = -hi(3)
    a0 = max(lo(ah), hlo) ; a1 = min(hi(ah), hhi)
    b0r = max(lo(av), vlo) ; b1 = min(hi(av), vhi)
    dxp = fm%coarse%dx
    nw = min(2048, max(2, nint((a1-a0)/dxp) + 1))
    nh = min(2048, max(2, nint((b1-b0r)/dxp) + 1))
    allocate(val(nw,nh), img(4,nw,nh))
    !$omp parallel do private(i,p,v,x)
    do j = 1, nh
      do i = 1, nw
        p = 0.0 ; p(an) = pos
        p(ah) = a0 + (a1-a0)*(real(i)-0.5)/real(nw)
        p(av) = b0r + (b1-b0r)*(real(j)-0.5)/real(nh)
        if (vals_at(p, v)) then
          val(i,j) = qty_eval(q, v)
        else
          val(i,j) = nan32()
        end if
      end do
    end do
    !$omp end parallel do
    dmin = minval(val, mask=.not. ieee_is_nan(val)) ; dmax = maxval(val, mask=.not. ieee_is_nan(val))
    nvalid = count(.not. ieee_is_nan(val))
    lo_v = vmin ; hi_v = vmax
    if (lo_v == hi_v) then
      lo_v = dmin ; hi_v = dmax
      if (hi_v <= lo_v) hi_v = lo_v + 1.0
    end if
    do j = 1, nh
      do i = 1, nw
        if (ieee_is_nan(val(i,j))) then
          img(:,i,j) = 0.0 ; cycle
        end if
        x = min(max((val(i,j) - lo_v)/(hi_v - lo_v), 0.0), 1.0)
        call palette_rgb(x, ico, img(1,i,j), img(2,i,j), img(3,i,j))
        if (amode == 'value') then
          img(4,i,j) = alpha*min(1.0, min(x*tra, 0.5*tra))     ! alpha: maximum opacity
        else
          img(4,i,j) = alpha
        end if
      end do
    end do
    ! corners: LL, LR, UR, UL
    do i = 1, 4
      p = 0.0 ; p(an) = pos
      p(ah) = merge(a1, a0, i == 2 .or. i == 3)
      p(av) = merge(b1, b0r, i >= 3)
      cn(:,i) = dp(p)
    end do
    call add_image_quad(sc, name, cn, nw, nh, img)
    call add_legend(sc, legt, trim(qn)//' ['//trim(qty_unit(q))//']', lo_v, hi_v, pal33(ico))
    write(*,'(a,a,a,a,a,2es11.3,a,2es11.3)') '  plane ', trim(name), ' [', trim(qty_unit(q)), &
         ']  color range', lo_v, hi_v, '   data range', dmin, dmax
    if (plane_labels .and. index(name, 'plane') == 1) then
      write(lab,'(a,a,a,a,es9.2,a,es9.2)') trim(qn), ' [', trim(qty_unit(q)), ']  ', lo_v, ' - ', hi_v
      call add_label(sc, trim(lab), cn(:,4) + [0.0, 0.0, 0.8*label_size], 0.7*label_size, &
                     [1.0, 1.0, 1.0])
    end if
    if (nvalid == 0) then
      write(*,'(a)') '    WARNING: no data on this plane (outside the data / box region?)'
    else if (dmax <= dmin) then
      write(*,'(a)') '    WARNING: the quantity is constant on this plane (variable is 0 in the file?)'
    else if (dmax < lo_v .or. dmin > hi_v) then
      write(*,'(a)') '    WARNING: data range is outside vmin..vmax (the plane is drawn in a single colour)'
    end if
  end subroutine one_plane

  !============================================================================
  ! isosurfaces (marching tetrahedra on the coarse grid)
  !============================================================================
  subroutine make_isosurfaces()
    integer :: k, ims, im, nh, nv, nt
    integer(int64) :: ntot
    type(qty_t) :: q
    real(real32), allocatable :: s(:,:,:)
    real, allocatable :: pv(:,:), nv3(:,:)
    integer, allocatable :: tri(:,:)
    nh = merge(2, 1, half_domain)
    do k = 1, mi
      if (len_trim(iso_qty(k)) == 0) cycle
      call qty_parse(iso_qty(k), q)
      call check_loaded(q)
      im = add_mesh_begin(sc, 'iso'//trim(zpad(k,2))//'_'//trim(iso_qty(k)), &
                          [1.0, 1.0, 1.0, iso_alpha(k)], .true.)
      ! more than iso_max_tri triangles -> increase the thinning and extract again (count ~ 1/thinning^2)
      ist = max(1, iso_stride)
      do
        call node_field(q, .false., s)
        call smooth3(s, iso_smooth)
        call iso_extract(s, fm%coarse%org, fm%coarse%dx*real(ist), iso_level(k), nv, pv, nv3, nt, tri)
        ntot = int(nt, int64)*nh
        if (iso_max_tri <= 0 .or. ntot <= iso_max_tri .or. ist >= 16) exit
        write(*,'(a,i0,a,i0,a,i0)') '    isosurface: ', ntot, ' triangles > iso_max_tri = ', iso_max_tri, &
             ' -> stride x', ist + max(1, ceiling(real(ist)*(sqrt(real(ntot)/real(iso_max_tri)) - 1.0)))
        ist = ist + max(1, ceiling(real(ist)*(sqrt(real(ntot)/real(iso_max_tri)) - 1.0)))
      end do
      call isosurface(nv, pv, nv3, nt, tri, .false., iso_color(:,k), im)
      do ims = 2, nh
        call node_field(q, .true., s)
        call smooth3(s, iso_smooth)
        call iso_extract(s, fm%coarse%org, fm%coarse%dx*real(ist), iso_level(k), nv, pv, nv3, nt, tri)
        call isosurface(nv, pv, nv3, nt, tri, .true., iso_color(:,k), im)
      end do
      write(*,'(a,a,a,es11.3,a,a,a,i0,a,f6.2,a)') '  isosurface ', trim(iso_qty(k)), ' =', iso_level(k), &
           ' ', trim(qty_unit(q)), ' : ', sc%m(im)%nt, ' triangles (grid spacing ', &
           fm%coarse%dx*real(ist), ' '//trim(ulen())//')'
    end do
  end subroutine make_isosurfaces

  ! 1-2-1 smoothing n times in each direction (edges unchanged)
  subroutine smooth3(s, n)
    real(real32), intent(inout) :: s(:,:,:)
    integer, intent(in) :: n
    real(real32), allocatable :: t(:,:,:)
    integer :: it, d, i, j, k, m(3)
    m = shape(s)
    if (n <= 0 .or. any(m < 3)) return
    allocate(t, mold=s)
    do it = 1, n
      do d = 1, 3
        t = s
        !$omp parallel do private(i,j)
        do k = 1, m(3)
          do j = 1, m(2)
            do i = 1, m(1)
              select case (d)
              case (1)
                if (i > 1 .and. i < m(1)) s(i,j,k) = 0.25*t(i-1,j,k) + 0.5*t(i,j,k) + 0.25*t(i+1,j,k)
              case (2)
                if (j > 1 .and. j < m(2)) s(i,j,k) = 0.25*t(i,j-1,k) + 0.5*t(i,j,k) + 0.25*t(i,j+1,k)
              case (3)
                if (k > 1 .and. k < m(3)) s(i,j,k) = 0.25*t(i,j,k-1) + 0.5*t(i,j,k) + 0.25*t(i,j,k+1)
              end select
            end do
          end do
        end do
        !$omp end parallel do
      end do
    end do
  end subroutine smooth3

  ! scalar on coarse grid nodes (mirrored = values seen at z -> -z)
  subroutine node_field(q, mirrored, s)
    type(qty_t), intent(in) :: q
    logical, intent(in) :: mirrored
    real(real32), allocatable, intent(out) :: s(:,:,:)
    integer :: i, j, k, m, st, n(3)
    real :: v(32)
    st = max(1, ist)
    n = (fm%coarse%n - 1)/st + 1
    allocate(s(n(1), n(2), n(3)))
    !$omp parallel do private(i,j,m,v)
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          v = 0.0
          do m = 1, fm%coarse%nv
            v(fm%coarse%var(m)) = fm%coarse%a(1+(i-1)*st, 1+(j-1)*st, 1+(k-1)*st, m)
            if (mirrored .and. fm%coarse%var(m) <= 8) &
              v(fm%coarse%var(m)) = v(fm%coarse%var(m))*real(PARITY8(fm%coarse%var(m)))
          end do
          s(i,j,k) = qty_eval(q, v)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine node_field

  subroutine isosurface(nv, pv, nv3, nt, tri, mirrored, rgb, im)
    integer, intent(in) :: nv, nt, im
    real, intent(inout) :: pv(:,:), nv3(:,:)
    integer, intent(in) :: tri(:,:)
    real, intent(in) :: rgb(3)
    logical, intent(in) :: mirrored
    real, allocatable :: col(:,:)
    integer :: i
    if (nt == 0) return
    allocate(col(4,nv))
    do i = 1, nv
      if (mirrored) then
        pv(3,i) = -pv(3,i) ; nv3(3,i) = -nv3(3,i)
      end if
      pv(:,i) = dp(pv(:,i)) ; nv3(1:2,i) = dsgn*nv3(1:2,i)
      col(:,i) = [rgb(1), rgb(2), rgb(3), 1.0]      ! transparency is set in the material
    end do
    call add_triangles(sc, im, nv, pv, nv3, col, nt, tri)
  end subroutine isosurface

  !============================================================================
  ! arrows on a plane
  !============================================================================
  subroutine make_arrows()
    integer :: ah, av, an, nw, nh, i, j, vv(3)
    real :: lo(3), hi(3), a0, a1, b0r, b1, p(3), v(32), w(3), wm, scl, vmx, q(3,2), c(4), side(3), hd(3,3)
    real :: x, afac
    if (len_trim(arr_plane) == 0) return
    select case (arr_plane)
    case ('xy') ; ah = 1 ; av = 2 ; an = 3
    case ('xz') ; ah = 1 ; av = 3 ; an = 2
    case default; ah = 2 ; av = 3 ; an = 1
    end select
    vv = merge([VAR_VX, VAR_VY, VAR_VZ], [VAR_BX, VAR_BY, VAR_BZ], arr_vec == 'v')
    afac = real(ufac(merge('v', 'b', arr_vec == 'v')))
    lo = max(fm%coarse%org, fm%g%xlo)
    hi = min(fm%coarse%org + fm%coarse%dx*real(fm%coarse%n - 1), fm%g%xhi)
    if (half_domain) lo(3) = -hi(3)
    a0 = max(lo(ah), arr_h_lo) ; a1 = min(hi(ah), arr_h_hi)
    b0r = max(lo(av), arr_v_lo) ; b1 = min(hi(av), arr_v_hi)
    nw = arr_n ; nh = max(2, nint(real(arr_n)*(b1-b0r)/(a1-a0)))
    vmx = arr_vmax
    if (vmx <= 0.0) then
      do j = 1, nh
        do i = 1, nw
          p = 0.0 ; p(an) = arr_pos
          p(ah) = a0 + (a1-a0)*(real(i)-0.5)/real(nw) ; p(av) = b0r + (b1-b0r)*(real(j)-0.5)/real(nh)
          if (vals_at(p, v)) vmx = max(vmx, afac*sqrt(sum(v(vv)**2)))
        end do
      end do
    end if
    if (vmx <= 0.0) return
    write(*,'(a,a,es11.3,1x,a)') '  arrows ', trim(arr_vec)//'  max', vmx, &
         trim(ulabel(merge('v', 'b', arr_vec == 'v')))
    call add_legend(sc, 'arrows', 'arrows |'//trim(arr_vec)//'| ['//trim(ulabel(merge('v', 'b', arr_vec == 'v')))// &
                    ']', 0.0, vmx, pal33(3))
    scl = arr_scale
    if (scl <= 0.0) scl = 0.9*(a1-a0)/real(nw)/vmx
    do j = 1, nh
      do i = 1, nw
        p = 0.0 ; p(an) = arr_pos
        p(ah) = a0 + (a1-a0)*(real(i)-0.5)/real(nw) ; p(av) = b0r + (b1-b0r)*(real(j)-0.5)/real(nh)
        if (.not. vals_at(p, v)) cycle
        w = afac*v(vv) ; wm = sqrt(sum(w**2))
        if (wm <= 0.0) cycle
        if (wm > vmx) w = w*vmx/wm
        x = min(wm/vmx, 1.0)
        call palette_rgb(x, 3, c(1), c(2), c(3)) ; c(4) = 1.0
        q(:,1) = dp(p) ; q(:,2) = dp(p + scl*w)
        call add_polyline(sc, 'arrows', 2, q, reshape(c, [4,1]))
        ! arrow head
        side = 0.0 ; side(an) = 1.0
        side = cross(w, side)
        if (sum(side**2) <= 0.0) side = cross(w, [1.0, 0.0, 0.0])
        side = side/max(sqrt(sum(side**2)), 1.0e-30)*0.12*scl*sqrt(sum(w**2))
        hd(:,1) = dp(p + scl*w*0.7 + side) ; hd(:,2) = q(:,2) ; hd(:,3) = dp(p + scl*w*0.7 - side)
        call add_polyline(sc, 'arrows', 3, hd, reshape(c, [4,1]))
      end do
    end do
  end subroutine make_arrows

  pure function cross(a, b) result(c)
    real, intent(in) :: a(3), b(3)
    real :: c(3)
    c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
  end function cross


end program mag3d_main
