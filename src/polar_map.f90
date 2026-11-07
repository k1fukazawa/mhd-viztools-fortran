!==============================================================================
! polar_map : computation and plots of quantities mapped to the polar ionosphere (FAC, convection, open-closed boundary)
!
!   Replaces prodp3_ju02nt.f (open-closed boundary), prodpc5_ju02nt.f (convection, energy
!   flux) and polar_ju02nt.f (plots) of the legacy fig_h/polar; one program from computation to plots.
!
! Computation (colatitude x MLT grid on the ionosphere R = 1 planetary radius)
!   1. mapping point: from the ionospheric foot point (colatitude th0) along the dipole field line
!      to radius r_m = min(r_map, L) (L = 1/sin^2 th0). r_map is a little outside the inner
!      boundary (7 RJ) (default 10 RJ). Inside it the field is regarded as a dipole.
!   2. field-aligned current (MHD): J = rot(B - B_dip) at the mapping point, mapped to the
!      ionosphere using the conservation of J_par/B along the field line:
!         j_up = (J.B / B^2)_map * B_r,ion        (upward positive, B_r,ion = -2 b0 cos th0)
!   3. convection: v at the mapping point (corotation optionally subtracted) mapped to the ionosphere by the dipole mapping
!         u_th = d th0/dt = -(dL/dt) / (2 L^1.5 cos th0)
!         u_ph = sin th0 * v_ph / (r sin th)
!         dL/dt = v_r/sin^2 th - 2 cos th v_th / sin^3 th
!      electric field E = -u x B_ion; the potential phi from the Poisson equation on the sphere
!      lap phi = -div E solved directly (DFT in azimuth + tridiagonal in colatitude).
!      j_conv = -div(Sigma_P E) (field-aligned current for a uniform conductance Sigma_P = 1)
!   4. open-closed boundary: traced from the mapping point along B (RK2). Closed if it reaches the
!      other hemisphere (the z=0 plane for half-domain data), open if it leaves the domain.
!      Also the field-line length, equatorial crossing distance and the weighted-average
!      "energy flux" proxy <p sqrt(p/rho)> of the legacy code are computed.
!
! Data
!   * reads the m.data format (output of combine3d) directly. For half-domain data (z>=0, symmetric
!     about the equator) z<0 is filled using parity (Bx,By,vz antisymmetric). Full-domain data also work.
!   * near the mapping points (|x|,|y|,|z| <= r_map + a few cells) at full resolution; the global
!     field for tracing (B, rho, p) thinned by stride (standard: stride=3, about 4.5 GB).
!
! Output: PostScript (4 panels as the legacy plot by default) and raw binary of all maps (.bin/.txt)
! Usage:  ./polar_map [namelist file]     (default: polar.nml)
!==============================================================================
program polar_map
  use mhd_common
  use psplot
  use mhdunits
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: maxf = 1000
  real(real64), parameter :: pi = 3.14159265358979_real64
  !---- namelist &polar ---------------------------------------------------------
  character(len=512) :: datafile(maxf) = '', diri = './'
  character(len=64)  :: prefix = 'ju02nt', suffix = 'm.data'
  integer :: ntap = 171, nfiles = 1, nstep = 1, ndigit = 3
  character(len=128) :: title(maxf) = '', title_head = ''
  real :: t0 = 0.0, dt = 0.0
  character(len=16)  :: t_unit = ' h', t_fmt = '(f0.2)'
  character(len=128) :: header = ''      ! empty: chosen from the planet name
  ! grid of the data file
  integer :: nx2 = 3002, ny2 = 2002, nz2 = 1002, nvar = 8, nxp = 500
  real    :: h = 0.15
  real    :: i0 = -1.0, j0 = -1.0, k0 = -1.0
  logical :: half_domain = .true.          ! z>=0 only (symmetric about the equator)
  integer :: iv_rho = 1, iv_v = 2, iv_p = 5, iv_b = 6   ! variable numbers (x components for v, B)
  ! physics
  real    :: b0 = -1.0                     ! dipole (b0 of the simulation)
  real    :: r_in = 7.0                    ! inner boundary radius (weight of the energy flux)
  real    :: omega = 1.38e-4               ! planetary rotation rate (vo0 of the simulation)
  logical :: subtract_corot = .true.       ! subtract rigid corotation from the convection
  real    :: r_map = 10.0                  ! mapping radius [planetary radii]
  real    :: r_valid = 8.0                 ! mapping points inside this are not used
  ! J_par: difference width (multiple of the grid spacing h) and averaging length along the field line
  ! (suppresses amplification of 2-cell oscillations of the MHD grid by rot; J_par/B is conserved along B)
  real    :: jpar_delta = 2.0, jpar_avg = 4.0
  integer :: nsmooth = 2                   ! smoothing passes of the plots (5-point as legacy diff2; not applied to the data)
  ! ionospheric grid
  character(len=8) :: hemisphere = 'south' ! 'north' or 'south'
  character(len=8) :: view = ''            ! 'above' / 'below' (default: north=above, south=below)
  real    :: colat_max = 20.0              ! plotted range (angle from the pole) [deg]
  integer :: ncolat = 80, nmlt = 96
  ! field-line tracing
  logical :: do_trace = .true.
  integer :: stride = 3                    ! thinning of the global field
  real    :: trace_ds = 0.0                ! step (0: 0.5*h*stride)
  real    :: trace_maxlen = 5000.0         ! maximum tracing length
  real    :: box_lo(3) = -1.0e30, box_hi(3) = 1.0e30   ! region of the global field to read
  ! output / plot
  character(len=512) :: psfile = 'polar.ps', outbase = ''
  logical :: ps_per_file = .false.
  real    :: jmhd_max = 0.0, jconv_max = 0.0, eflux_max = 0.0, umax = 0.0
  integer :: nphi_cont = 16
  real    :: jfac = 1.0                    ! extra factor for the currents (units are converted by &units)
  real    :: sigma_p = 1.0                 ! uniform Pedersen conductance for j_conv [S] (1 for normalized units)
  integer :: vec_ncolat = 10, vec_nmlt = 24
  namelist /polar/ datafile, diri, prefix, suffix, ntap, nfiles, nstep, ndigit,  &
       title, title_head, t0, dt, t_unit, t_fmt, header, nx2, ny2, nz2, nvar,    &
       nxp, h, i0, j0, k0, half_domain, iv_rho, iv_v, iv_p, iv_b, b0, r_in,      &
       omega, subtract_corot, r_map, r_valid, hemisphere, view, colat_max,       &
       ncolat, nmlt, do_trace, stride, trace_ds, trace_maxlen, box_lo, box_hi,   &
       psfile, outbase, ps_per_file, jmhd_max, jconv_max, eflux_max, umax,       &
       nphi_cont, jfac, sigma_p, vec_ncolat, vec_nmlt, jpar_delta, jpar_avg, nsmooth

  !---- data blocks -----------------------------------------------------------------
  ! loc : full resolution near the mapping points (8 variables) / glo : thinned field for tracing (Bx,By,Bz,rho,p)
  type :: block
    integer :: n(3) = 0
    real    :: org(3) = 0.0, dx = 1.0      ! position of grid point 1 and spacing
    real(real32), allocatable :: a(:,:,:,:)
  end type block
  type(block) :: loc, glo

  !---- maps on the ionospheric grid (ncolat, nmlt) ----------------------------------
  integer, parameter :: nmap = 12
  ! 1 j_mhd  2 j_conv  3 phi  4 u_th(away from the pole)  5 u_ph  6 topology
  ! 7 eflux  8 length  9 r_eq 10 r_m 11 upot_th 12 upot_ph
  character(len=8), parameter :: mapname(nmap) = [character(len=8) :: &
       'j_mhd', 'j_conv', 'phi', 'u_colat', 'u_mlt', 'topology', 'eflux', &
       'length', 'r_eq', 'r_map', 'up_colat', 'up_mlt']
  real(real32), allocatable :: mp(:,:,:)
  real(real64) :: dth, dph, hs
  real :: xlo(3), xhi(3)
  integer :: u, ios, ifile, nf, ino
  real(real64) :: ta
  real(real32) :: rnan
  real :: fu = 1.0, fphi = 1.0, fj = 1.0, fjc = 1.0, fef = 1.0   ! unit conversion factors
  character(len=:), allocatable :: fni, nmlfile, psname, obase
  character(len=512) :: arg
  character(len=256) :: msg
  character(len=128) :: ttl
  type(seqfile) :: fin

  !---- namelist -------------------------------------------------------------------
  nmlfile = 'polar.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=u, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(u, nml=polar, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &polar : '//trim(msg))
  call units_read(u)
  close(u)
  if (len_trim(header) == 0) header = default_header()
  if (i0 < 0.0) i0 = real(nx2+1)/2.0 - real(nxp)
  if (j0 < 0.0) j0 = real(ny2+1)/2.0
  if (k0 < 0.0) then
    if (half_domain) then
      k0 = 1.5
    else
      k0 = real(nz2+1)/2.0
    end if
  end if
  if (len_trim(view) == 0) view = merge('above', 'below', hemisphere == 'north')
  hs = merge(1.0_real64, -1.0_real64, hemisphere == 'north')
  if (trace_ds <= 0.0) trace_ds = 0.5*h*real(stride)
  rnan = ieee_value(rnan, ieee_quiet_nan)
  dth = real(colat_max,real64)*pi/180.0_real64/real(ncolat,real64)
  dph = 2.0_real64*pi/real(nmlt,real64)
  ! data domain [RJ]
  xlo = [(1.0-i0)*h, (1.0-j0)*h, (1.0-k0)*h]
  xhi = [(real(nx2)-i0)*h, (real(ny2)-j0)*h, (real(nz2)-k0)*h]
  if (half_domain) xlo(3) = -xhi(3)
  nf = count(len_trim(datafile) > 0)
  if (nf == 0) nf = nfiles
  if (nf < 1) call fatal('no input file: nfiles must be >= 1 (or give datafile(1), ...)')
  allocate(mp(ncolat, nmlt, nmap))

  write(*,'(a)') '----- polar_map -----'
  write(*,'(a,3i6,a,l1)') ' grid: ', nx2, ny2, nz2, '  half_domain=', half_domain
  write(*,'(a,3f9.2)')    ' planet center (i0,j0,k0): ', i0, j0, k0
  write(*,'(a,a,a,f6.2,a,f5.1,a)') ' hemisphere: ', trim(hemisphere), '  r_map=', r_map, &
                                   '  colat_max=', colat_max, ' deg'
  write(*,'(a,i0,a,i0)') ' ionospheric grid: ncolat=', ncolat, ' nmlt=', nmlt
  call units_report()
  ! normalized -> display units
  fu   = real(ufac('v'))                 ! ionospheric speed [km/s]
  fphi = real(ufac('phi'))               ! potential [kV]
  fj   = real(ufac('j'))*jfac            ! current density [uA/m^2]
  fjc  = jfac                            ! -div E (Sigma_P = 1)
  if (UC%phys) fjc = real(sigma_p*UC%vs*UC%bs/UC%xs*1.0e6_real64)*jfac   ! Sigma_P div E [uA/m^2]
  fef  = real(ufac('eflux'))             ! [mW/m^2]

  ino = 0
  do ifile = 1, nf
    if (len_trim(datafile(ifile)) > 0) then
      fni = trim(datafile(ifile))
    else
      fni = trim(make_fname(diri, prefix, ntap + (ifile-1)*nstep, ndigit, suffix))
    end if
    write(*,'(a)') ' input: '//fni
    ta = wtime()
    call seq_open_read(fin, fni, 4_int64*nx2*ny2)
    call seq_check_size(fin, int(nvar,int64)*nz2)
    call load_local()
    if (do_trace) call load_global()
    call seq_close(fin)
    write(*,'(a,f8.1,a)') '  read: ', wtime()-ta, ' s'

    call compute_maps()
    call to_units()
    write(*,'(a,f8.1,a)') '  maps: ', wtime()-ta, ' s'
    call map_summary()

    ! ---- outputs
    if (len_trim(outbase) > 0) then
      obase = trim(outbase)//'_'//trim(zpad(ifile,3))
    else
      obase = trim(strip_ext(psfile))//'_'//trim(zpad(ifile,3))
    end if
    call write_maps(obase)
    if (ps_per_file .or. ifile == 1) then
      if (ifile > 1) call ps_close()
      if (ps_per_file) then
        psname = trim(strip_ext(psfile))//'_'//trim(zpad(ifile,3))//'.ps'
      else
        psname = trim(psfile)
      end if
      call ps_open(psname) ; ino = 0
      write(*,'(a)') ' PostScript: '//psname
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
    ino = ino + 1
    call draw_page(ino, ttl)
    write(*,'(a,f8.1,a)') '  total: ', wtime()-ta, ' s'
  end do
  call ps_close()
  write(*,'(a)') ' done.'

contains

  !================================================================================
  ! normalized -> display units (NaN unchanged)
  !================================================================================
  subroutine to_units()
    mp(:,:,1) = mp(:,:,1)*fj
    mp(:,:,2) = mp(:,:,2)*fjc
    mp(:,:,3) = mp(:,:,3)*fphi
    mp(:,:,4:5) = mp(:,:,4:5)*fu
    mp(:,:,7) = mp(:,:,7)*fef
    mp(:,:,11:12) = mp(:,:,11:12)*fu
  end subroutine to_units

  subroutine map_summary()
    integer :: m
    real :: a, b
    do m = 1, nmap
      if (any(m == [6, 8, 9, 10])) cycle
      if (count(.not. ieee_is_nan(mp(:,:,m))) == 0) cycle
      a = minval(mp(:,:,m), mask=.not. ieee_is_nan(mp(:,:,m)))
      b = maxval(mp(:,:,m), mask=.not. ieee_is_nan(mp(:,:,m)))
      write(*,'(4x,a8,a,2es11.3,2x,a)') mapname(m), ' range', a, b, trim(map_unit(m))
    end do
  end subroutine map_summary

  function map_unit(m) result(s)
    integer, intent(in) :: m
    character(len=12) :: s
    select case (m)
    case (1, 2)       ; s = ulabel('j')
    case (3)          ; s = ulabel('phi')
    case (4, 5, 11, 12) ; s = ulabel('v')
    case (7)          ; s = ulabel('eflux')
    case (8, 9, 10)   ; s = ulen()
    case default      ; s = ''
    end select
  end function map_unit

  !================================================================================
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

  pure integer(int64) function rec(m, k)
    integer, intent(in) :: m, k
    rec = int(m-1,int64)*nz2 + k
  end function rec

  !================================================================================
  ! full-resolution block near the mapping points (all 8 variables)
  !================================================================================
  subroutine load_local()
    integer :: ia(3), ib(3), m, k, j, nxl
    real :: rr
    real(real32), allocatable :: row(:)
    rr = r_map + max(jpar_avg, 0.0) + (jpar_delta + 3.0)*h
    ia(1) = max(1,   floor(i0 - rr/h)) ; ib(1) = min(nx2, ceiling(i0 + rr/h))
    ia(2) = max(1,   floor(j0 - rr/h)) ; ib(2) = min(ny2, ceiling(j0 + rr/h))
    if (half_domain) then
      ia(3) = 1
    else
      ia(3) = max(1, floor(k0 - rr/h))
    end if
    ib(3) = min(nz2, ceiling(k0 + rr/h))
    loc%n = ib - ia + 1
    loc%dx = h
    loc%org = [(real(ia(1))-i0)*h, (real(ia(2))-j0)*h, (real(ia(3))-k0)*h]
    if (allocated(loc%a)) deallocate(loc%a)
    allocate(loc%a(loc%n(1), loc%n(2), loc%n(3), 8))
    nxl = loc%n(1)
    allocate(row(nxl))
    do m = 1, 8
      do k = ia(3), ib(3)
        do j = ia(2), ib(2)
          call seq_read_part(fin, rec(src_var(m),k), int(j-1,int64)*nx2 + ia(1) - 1, &
                             row, int(nxl,int64))
          loc%a(:, j-ia(2)+1, k-ia(3)+1, m) = row
        end do
      end do
    end do
  end subroutine load_local

  ! local block var m (1 rho, 2-4 v, 5 p, 6-8 B) -> file variable number
  integer function src_var(m)
    integer, intent(in) :: m
    select case (m)
    case (1)    ; src_var = iv_rho
    case (2:4)  ; src_var = iv_v + m - 2
    case (5)    ; src_var = iv_p
    case default; src_var = iv_b + m - 6
    end select
  end function src_var

  !================================================================================
  ! global field for tracing (Bx,By,Bz,rho,p) read with stride
  !================================================================================
  subroutine load_global()
    integer :: ia(3), ib(3), m, k, kk, mv, s, d
    real :: c(3), lo, hi
    integer :: nn(3)
    real(real32), allocatable :: pl(:,:)
    s = max(1, stride)
    c = [i0, j0, k0] ; nn = [nx2, ny2, nz2]
    do d = 1, 3
      lo = max(box_lo(d), -1.0e8) ; hi = min(box_hi(d), 1.0e8)
      ia(d) = max(1,     floor(c(d) + lo/h))
      ib(d) = min(nn(d), ceiling(c(d) + hi/h))
    end do
    if (half_domain) ia(3) = 1
    if (any(ib - ia < 2*s)) call fatal('global box too small')
    glo%n = (ib - ia)/s + 1
    glo%dx = h*real(s)
    glo%org = [(real(ia(1))-i0)*h, (real(ia(2))-j0)*h, (real(ia(3))-k0)*h]
    write(*,'(a,3i6,a,f7.2,a)') '  global field (stride): ', glo%n, '   memory ', &
         5.0*4.0*real(glo%n(1))*real(glo%n(2))*real(glo%n(3))/1024.0**3, ' GiB'
    if (allocated(glo%a)) deallocate(glo%a)
    allocate(glo%a(glo%n(1), glo%n(2), glo%n(3), 5))
    allocate(pl(nx2, ny2))
    do m = 1, 5
      mv = merge(iv_b + m - 1, merge(iv_rho, iv_p, m == 4), m <= 3)
      do kk = 1, glo%n(3)
        k = ia(3) + (kk-1)*s
        call seq_read_rec(fin, rec(mv,k), pl, int(nx2,int64)*ny2)
        glo%a(:,:,kk,m) = pl(ia(1):ia(1)+(glo%n(1)-1)*s:s, ia(2):ia(2)+(glo%n(2)-1)*s:s)
      end do
    end do
    deallocate(pl)
  end subroutine load_global

  !================================================================================
  ! trilinear sampling. z<0 in half-domain: mirror with parity.
  ! list(:) = variable slots in the block; par(:) = parity under z -> -z
  ! returns .false. if outside the block
  !================================================================================
  logical function sample(bk, pos, list, par, val)
    type(block), intent(in) :: bk
    real, intent(in) :: pos(3)
    integer, intent(in) :: list(:), par(:)
    real, intent(out) :: val(:)
    real :: p(3), f(3), w(3)
    integer :: i(3), m
    logical :: mir
    p = pos ; mir = .false.
    if (half_domain .and. p(3) < 0.0) then
      p(3) = -p(3) ; mir = .true.
    end if
    f = (p - bk%org)/bk%dx + 1.0
    i = floor(f)
    sample = all(i >= 1) .and. all(i <= bk%n - 1)
    if (.not. sample) return
    w = f - real(i)
    do m = 1, size(list)
      val(m) = (1-w(3))*((1-w(2))*((1-w(1))*bk%a(i(1),i(2),i(3),list(m))     &
                                  +   w(1) *bk%a(i(1)+1,i(2),i(3),list(m)))  &
                        +   w(2) *((1-w(1))*bk%a(i(1),i(2)+1,i(3),list(m))   &
                                  +   w(1) *bk%a(i(1)+1,i(2)+1,i(3),list(m))))&
             +   w(3) *((1-w(2))*((1-w(1))*bk%a(i(1),i(2),i(3)+1,list(m))    &
                                  +   w(1) *bk%a(i(1)+1,i(2),i(3)+1,list(m)))&
                        +   w(2) *((1-w(1))*bk%a(i(1),i(2)+1,i(3)+1,list(m)) &
                                  +   w(1) *bk%a(i(1)+1,i(2)+1,i(3)+1,list(m))))
      if (mir) val(m) = val(m)*real(par(m))
    end do
  end function sample

  ! dipole field (simulation definition)
  pure function bdip(p) result(b)
    real, intent(in) :: p(3)
    real :: b(3), r2, r5
    r2 = sum(p**2) ; r5 = r2*r2*sqrt(r2)
    b(1) = -3.0*b0*p(1)*p(3)/r5
    b(2) = -3.0*b0*p(2)*p(3)/r5
    b(3) =  b0*(p(1)**2 + p(2)**2 - 2.0*p(3)**2)/r5
  end function bdip

  !================================================================================
  ! all maps
  !================================================================================
  subroutine compute_maps()
    integer :: i, j
    real :: q(16), p(3), bm(3), bb2, jb
    real :: th0, phj, ell, rm, thm, xj, yj, sth, cth, br_ion
    real :: vr, vt, vp, dldt, rs, v(3), er(3), et(3), ep(3)
    real :: res(5)
    logical :: ok
    mp = rnan
    !$omp parallel do schedule(dynamic,1) collapse(2) &
    !$omp private(i,j,q,p,bm,bb2,jb,th0,phj,ell,rm,thm,xj,yj,sth,cth, &
    !$omp         br_ion,vr,vt,vp,dldt,rs,v,er,et,ep,res,ok)
    do j = 1, nmlt
      do i = 1, ncolat
        th0 = real((real(i,real64)-0.5_real64)*dth)           ! angle from the pole
        phj = real(real(j-1,real64)*dph - pi)                 ! angle from noon (j=1: MLT 0)
        ell = 1.0/sin(th0)**2
        rm  = min(r_map, ell)
        mp(i,j,10) = rm
        if (rm < r_valid) cycle
        thm = asin(min(sqrt(rm/ell), 1.0))                   ! angle of the mapping point from the pole
        ! JSO: X = sin cos(phi), Y = sin sin(phi) ; sim: x=-X, y=-Y
        xj = rm*sin(thm)*cos(phj) ; yj = rm*sin(thm)*sin(phj)
        p = [-xj, -yj, real(hs)*rm*cos(thm)]
        ! ---- B, v at mapping point
        ok = sample(loc, p, [6,7,8,2,3,4], [-1,-1,1,1,1,-1], q(1:6))
        if (.not. ok) cycle
        ! ---- J_par/B = (J.B)/B^2 (J = rot(B - Bdip)) averaged along the field line
        jb = jpar_line(p, ok)
        if (.not. ok) cycle
        bm  = q(1:3)
        bb2 = sum(bm**2)
        if (bb2 <= 0.0) cycle
        ! B_r = -2 b0 cos(theta) at the ionosphere (true colatitude theta)
        cth = real(hs)*cos(th0)
        br_ion = -2.0*b0*cth
        mp(i,j,1) = jb*br_ion
        ! ---- convection mapping (spherical basis of the simulation frame)
        v = q(4:6)
        if (subtract_corot) v = v - omega*[-p(2), p(1), 0.0]
        rs  = sqrt(p(1)**2 + p(2)**2)
        er  = p/rm
        et  = [p(1)*p(3)/(rm*rs), p(2)*p(3)/(rm*rs), -rs/rm]   ! theta-hat (true colatitude)
        ep  = [-p(2)/rs, p(1)/rs, 0.0]
        vr = sum(v*er) ; vt = sum(v*et) ; vp = sum(v*ep)
        sth = rs/rm                                            ! sin(theta) at map point
        cth = p(3)/rm
        dldt = vr/sth**2 - 2.0*cth*vt/sth**3
        ! d theta0/dt (true colatitude) ; cos(theta0) = hs*cos(th0)
        ! d(theta_true)/dt = -(dL/dt)/(2 L^1.5 cos theta_true), cos theta_true = hs cos(th0)
        ! away from the pole = hs * d(theta_true)/dt
        mp(i,j,4) = -dldt/(2.0*ell**1.5*cos(th0))                        ! positive away from the pole
        mp(i,j,5) = sin(th0)*vp/(rm*sth)                                 ! positive toward increasing MLT
        ! ---- field-line tracing
        if (do_trace) then
          call trace(p, res)
          mp(i,j,6) = res(1) ; mp(i,j,7) = res(2) ; mp(i,j,8) = res(3) ; mp(i,j,9) = res(4)
        end if
      end do
    end do
    !$omp end parallel do
    call potential()
  end subroutine compute_maps

  !================================================================================
  ! (J.B)/B^2 at point p.  J = rot(B - B_dip), central differences of width
  ! jpar_delta*h (wider than the grid spacing to suppress 2-cell oscillations)
  !================================================================================
  real function jpar_point(p, ok)
    real, intent(in) :: p(3)
    logical, intent(out) :: ok
    real :: d, pp(3), pm(3), bp(3), bmm(3), bm(3)
    real(real64) :: jc(3), db(3)
    integer :: ax
    jpar_point = 0.0
    ok = sample(loc, p, [6,7,8], [-1,-1,1], bm)
    if (.not. ok) return
    d = max(jpar_delta, 0.5)*h
    jc = 0.0
    do ax = 1, 3
      pp = p ; pp(ax) = p(ax) + d
      pm = p ; pm(ax) = p(ax) - d
      ok = sample(loc, pp, [6,7,8], [-1,-1,1], bp)
      if (ok) ok = sample(loc, pm, [6,7,8], [-1,-1,1], bmm)
      if (.not. ok) return
      ! d(B-Bdip)/dx_ax (difference taken in double precision)
      db = ((real(bp,real64) - real(bmm,real64)) - (real(bdip(pp),real64) - real(bdip(pm),real64))) &
           /(2.0_real64*real(d,real64))
      select case (ax)
      case (1) ; jc(2) = jc(2) - db(3) ; jc(3) = jc(3) + db(2)
      case (2) ; jc(1) = jc(1) + db(3) ; jc(3) = jc(3) - db(1)
      case (3) ; jc(1) = jc(1) - db(2) ; jc(2) = jc(2) + db(1)
      end select
    end do
    if (sum(bm**2) <= 0.0) then
      ok = .false. ; return
    end if
    jpar_point = real(sum(jc*real(bm,real64))/sum(real(bm,real64)**2))
  end function jpar_point

  ! trace jpar_avg outward along the field line from the mapping point and average (J.B)/B^2
  real function jpar_line(p0, ok)
    real, intent(in) :: p0(3)
    logical, intent(out) :: ok
    real :: p(3), b(3), sgn, ds, s, acc, x
    integer :: n
    logical :: ok2
    jpar_line = jpar_point(p0, ok)
    if (.not. ok .or. jpar_avg <= 0.0) return
    ds = 0.5*h
    acc = jpar_line ; n = 1
    p = p0 ; s = 0.0
    if (.not. sample(loc, p, [6,7,8], [-1,-1,1], b)) return
    sgn = sign(1.0, sum(b*p))
    do while (s < jpar_avg)
      if (.not. sample(loc, p, [6,7,8], [-1,-1,1], b)) exit
      if (sum(b**2) <= 0.0) exit
      p = p + ds*sgn*b/sqrt(sum(b**2))
      s = s + ds
      x = jpar_point(p, ok2)
      if (.not. ok2) exit
      acc = acc + x ; n = n + 1
    end do
    jpar_line = acc/real(n)
  end function jpar_line

  !================================================================================
  ! field-line tracing from p (outward). res = topology, eflux, length, r_eq
  !   topology: 0 closed, 1 open, 2 undetermined
  !================================================================================
  subroutine trace(pstart, res)
    real, intent(in) :: pstart(3)
    real, intent(out) :: res(5)
    real :: p(3), pm(3), b(5), bh(3), sgn, s, ds, r, w, ws, we, zs, rmin
    integer :: n, nmax
    logical :: ok
    res = [2.0, rnan, 0.0, rnan, 0.0]
    p = pstart
    ds = trace_ds
    nmax = int(trace_maxlen/ds)
    if (.not. getb(p, b)) return
    sgn = sign(1.0, sum(b(1:3)*p))           ! B.r > 0: +B points outward
    zs  = sign(1.0, p(3))
    s = 0.0 ; ws = 0.0 ; we = 0.0 ; rmin = huge(1.0)
    do n = 1, nmax
      bh = sgn*b(1:3)/sqrt(sum(b(1:3)**2))
      pm = p + 0.5*ds*bh
      ok = getb(pm, b)
      if (.not. ok) then
        res(1) = 1.0 ; exit                     ! left the domain -> open
      end if
      if (sum(b(1:3)**2) <= 0.0) exit
      bh = sgn*b(1:3)/sqrt(sum(b(1:3)**2))
      p = p + ds*bh
      s = s + ds
      r = sqrt(sum(p**2))
      if (p(3)*zs <= 0.0) then                  ! crossed the equatorial plane
        res(4) = sqrt(p(1)**2 + p(2)**2)
        if (half_domain) then
          res(1) = 0.0 ; exit                   ! mirror symmetry -> closed
        end if
        zs = -zs
      end if
      if (.not. getb(p, b)) then
        res(1) = 1.0 ; exit
      end if
      if (r < r_map - 0.5*h .and. s > 2.0*r_map) then
        res(1) = 0.0 ; exit                     ! returned to the planet
      end if
      ! line average with the weight r^2/(r_in^2+r^2) of the legacy prodpc5 (up to the equator)
      if (ieee_is_nan(res(4))) then
        w = ds*r*r/(r_in*r_in + r*r)
        ws = ws + w
        we = we + w*sqrt(max(b(5)/max(b(4),1.0e-30), 0.0))*b(5)
      end if
    end do
    res(3) = s
    if (ws > 0.0) res(2) = we/ws
  end subroutine trace

  ! B (3), rho, p at p from the finest block containing it
  logical function getb(p, b)
    real, intent(in) :: p(3)
    real, intent(out) :: b(5)
    getb = sample(loc, p, [6,7,8,1,5], [-1,-1,1,1,1], b)
    if (.not. getb) getb = sample(glo, p, [1,2,3,4,5], [-1,-1,1,1,1], b)
    if (getb) then
      getb = all(p >= xlo) .and. all(p <= xhi)
    end if
  end function getb

  !================================================================================
  ! ionospheric electric field, its divergence and the potential
  !   coordinates: th = angle from the pole (north: colatitude, south: pi - colatitude), ph = MLT angle
  !   E = -u x B_ion  (B_ion ~ B_r r_hat)  ->  E_th = -u_ph B_r , E_ph = u_th B_r
  !   (u_th, E_th positive away from the pole; signs arranged so that the same formulae hold in the south)
  !================================================================================
  subroutine potential()
    real(real64), allocatable :: eth(:,:), eph(:,:), dv(:,:), phi(:,:), st(:), sf(:)
    real(real64) :: br, a, b, f1, f2, ef
    integer :: i, j, jp, jm
    allocate(eth(ncolat,nmlt), eph(ncolat,nmlt), dv(ncolat,nmlt), phi(ncolat,nmlt))
    allocate(st(ncolat), sf(0:ncolat))
    do i = 1, ncolat
      st(i) = sin((real(i,real64)-0.5_real64)*dth)
    end do
    do i = 0, ncolat
      sf(i) = sin(real(i,real64)*dth)
    end do
    do j = 1, nmlt
      do i = 1, ncolat
        br = -2.0_real64*b0*cos((real(i,real64)-0.5_real64)*dth)*hs   ! B_r at ionosphere
        if (ieee_is_nan(mp(i,j,4))) then
          eth(i,j) = 0.0 ; eph(i,j) = 0.0
        else
          ! the unit vector away from the pole is +theta_hat in the north and -theta_hat in the south;
          ! increasing MLT is +phi_hat in both. Components of E_t = -u x B_r r_hat:
          !   E_away = -hs * u_phi * B_r * hs ... which simplifies to
          eth(i,j) = -real(mp(i,j,5),real64)*br*hs
          eph(i,j) =  real(mp(i,j,4),real64)*br*hs
        end if
      end do
    end do
    ! div E (finite volume in th, central in ph)
    do j = 1, nmlt
      jp = mod(j, nmlt) + 1 ; jm = mod(j-2+nmlt, nmlt) + 1
      do i = 1, ncolat
        if (i < ncolat) then
          f2 = sf(i)*0.5_real64*(eth(i,j)+eth(i+1,j))
        else
          f2 = sf(i)*(1.5_real64*eth(i,j) - 0.5_real64*eth(i-1,j))
        end if
        if (i > 1) then
          f1 = sf(i-1)*0.5_real64*(eth(i,j)+eth(i-1,j))
        else
          f1 = 0.0_real64
        end if
        dv(i,j) = (f2 - f1)/(dth*st(i)) + (eph(i,jp)-eph(i,jm))/(2.0_real64*dph*st(i))
      end do
    end do
    ! j_conv = -div E (Sigma_P = 1), upward positive
    do j = 1, nmlt
      do i = 1, ncolat
        if (.not. ieee_is_nan(mp(i,j,4))) mp(i,j,2) = real(-dv(i,j))
      end do
    end do
    ! lap phi = -div E , phi = 0 at th = colat_max
    call poisson_cap(-dv, phi)
    mp(:,:,3) = real(phi)
    ! potential flow u = E x B / B^2 from E = -grad phi
    do j = 1, nmlt
      jp = mod(j, nmlt) + 1 ; jm = mod(j-2+nmlt, nmlt) + 1
      do i = 1, ncolat
        br = -2.0_real64*b0*cos((real(i,real64)-0.5_real64)*dth)*hs
        if (i == 1) then
          a = (phi(2,j) - phi(1,j))/dth
        else if (i == ncolat) then
          a = (-phi(i,j) - phi(i-1,j))/(2.0_real64*dth)         ! phi=0 at the edge
        else
          a = (phi(i+1,j) - phi(i-1,j))/(2.0_real64*dth)
        end if
        b = (phi(i,jp) - phi(i,jm))/(2.0_real64*dph*st(i))
        ! E_th = -a, E_ph = -b ; invert eth/eph relations above
        ef = -a ; mp(i,j,12) = real(-ef/(br*hs))      ! u_ph  from E_th = -u_ph br hs
        ef = -b ; mp(i,j,11) = real( ef/(br*hs))      ! u_th  from E_ph =  u_th br hs
      end do
    end do
  end subroutine potential

  !================================================================================
  ! Poisson equation on the polar cap of the unit sphere (direct solver)
  !   (1/s) d/dth (s dphi/dth) + 1/s^2 d2phi/dph2 = rhs ,  phi(th_max) = 0
  !   azimuthal real DFT + tridiagonal (Thomas) for each wavenumber
  !================================================================================
  subroutine poisson_cap(rhs, phi)
    real(real64), intent(in)  :: rhs(:,:)
    real(real64), intent(out) :: phi(:,:)
    real(real64), allocatable :: cr(:,:), ci(:,:), ar(:), br(:), cc(:), dr(:), di(:), &
                                 xr(:), xi(:), st(:), sf(:), cs(:,:), sn(:,:)
    integer :: i, j, m, mmax
    real(real64) :: ang, w
    allocate(st(ncolat), sf(0:ncolat))
    do i = 1, ncolat
      st(i) = sin((real(i,real64)-0.5_real64)*dth)
    end do
    do i = 0, ncolat
      sf(i) = sin(real(i,real64)*dth)
    end do
    mmax = nmlt/2
    allocate(cr(ncolat,0:mmax), ci(ncolat,0:mmax), cs(nmlt,0:mmax), sn(nmlt,0:mmax))
    do m = 0, mmax
      do j = 1, nmlt
        ang = real(m,real64)*real(j-1,real64)*dph
        cs(j,m) = cos(ang) ; sn(j,m) = sin(ang)
      end do
    end do
    ! forward DFT
    do m = 0, mmax
      do i = 1, ncolat
        cr(i,m) = sum(rhs(i,:)*cs(:,m))/real(nmlt,real64)
        ci(i,m) = sum(rhs(i,:)*sn(:,m))/real(nmlt,real64)
      end do
    end do
    allocate(ar(ncolat), br(ncolat), cc(ncolat), dr(ncolat), di(ncolat), xr(ncolat), xi(ncolat))
    do m = 0, mmax
      do i = 1, ncolat
        ar(i) = sf(i-1)/(dth*dth*st(i))
        cc(i) = sf(i)  /(dth*dth*st(i))
        br(i) = -(ar(i) + cc(i)) - real(m*m,real64)/st(i)**2
      end do
      ar(1) = 0.0_real64                         ! pole: sf(0) = 0
      br(ncolat) = br(ncolat) - cc(ncolat)       ! Dirichlet 0 at the edge face (ghost = -phi)
      cc(ncolat) = 0.0_real64
      dr = cr(:,m) ; di = ci(:,m)
      call thomas(ar, br, cc, dr, xr)
      call thomas(ar, br, cc, di, xi)
      cr(:,m) = xr ; ci(:,m) = xi
    end do
    ! inverse
    do j = 1, nmlt
      do i = 1, ncolat
        w = cr(i,0)
        do m = 1, mmax
          if (2*m == nmlt) then
            w = w + cr(i,m)*cs(j,m)
          else
            w = w + 2.0_real64*(cr(i,m)*cs(j,m) + ci(i,m)*sn(j,m))
          end if
        end do
        phi(i,j) = w
      end do
    end do
  end subroutine poisson_cap

  subroutine thomas(a, b, c, d, x)
    real(real64), intent(in) :: a(:), b(:), c(:), d(:)
    real(real64), intent(out) :: x(:)
    real(real64), allocatable :: cp(:), dp(:)
    integer :: n, i
    n = size(b)
    allocate(cp(n), dp(n))
    cp(1) = c(1)/b(1) ; dp(1) = d(1)/b(1)
    do i = 2, n
      cp(i) = c(i)/(b(i) - a(i)*cp(i-1))
      dp(i) = (d(i) - a(i)*dp(i-1))/(b(i) - a(i)*cp(i-1))
    end do
    x(n) = dp(n)
    do i = n-1, 1, -1
      x(i) = dp(i) - cp(i)*x(i+1)
    end do
  end subroutine thomas

  !================================================================================
  subroutine write_maps(base)
    character(len=*), intent(in) :: base
    integer :: uu, m
    open(newunit=uu, file=base//'.bin', access='stream', form='unformatted', status='replace')
    write(uu) mp
    close(uu)
    open(newunit=uu, file=base//'.txt', status='replace')
    write(uu,'(a)') '# polar_map output: float32 native, Fortran order map(ncolat,nmlt,nmap)'
    write(uu,'(a)') 'file      = '//fni
    write(uu,'(a)') 'hemisphere= '//trim(hemisphere)
    write(uu,'(a,i0,1x,i0,1x,i0)') 'ncolat nmlt nmap = ', ncolat, nmlt, nmap
    write(uu,'(a,f10.4,a)') 'colat(i) = (i-0.5)*', colat_max/real(ncolat), ' deg from the pole'
    write(uu,'(a,f10.4,a)') 'mlt(j)   = (j-1)*', 24.0/real(nmlt), ' h'
    write(uu,'(a,f8.3,a,f8.3,a,l1)') 'r_map = ', r_map, '  r_valid = ', r_valid, &
                                     '  subtract_corot = ', subtract_corot
    if (UC%phys) then
      write(uu,'(a,a)') 'units     = physical, planet = ', trim(UC%planet)
    else
      write(uu,'(a)') 'units     = normalized (simulation units)'
    end if
    do m = 1, nmap
      write(uu,'(i3,2x,a8,2x,a)') m, mapname(m), trim(map_unit(m))
    end do
    write(uu,'(a)') 'j_mhd : (J.B/B^2)_map * B_r,ion  (upward +, J = rot B / mu0, x jfac)'
    if (UC%phys) then
      write(uu,'(a,g0,a)') 'j_conv: -Sigma_P div(E), Sigma_P = ', sigma_p, ' S (upward +)'
    else
      write(uu,'(a)') 'j_conv: -div(E) with Sigma_P = 1 (upward +)'
    end if
    write(uu,'(a)') 'phi   : potential, lap phi = -div E, phi = 0 at colat_max'
    write(uu,'(a)') 'u_colat/u_mlt: mapped ionospheric flow (speed at R = 1 planetary radius),'// &
                    ' + = away from pole / + = increasing MLT'
    write(uu,'(a)') 'topology: 0 closed, 1 open, 2 undetermined'
    write(uu,'(a)') 'eflux : <p sqrt(p/rho)> along the line (weight r^2/(r_in^2+r^2), up to the equator)'
    write(uu,'(a)') 'up_colat/up_mlt: E x B flow from -grad(phi)'
    close(uu)
  end subroutine write_maps

  !================================================================================
  ! plotting
  !================================================================================
  ! polar grid -> square image (n x n), radius = angle from the pole
  subroutine to_image(fld, n, img)
    real(real32), intent(in) :: fld(:,:)
    integer, intent(in) :: n
    real(real32), intent(out) :: img(n,n)
    integer :: ix, iy, i1, j1, j2
    real :: xd, yd, rho, psi, rmax, fi, fj, wi, wj, sgx, v11, v12, v21, v22
    rmax = colat_max*real(pi)/180.0
    sgx = merge(-1.0, 1.0, view == 'above')
    do iy = 1, n
      do ix = 1, n
        xd = rmax*(2.0*real(ix-1)/real(n-1) - 1.0)
        yd = rmax*(2.0*real(iy-1)/real(n-1) - 1.0)
        rho = sqrt(xd*xd + yd*yd)
        if (rho > rmax) then
          img(ix,iy) = rnan ; cycle
        end if
        psi = atan2(sgx*xd, yd)
        if (psi < 0.0) psi = psi + 2.0*real(pi)
        fi = rho/real(dth) + 0.5
        fj = (psi + real(pi))/real(dph) + 1.0
        i1 = floor(fi) ; wi = fi - real(i1)
        if (i1 < 1) then
          i1 = 1 ; wi = 0.0
        else if (i1 >= ncolat) then
          i1 = ncolat - 1 ; wi = 1.0
        end if
        j1 = floor(fj) ; wj = fj - real(j1)
        j1 = mod(j1-1+nmlt, nmlt) + 1 ; j2 = mod(j1, nmlt) + 1
        v11 = fld(i1,j1) ; v21 = fld(i1+1,j1) ; v12 = fld(i1,j2) ; v22 = fld(i1+1,j2)
        img(ix,iy) = (1-wi)*((1-wj)*v11 + wj*v12) + wi*((1-wj)*v21 + wj*v22)
      end do
    end do
  end subroutine to_image

  subroutine draw_page(ino, ttl)
    integer, intent(in) :: ino
    character(len=*), intent(in) :: ttl
    real, parameter :: px0 = 2.0, py0 = 2.0, sz = 7.0, gap = 2.8
    real :: xb(4), yb(4)
    character(len=64) :: s
    xb = [px0, px0+sz+gap, px0, px0+sz+gap]
    yb = [py0+sz+gap, py0+sz+gap, py0, py0]
    call ps_page(2, 1, 1, ino)                 ! A4 portrait (as the legacy plot)
    call ps_pen(2)
    call ps_text(px0, py0+2*sz+gap+2.9, 0.45, trim(header))
    write(s,'(a)') merge('Northern Hemisphere of Ionosphere', 'Southern Hemisphere of Ionosphere', &
                         hemisphere == 'north')
    call ps_text(px0, py0+2*sz+gap+2.2, 0.45, trim(s))
    call ps_text(px0, py0+2*sz+gap+1.5, 0.42, trim(ttl))
    ! 1: FAC (MHD) + OCB
    call panel_color(1, xb(1), yb(1), sz, 'Field-Aligned Current (MHD, up +)'//ubr(1), 8, jmhd_max, .true.)
    ! 2: potential + flow
    call panel_potential(xb(2), yb(2), sz)
    ! 3: FAC from convection (uniform conductance)
    call panel_color(2, xb(3), yb(3), sz, 'FAC from convection (div E, up +)'//ubr(2), 8, jconv_max, .true.)
    ! 4: energy flux proxy + OCB
    call panel_color(7, xb(4), yb(4), sz, 'Energy flux proxy <p sqrt(T)>'//ubr(7), 6, eflux_max, .false.)
    call ps_endpage()
  end subroutine draw_page

  function ubr(m) result(s)
    integer, intent(in) :: m
    character(len=:), allocatable :: s
    s = ''
    if (len_trim(map_unit(m)) > 0) s = ' ['//trim(map_unit(m))//']'
  end function ubr

  subroutine frame(xb, yb, sz)
    real, intent(in) :: xb, yb, sz
    real :: cx, cy, rr, lat, sgx
    integer :: k
    character(len=8) :: s
    cx = xb + 0.5*sz ; cy = yb + 0.5*sz
    sgx = merge(-1.0, 1.0, view == 'above')
    call ps_pen(1) ; call ps_rgb(0.0, 0.0, 0.0)
    do k = 1, 8
      lat = 90.0 - 5.0*real(k)
      if (90.0 - lat > colat_max + 0.01) exit
      rr = 0.5*sz*(90.0-lat)/colat_max
      call ps_circle(cx, cy, rr)
      write(s,'(i0)') nint(lat)
      call ps_text(cx + rr*0.707 + 0.05, cy - rr*0.707 - 0.3, 0.22, trim(s)//'\312')
    end do
    call ps_line(xb, cy, xb+sz, cy) ; call ps_line(cx, yb, cx, yb+sz) ; call ps_stroke()
    call ps_text(cx, yb+sz+0.15, 0.3, '12', align='c')
    call ps_text(cx, yb-0.5, 0.3, '00', align='c')
    call ps_text(xb-0.15, cy-0.12, 0.3, merge('18', '06', sgx < 0.0), align='r')
    call ps_text(xb+sz+0.15, cy-0.12, 0.3, merge('06', '18', sgx < 0.0))
  end subroutine frame

  subroutine ocb(xb, yb, sz, n)
    real, intent(in) :: xb, yb, sz
    integer, intent(in) :: n
    real(real32), allocatable :: img(:,:), t(:,:)
    if (.not. do_trace) return
    allocate(img(n,n), t(ncolat,nmlt))
    t = mp(:,:,6)
    where (ieee_is_nan(t) .or. t > 1.5) t = rnan          ! undetermined -> skip
    call to_image(t, n, img)
    call ps_pen(3) ; call ps_rgb(0.0, 0.0, 0.0)
    call ps_contour(n, n, img, xb, yb, sz, sz, 0.5)
    call ps_pen(1)
  end subroutine ocb

  !--------------------------------------------------------------------------------
  ! smoothing for the plots (5-point smoothing with coefficient 1/8 as legacy diff2, nsmooth times).
  ! periodic in MLT; the poleward and equatorward edges use the point itself. NaN is not smoothed.
  !--------------------------------------------------------------------------------
  function smooth_map(a) result(b)
    real(real32), intent(in) :: a(:,:)
    real(real32) :: b(size(a,1), size(a,2)), c(size(a,1), size(a,2)), nb(4)
    integer :: it, i, j, ni, nj, k
    real :: w
    ni = size(a,1) ; nj = size(a,2)
    b = a
    do it = 1, nsmooth
      c = b
      do j = 1, nj
        do i = 1, ni
          if (ieee_is_nan(c(i,j))) cycle
          nb = [c(max(i-1,1),j), c(min(i+1,ni),j), c(i,mod(j-2+nj,nj)+1), c(i,mod(j,nj)+1)]
          w = 0.5 ; b(i,j) = 0.5*c(i,j)
          do k = 1, 4
            if (.not. ieee_is_nan(nb(k))) then
              b(i,j) = b(i,j) + 0.125*nb(k) ; w = w + 0.125
            end if
          end do
          b(i,j) = b(i,j)/w
        end do
      end do
    end do
  end function smooth_map

  subroutine panel_color(m, xb, yb, sz, lab, ico, vmx, symm)
    integer, intent(in) :: m, ico
    real, intent(in) :: xb, yb, sz, vmx
    character(len=*), intent(in) :: lab
    logical, intent(in) :: symm
    integer, parameter :: n = 401
    real(real32), allocatable :: img(:,:)
    real :: lo, hi
    character(len=32) :: s1, s2
    allocate(img(n,n))
    call to_image(smooth_map(mp(:,:,m)), n, img)
    hi = vmx
    if (hi <= 0.0) hi = maxval(abs(img), mask=(.not. ieee_is_nan(img)))
    if (hi <= 0.0) hi = 1.0
    lo = merge(-hi, 0.0, symm)
    call ps_image(n, n, xb, yb, sz, sz, img, lo, hi, ico, 64, 8)
    call ocb(xb, yb, sz, n)
    call frame(xb, yb, sz)
    call ps_text(xb, yb+sz+0.6, 0.3, lab)
    call ps_colorbar(xb+sz*0.25, yb-1.2, sz*0.5, 0.2, ico)
    write(s1,'(es10.2)') lo ; write(s2,'(es10.2)') hi
    call ps_text(xb+sz*0.25-0.1, yb-1.18, 0.22, trim(adjustl(s1)), align='r')
    call ps_text(xb+sz*0.75+0.1, yb-1.18, 0.22, trim(adjustl(s2)))
  end subroutine panel_color

  subroutine panel_potential(xb, yb, sz)
    real, intent(in) :: xb, yb, sz
    integer, parameter :: n = 401
    real(real32), allocatable :: img(:,:)
    real :: pmx, lev, cx, cy, rmax, sgx, uv, ua, um, ur, psi, rho, x1, y1, x2, y2, sc, emx
    integer :: k, i, j, ii, jj
    character(len=96) :: s
    allocate(img(n,n))
    call to_image(mp(:,:,3), n, img)
    pmx = maxval(abs(img), mask=(.not. ieee_is_nan(img)))
    if (pmx <= 0.0) pmx = 1.0
    ! no contours if the potential is negligible compared with |E| x domain width
    emx = maxval(sqrt(mp(:,:,4)**2 + mp(:,:,5)**2), mask=.not. ieee_is_nan(mp(:,:,4))) &
          * 2.0*abs(b0)*colat_max*real(pi)/180.0 * (fphi/fu)
    if (pmx < 1.0e-4*emx) pmx = -1.0
    if (pmx > 0.0) then
    do k = -nphi_cont, nphi_cont
      if (k == 0) cycle
      lev = pmx*real(k)/real(nphi_cont+1)
      call ps_pen(2)
      if (k > 0) then
        call ps_rgb(1.0, 0.0, 0.0)
      else
        call ps_rgb(0.0, 0.0, 1.0)
      end if
      call ps_contour(n, n, img, xb, yb, sz, sz, lev)
    end do
    end if
    call ps_rgb(0.0, 0.0, 0.0)
    ! flow arrows (mapped flow)
    cx = xb + 0.5*sz ; cy = yb + 0.5*sz
    rmax = colat_max*real(pi)/180.0
    sgx = merge(-1.0, 1.0, view == 'above')
    um = umax
    if (um <= 0.0) um = maxval(sqrt(mp(:,:,4)**2 + mp(:,:,5)**2), &
                               mask=(.not. ieee_is_nan(mp(:,:,4))))
    if (um <= 0.0) um = 1.0
    sc = 0.9*(0.5*sz/real(vec_ncolat))/um               ! maximum speed = 0.9 x grid spacing
    call ps_pen(1) ; call ps_rgb(0.0, 0.5, 0.0)
    do ii = 1, vec_ncolat
      i = nint((real(ii)-0.5)*real(ncolat)/real(vec_ncolat) + 0.5)
      i = max(1, min(ncolat, i))
      do jj = 1, vec_nmlt
        j = nint(real(jj-1)*real(nmlt)/real(vec_nmlt)) + 1
        j = max(1, min(nmlt, j))
        ur = mp(i,j,4) ; ua = mp(i,j,5)
        if (ieee_is_nan(ur)) cycle
        rho = real((real(i)-0.5)*dth)
        psi = real((j-1)*dph - pi)
        x1 = cx + 0.5*sz*sgx*rho*sin(psi)/rmax
        y1 = cy + 0.5*sz*rho*cos(psi)/rmax
        ! e_rho = (sgx sin, cos), e_psi = (sgx cos, -sin)
        uv = sqrt(ur*ur + ua*ua)
        if (uv > um) then
          ur = ur*um/uv ; ua = ua*um/uv
        end if
        x2 = x1 + sc*(ur*sgx*sin(psi) + ua*sgx*cos(psi))
        y2 = y1 + sc*(ur*cos(psi) - ua*sin(psi))
        if (abs(x2-x1) + abs(y2-y1) < 0.02) cycle
        call ps_arrowhead(x1, y1, x2, y2, 0.12*min(1.0,uv/um)+0.03, 0.4)
      end do
    end do
    call ps_stroke()
    call ocb(xb, yb, sz, n)
    call frame(xb, yb, sz)
    call ps_text(xb, yb+sz+0.6, 0.3, 'Potential & mapped flow')
    if (pmx > 0.0) then
      write(s,'(a,es9.2,1x,a,a,es9.2,1x,a)') 'phi max ', pmx, trim(map_unit(3)), ',  |u| max ', um, &
                                             trim(map_unit(4))
    else
      write(s,'(a,es9.2,1x,a)') 'phi ~ 0,  |u| max ', um, trim(map_unit(4))
    end if
    call ps_text(xb+sz*0.5, yb-1.18, 0.22, trim(s), align='c')
  end subroutine panel_potential

end program polar_map
