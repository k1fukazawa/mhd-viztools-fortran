!==============================================================================
! vorticity : components of the vorticity parallel and perpendicular to B
!
!   Replaces the legacy vorticity code. Same input/output formats:
!     input   ju02nt###m.data : sequential unformatted, nb variables x nz2 planes
!     output  ju02nt###.vort  : same format, variable 1 = omega_par, 2 = omega_perp
!       omega      = rot v   (2nd-order central differences)
!       omega_par  = omega . B / |B|
!       omega_perp = |omega x B| / |B| = sqrt(|omega|^2 - omega_par^2)
!
! Changes from the legacy code
!   (1) legacy: f(nx2,ny2,nz2,8)+pp(...,3) were all kept in memory (about 260 GB)
!       new   : rolling buffer of only 3 z planes (about 0.35 GB on the standard grid)
!   (2) only the records of the variables needed (vx,vy,vz,Bx,By,Bz) are read by stream access
!       (the transposing copy of the legacy code is gone)
!   (3) the computation is OpenMP parallel
!   (4) NaN at |B|=0 or omega=0 fixed
!       (division by zero and square roots of negative numbers are avoided)
!   (5) boundary planes (i,j,k = 1, n2) are output as 0
!       (the legacy code did not compute them and wrote the values of rho and vx there)
!   (6) additional outputs wx, wy, wz, |w| selectable in the namelist
!   (7) computed in double precision (input/output in single). omega_perp = |w| sqrt(1 - cos^2)
!       lost digits in single precision when w and B were nearly parallel
!   The output stays in normalized units (1/ts). Conversion to [1/s] is done by slice_plot /
!   bin2hdf when they read it (&units, varunit='omega')
!
! Usage:  ./vorticity [namelist file]     (default: vort.nml)
!==============================================================================
program vorticity
  use mhd_common
  implicit none

  !---- namelist --------------------------------------------------------------
  integer :: nx2 = 3002, ny2 = 2002, nz2 = 1002, nb = 8
  real    :: hx = 0.15, hy = 0.15, hz = 0.15
  integer :: iv(3) = [2, 3, 4]          ! variable numbers of vx, vy, vz
  integer :: ibf(3) = [6, 7, 8]         ! variable numbers of Bx, By, Bz
  character(len=512) :: diri = './', diro = './'
  character(len=64)  :: prefix = 'ju02nt', prefix_out = ''
  character(len=64)  :: suffix_in = 'm.data', suffix_out = '.vort'
  integer :: ntap = 186, nfiles = 1, nstep = 1, ndigit = 3
  ! quantities to output (they become variables 1,2,... of the output file in this order)
  logical :: out_par = .true., out_perp = .true.
  logical :: out_wx = .false., out_wy = .false., out_wz = .false., out_wmag = .false.
  namelist /vort/ nx2, ny2, nz2, nb, hx, hy, hz, iv, ibf, diri, diro, prefix, &
                  prefix_out, suffix_in, suffix_out, ntap, nfiles, nstep,     &
                  ndigit, out_par, out_perp, out_wx, out_wy, out_wz, out_wmag

  real(real32), allocatable :: v(:,:,:,:)     ! (nx2,ny2,3 comp,3 slots)
  real(real32), allocatable :: b(:,:,:)       ! (nx2,ny2,3)
  real(real32), allocatable :: w(:,:,:)       ! (nx2,ny2,nout)
  type(seqfile) :: fin, fout
  integer :: nout, ifile, ntap1, k, c, u, ios, sm, s0, sp, o
  integer(int64) :: nplane, recb
  real(real64) :: t0, tr, tc, tw, ta
  character(len=:), allocatable :: fni, fno, nmlfile
  character(len=512) :: arg
  character(len=256) :: msg
  character(len=8) :: names(6)
  logical :: sel(6)

  nmlfile = 'vort.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=u, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(u, nml=vort, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &vort : '//trim(msg))
  close(u)
  if (len_trim(prefix_out) == 0) prefix_out = prefix

  names = ['par     ', 'perp    ', 'wx      ', 'wy      ', 'wz      ', 'wmag    ']
  sel   = [out_par, out_perp, out_wx, out_wy, out_wz, out_wmag]
  nout  = count(sel)
  if (nout == 0) call fatal('no output selected')

  nplane = int(nx2,int64)*ny2
  recb   = 4_int64*nplane
  write(*,'(a)') '----- vorticity -----'
  write(*,'(a,3i6)') ' grid (nx2,ny2,nz2): ', nx2, ny2, nz2
  write(*,'(a)') ' outputs:'
  do o = 1, 6
    if (sel(o)) write(*,'(a,a)') '   - ', trim(names(o))
  end do
  write(*,'(a,f8.2,a)') ' work memory: ', 4.0d0*nplane*(9+3+nout)/1024d0**3, ' GiB'

  allocate(v(nx2,ny2,3,3), b(nx2,ny2,3), w(nx2,ny2,nout))

  if (nfiles < 1) call fatal('no input file: nfiles must be >= 1')
  do ifile = 1, nfiles
    ntap1 = ntap + (ifile-1)*nstep
    fni = trim(make_fname(diri, prefix, ntap1, ndigit, suffix_in))
    fno = trim(make_fname(diro, prefix_out, ntap1, ndigit, suffix_out))
    write(*,'(a)') ' input : '//fni
    write(*,'(a)') ' output: '//fno
    ta = wtime() ; tr = 0 ; tc = 0 ; tw = 0

    call seq_open_read(fin, fni, recb)
    call seq_check_size(fin, int(nb,int64)*nz2)
    ! the output uses the same record markers and endianness as the input
    call seq_open_write(fout, fno, recb, fin%marker, fin%swap)

    ! k = 1 and k = nz2 : boundary -> 0
    w = 0.0
    do o = 1, nout
      call seq_write_rec(fout, rec(o,1),   w(:,:,o), nplane)
      call seq_write_rec(fout, rec(o,nz2), w(:,:,o), nplane)
    end do

    ! preload v at k = 1, 2
    t0 = wtime()
    do c = 1, 3
      call seq_read_rec(fin, rec(iv(c),1), v(:,:,c,slot(1)), nplane)
      call seq_read_rec(fin, rec(iv(c),2), v(:,:,c,slot(2)), nplane)
    end do
    tr = tr + (wtime()-t0)

    do k = 2, nz2-1
      t0 = wtime()
      sm = slot(k-1) ; s0 = slot(k) ; sp = slot(k+1)
      do c = 1, 3
        call seq_read_rec(fin, rec(iv(c),k+1), v(:,:,c,sp), nplane)
        call seq_read_rec(fin, rec(ibf(c),k),  b(:,:,c),    nplane)
      end do
      tr = tr + (wtime()-t0)

      t0 = wtime()
      call calc_plane(v(:,:,:,sm), v(:,:,:,s0), v(:,:,:,sp), b, w)
      tc = tc + (wtime()-t0)

      t0 = wtime()
      do o = 1, nout
        call seq_write_rec(fout, rec(o,k), w(:,:,o), nplane)
      end do
      tw = tw + (wtime()-t0)
    end do

    call seq_close(fin)
    call seq_close(fout)
    call report_rate('  read ', 6_int64*recb*(nz2-2), tr)
    call report_rate('  calc ', 0_int64, tc)
    call report_rate('  write', int(nout,int64)*recb*nz2, tw)
    call report_rate('  total', (6_int64+nout)*recb*nz2, wtime()-ta)
  end do

  deallocate(v, b, w)

contains

  pure integer function slot(k)
    integer, intent(in) :: k
    slot = mod(k-1, 3) + 1
  end function slot

  pure integer(int64) function rec(m, k)
    integer, intent(in) :: m, k
    rec = int(m-1,int64)*nz2 + k
  end function rec

  !----------------------------------------------------------------------------
  ! one z-plane: vm, v0, vp = v at k-1, k, k+1 ; bb = B at k
  !----------------------------------------------------------------------------
  subroutine calc_plane(vm, v0, vp, bb, ww)
    real(real32), intent(in)  :: vm(:,:,:), v0(:,:,:), vp(:,:,:), bb(:,:,:)
    real(real32), intent(out) :: ww(:,:,:)
    real(real64) :: xh, yh, zh, w1, w2, w3, bmag, omag2, dotp, par, perp, b1, b2, b3
    integer :: i, j, o
    xh = 1.0_real64/real(hx,real64)/2.0_real64
    yh = 1.0_real64/real(hy,real64)/2.0_real64
    zh = 1.0_real64/real(hz,real64)/2.0_real64
    !$omp parallel do private(i,w1,w2,w3,bmag,omag2,dotp,par,perp,o,b1,b2,b3) schedule(static)
    do j = 1, ny2
      if (j == 1 .or. j == ny2) then
        ww(:,j,:) = 0.0
        cycle
      end if
      ww(1,j,:)   = 0.0
      ww(nx2,j,:) = 0.0
      do i = 2, nx2-1
        ! (same order of operations as the legacy code)
        w1 = (real(v0(i,j+1,3),real64)-v0(i,j-1,3))*yh - (real(vp(i,j,2),real64)-vm(i,j,2))*zh
        w2 = (real(vp(i,j,1),real64)-vm(i,j,1))*zh     - (real(v0(i+1,j,3),real64)-v0(i-1,j,3))*xh
        w3 = (real(v0(i+1,j,2),real64)-v0(i-1,j,2))*xh - (real(v0(i,j+1,1),real64)-v0(i,j-1,1))*yh
        b1 = bb(i,j,1) ; b2 = bb(i,j,2) ; b3 = bb(i,j,3)
        bmag  = sqrt(b1*b1 + b2*b2 + b3*b3)
        omag2 = w1*w1 + w2*w2 + w3*w3
        dotp  = b1*w1 + b2*w2 + b3*w3
        if (bmag > 0.0_real64) then
          par = dotp/bmag
        else
          par = 0.0_real64
        end if
        ! omega_perp = |omega x B| / |B| (same as omag*sinth of the legacy code, without loss of digits)
        if (bmag > 0.0_real64) then
          perp = sqrt((w2*b3-w3*b2)**2 + (w3*b1-w1*b3)**2 + (w1*b2-w2*b1)**2)/bmag
        else
          perp = sqrt(omag2)
        end if
        o = 0
        if (out_par)  then ; o = o + 1 ; ww(i,j,o) = real(par)          ; end if
        if (out_perp) then ; o = o + 1 ; ww(i,j,o) = real(perp)         ; end if
        if (out_wx)   then ; o = o + 1 ; ww(i,j,o) = real(w1)           ; end if
        if (out_wy)   then ; o = o + 1 ; ww(i,j,o) = real(w2)           ; end if
        if (out_wz)   then ; o = o + 1 ; ww(i,j,o) = real(w3)           ; end if
        if (out_wmag) then ; o = o + 1 ; ww(i,j,o) = real(sqrt(omag2))  ; end if
      end do
    end do
    !$omp end parallel do
  end subroutine calc_plane

end program vorticity
