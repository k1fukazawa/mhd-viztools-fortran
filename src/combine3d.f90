!==============================================================================
! combine3d : combine per-process (per-rank) output files into one 3-D data file
!
!   Replaces the legacy data_conbine3d_sxt_mem.f90 (K. FUKAZAWA, 2011).
!   The output file (e.g. ju03nt000m.data) has exactly the same format as before:
!     sequential unformatted, do m=1,nb ; do k=1,onz ; write(u) f(1:onx,1:ony,k)
!
! Changes from the legacy code
!   (1) only ranks that overlap the output region are read (2048 -> 450 ranks by default)
!       legacy: the whole 6000x4000x2000 grid (192 GB) was expanded in memory, then cut out
!       new   : processed in slabs of one rank layer in z (onx*ony*nzz, about 3 GB)
!   (2) rank files are read by stream access, only the needed part of the needed records
!   (3) OpenMP: ranks in a slab are read in parallel with tasks, and writing the
!       previous slab overlaps reading the next one (overlap=.true.)
!   (4) parameters are given in a namelist (default combine.nml); no recompilation
!
! Usage:  ./combine3d [namelist file]
!==============================================================================
program combine3d
  use mhd_common
  !$ use omp_lib
  implicit none

  !---- namelist parameters (defaults = settings of the legacy code) ----------
  integer :: npex = 8, npey = 16, npez = 16       ! number of ranks in x, y, z
  integer :: nxx = 750, nyy = 250, nzz = 125      ! grid points per rank (no ghost cells)
  integer :: nb  = 8                              ! number of variables (records)
  ! output sub-region: start index and size in the global grid (1..nxx*npex etc.)
  !   legacy code: ic = i + (nx-mx)/2 - (nxp-mxp) = i + 999 -> ox0 = 1000
  !             jc = j + (ny-my)/2           = j + 999 -> oy0 = 1000
  !             kc = k                                 -> oz0 = 1
  integer :: ox0 = 1000, oy0 = 1000, oz0 = 1
  integer :: onx = 3002, ony = 2002, onz = 1002
  character(len=512) :: diri = './', diro = './'
  character(len=64)  :: prefix_in  = 'ju03nt'
  character(len=64)  :: prefix_out = ''           ! empty: same as prefix_in
  character(len=64)  :: suffix_in  = '.tdata.'    ! text before the rank number
  character(len=64)  :: suffix_out = 'm.data'
  integer :: ntap = 0, nfiles = 1, nstep = 1      ! file numbers ntap, ntap+nstep, ...
  integer :: ndigit = 3, rank_digits = 4
  logical :: overlap = .true.                     ! overlap writing and reading
  namelist /combine/ npex, npey, npez, nxx, nyy, nzz, nb, ox0, oy0, oz0,  &
                     onx, ony, onz, diri, diro, prefix_in, prefix_out,     &
                     suffix_in, suffix_out, ntap, nfiles, nstep, ndigit,   &
                     rank_digits, overlap

  !---- work -----------------------------------------------------------------
  real(real32), allocatable :: slab(:,:,:,:)      ! (onx, ony, nzz, nbuf)
  integer :: rx0, rx1, ry0, ry1, rz0, rz1
  integer :: ifile, ntap1, m, kz, rx, ry, gk0, gk1, nk, ib, pb, pnk, nbuf
  integer :: uo, ios, nthreads
  logical :: have_prev
  integer(int64) :: rank_words, bytes_read
  real(real64) :: t0, t1, tf0
  character(len=:), allocatable :: fout, nmlfile
  character(len=512) :: arg
  character(len=256) :: msg

  !---- read namelist --------------------------------------------------------
  nmlfile = 'combine.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=uo, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(uo, nml=combine, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &combine : '//trim(msg))
  close(uo)
  if (len_trim(prefix_out) == 0) prefix_out = prefix_in

  if (ox0 < 1 .or. oy0 < 1 .or. oz0 < 1 .or.                  &
      ox0+onx-1 > nxx*npex .or. oy0+ony-1 > nyy*npey .or.     &
      oz0+onz-1 > nzz*npez) call fatal('output region exceeds the global grid')

  rx0 = (ox0-1)/nxx ; rx1 = (ox0+onx-2)/nxx
  ry0 = (oy0-1)/nyy ; ry1 = (oy0+ony-2)/nyy
  rz0 = (oz0-1)/nzz ; rz1 = (oz0+onz-2)/nzz
  rank_words = int(nxx,int64)*nyy*nzz
  nbuf = merge(2, 1, overlap)
  nthreads = 1
  !$ nthreads = omp_get_max_threads()

  write(*,'(a)')        '----- combine3d -----'
  write(*,'(a,3i6)')    ' ranks (npex,npey,npez) : ', npex, npey, npez
  write(*,'(a,3i6)')    ' rank grid (nxx,nyy,nzz): ', nxx, nyy, nzz
  write(*,'(a,3i6)')    ' output origin (global) : ', ox0, oy0, oz0
  write(*,'(a,3i6)')    ' output size            : ', onx, ony, onz
  write(*,'(a,i0,a,i0)') ' ranks to read per var  : ', &
       (rx1-rx0+1)*(ry1-ry0+1)*(rz1-rz0+1), ' / ', npex*npey*npez
  write(*,'(a,f8.2,a)') ' slab buffer            : ', &
       4.0d0*onx*ony*nzz*nbuf/1024d0**3, ' GiB'
  write(*,'(a,i0,a,l1)') ' OpenMP threads         : ', nthreads, '   overlap=', overlap

  allocate(slab(onx, ony, nzz, nbuf))

  !---- loop over output files ----------------------------------------------
  if (nfiles < 1) call fatal('no input file: nfiles must be >= 1')
  do ifile = 1, nfiles
    ntap1 = ntap + (ifile-1)*nstep
    fout  = trim(make_fname(diro, prefix_out, ntap1, ndigit, suffix_out))
    write(*,'(a)') ' output: '//fout
    open(newunit=uo, file=fout, access='sequential', form='unformatted', &
         status='replace', action='write', iostat=ios, iomsg=msg)
    if (ios /= 0) call fatal('cannot create '//fout//' : '//trim(msg))
    tf0 = wtime()

    do m = 1, nb
      t0 = wtime()
      bytes_read = 0
      have_prev  = .false.
      ib = 1 ; pb = 1 ; pnk = 0
      !$omp parallel default(shared)
      !$omp single
      do kz = rz0, rz1
        gk0 = max(oz0, kz*nzz+1)
        gk1 = min(oz0+onz-1, (kz+1)*nzz)
        nk  = gk1 - gk0 + 1
        if (overlap) then
          if (have_prev) ib = 3 - pb
          if (have_prev) then
            !$omp task default(shared) firstprivate(pb, pnk)
            call write_slab(pb, pnk)
            !$omp end task
          end if
        else
          if (have_prev) call write_slab(pb, pnk)
          ib = 1
        end if
        do ry = ry0, ry1
          do rx = rx0, rx1
            !$omp task default(shared) firstprivate(rx, ry, kz, gk0, nk, ib)
            call fill_block(rx, ry, kz, gk0, nk, ib)
            !$omp end task
          end do
        end do
        !$omp taskwait
        pb = ib ; pnk = nk ; have_prev = .true.
      end do
      if (have_prev) call write_slab(pb, pnk)
      !$omp end single
      !$omp end parallel
      t1 = wtime()
      write(msg,'(a,i0,a)') '  m=', m, ' read+write'
      call report_rate(trim(msg), bytes_read + 4_int64*onx*ony*onz, t1-t0)
    end do
    close(uo)
    call report_rate('  file total (read+write)', &
         int(nb,int64)*4_int64*onx*ony*onz, wtime()-tf0)
  end do

  deallocate(slab)
  write(*,'(a)') ' End data conversion'

contains

  !----------------------------------------------------------------------------
  ! read the part global z = gk0..gk0+nk-1 of variable m of rank (rx,ry,kz)
  ! and put it at the right place of slab(:,:,1:nk,ib)
  !----------------------------------------------------------------------------
  recursive subroutine fill_block(rx, ry, kz, gk0, nk, ib)
    integer, intent(in) :: rx, ry, kz, gk0, nk, ib
    type(seqfile) :: sf
    real(real32), allocatable :: buf(:,:,:)
    integer :: irank, gi0, gi1, gj0, gj1, lk0, li0, li1, lj0, lj1, nj, k
    integer(int64) :: off, nw

    irank = rx + npex*(ry + npey*kz)
    call seq_open_read(sf, make_fname(diri, prefix_in, ntap1, ndigit, &
                       trim(suffix_in)//trim(zpad(irank, rank_digits))), &
                       4_int64*rank_words)

    gi0 = max(ox0, rx*nxx+1) ; gi1 = min(ox0+onx-1, (rx+1)*nxx)
    gj0 = max(oy0, ry*nyy+1) ; gj1 = min(oy0+ony-1, (ry+1)*nyy)
    li0 = gi0 - rx*nxx ; li1 = gi1 - rx*nxx
    lj0 = gj0 - ry*nyy ; lj1 = gj1 - ry*nyy
    lk0 = gk0 - kz*nzz
    nj  = lj1 - lj0 + 1

    allocate(buf(nxx, nj, nk))
    if (nj == nyy) then
      ! the z planes needed are contiguous in the file -> read at once
      off = int(lk0-1,int64)*nxx*nyy
      nw  = int(nxx,int64)*nyy*nk
      call seq_read_part(sf, int(m,int64), off, buf, nw)
    else
      ! only part of y is needed -> read only the needed rows of each z plane
      do k = 1, nk
        off = int(lk0+k-2,int64)*nxx*nyy + int(lj0-1,int64)*nxx
        nw  = int(nxx,int64)*nj
        call seq_read_part(sf, int(m,int64), off, buf(:,:,k), nw)
      end do
    end if
    call seq_close(sf)

    slab(gi0-ox0+1:gi1-ox0+1, gj0-oy0+1:gj1-oy0+1, 1:nk, ib) = buf(li0:li1, :, :)
    deallocate(buf)

    !$omp atomic
    bytes_read = bytes_read + 4_int64*nxx*nj*nk
  end subroutine fill_block

  !----------------------------------------------------------------------------
  recursive subroutine write_slab(b, nk)
    integer, intent(in) :: b, nk
    integer :: k
    do k = 1, nk
      write(uo) slab(:, :, k, b)
    end do
  end subroutine write_slab

end program combine3d
