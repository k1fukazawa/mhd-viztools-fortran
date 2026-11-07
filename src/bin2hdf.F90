!==============================================================================
! bin2hdf : convert combined binaries (m.data / .vort) to HDF5 + XDMF for ParaView
!
!   Merges the legacy HDF5 converters for MHD data (8 variables) and vorticity (2 variables).
!   Number and names of variables and file names are given in a namelist (namelist/bin2hdf_*.nml).
!   HDF5 dataset names and array shape (nz2 x ny2 x nx2 in C order) are the same as before,
!   so existing Python (h5py) scripts and ParaView settings still work.
!
! Changes from the legacy code
!   (1) legacy: a whole variable (24 GB) was read into memory before writing
!       new   : written as hyperslabs of nslab planes (default 32 planes = 0.75 GB x 2 buffers)
!   (2) the file data type is native (little endian) by default
!       legacy: H5T_IEEE_F32BE -> byte swapping of all data on writing
!       filetype='be' selects the legacy type. h5py and ParaView convert automatically,
!       so readers do not need to care.
!   (3) OpenMP: reading the next slab overlaps the HDF5 write (overlap)
!       HDF5 is always called by one thread only (no thread-safe HDF5 needed)
!   (4) out_format='raw' : builds and runs without the HDF5 library.
!       Writes raw binary (.raw) and an XDMF file referring to it (readable by ParaView)
!   (5) origin='planet' writes XDMF with the planet center at the origin (in planetary radii)
!       (default 'legacy' keeps the Origin values of the legacy code)
!   (6) values are converted to physical units by &units (mhdunits) (default).
!       phys_units=.false. keeps the normalized values of the legacy code.
!       The unit of each variable is written to the HDF5 attributes units / scale_from_normalized
!       and to an XDMF Information element. The kind is varunit (guessed from varname if empty)
!
! Build: with HDF5    -> compile with -DUSE_HDF5 using h5fc etc. (see Makefile)
!        without HDF5 -> compile as is (raw mode only)
! Usage:  ./bin2hdf [namelist file]     (default: bin2hdf.nml)
!==============================================================================
program bin2hdf_main
  use mhd_common
  use mhdunits
#ifdef USE_HDF5
  use hdf5
#endif
  implicit none

  integer, parameter :: maxvar = 32
  !---- namelist --------------------------------------------------------------
  integer :: nx2 = 3002, ny2 = 2002, nz2 = 1002, nvar = 8
  character(len=32) :: varname(maxvar) = ''
  character(len=8)  :: varunit(maxvar) = ''     ! n, v, p, b, omega, ... (empty: guessed from varname)
  character(len=512) :: diri = './', diro = './'
  character(len=64)  :: prefix = 'ju02nt', prefix_out = ''
  character(len=64)  :: suffix_in = 'm.data', suffix_out = '', suffix_xmf = '.xmf'
  integer :: ntap = 181, nfiles = 1, nstep = 1, ndigit = 3
  real    :: h = 0.15
  integer :: nxp = 500
  character(len=16) :: origin = 'legacy'      ! 'legacy' or 'planet'
  real    :: i0 = -1.0, j0 = -1.0, k0 = -1.0  ! grid position of the planet center (negative: automatic)
  character(len=16) :: filetype = 'native'    ! 'native' or 'be' (data type in HDF5)
  character(len=16) :: out_format = 'hdf5'    ! 'hdf5' or 'raw'
  integer :: nslab = 32                        ! number of z planes written at once
  integer :: compress = 0                      ! gzip level (0 = none; needs chunking)
  logical :: overlap = .true.
  namelist /bin2hdf/ nx2, ny2, nz2, nvar, varname, varunit, diri, diro, prefix,    &
                     prefix_out, suffix_in, suffix_out, suffix_xmf, ntap, &
                     nfiles, nstep, ndigit, h, nxp, origin, i0, j0, k0,   &
                     filetype, out_format, nslab, compress, overlap

  real(real32), allocatable :: buf(:,:,:,:)          ! (nx2,ny2,nslab,2)
  real(real64) :: vfac(maxvar) = 1.0_real64          ! normalized -> physical units
  character(len=12) :: vlab(maxvar) = ''
  type(seqfile) :: fin
  integer :: u, ios, ifile, ntap1, nsl, ntask, t, cur, nxt, uraw
  integer(int64) :: nplane, recb
  real(real64) :: ta, tr, tw, t0
  character(len=:), allocatable :: fni, fno, fnx, fnh, nmlfile
  character(len=512) :: arg
  character(len=256) :: msg
  logical :: is_raw
#ifdef USE_HDF5
  integer(HID_T) :: file_id, dset_id, fspace, mspace, dcpl, ftype
  integer :: herr
  integer(HSIZE_T) :: dims(3), cnt(3), offs(3), chunk(3)
#endif

  !---- namelist ---------------------------------------------------------------
  varname(1:8) = ['RHO', 'VX ', 'VY ', 'VZ ', 'P  ', 'BX ', 'BY ', 'BZ ']
  nmlfile = 'bin2hdf.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg) ; nmlfile = trim(arg)
  end if
  open(newunit=u, file=nmlfile, status='old', action='read', iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('cannot open namelist file '//nmlfile//' : '//trim(msg))
  read(u, nml=bin2hdf, iostat=ios, iomsg=msg)
  if (ios /= 0) call fatal('namelist &bin2hdf : '//trim(msg))
  call units_read(u)
  close(u)
  if (nvar < 1 .or. nvar > maxvar) call fatal('nvar out of range')
  if (any(len_trim(varname(1:nvar)) == 0)) call fatal('varname(1:nvar) must be given')
  do u = 1, nvar
    if (len_trim(varunit(u)) == 0) varunit(u) = kind_of_name(varname(u))
    vfac(u) = ufac(varunit(u))
    vlab(u) = ulabel(varunit(u))
  end do
  if (len_trim(prefix_out) == 0) prefix_out = prefix
  is_raw = (trim(out_format) == 'raw')
  if (len_trim(suffix_out) == 0) suffix_out = merge('.raw', '.h5 ', is_raw)
  if (i0 < 0.0) i0 = real(nx2/2 - nxp) + 0.5
  if (j0 < 0.0) j0 = real(ny2/2) + 0.5
  if (k0 < 0.0) k0 = 1.5
  nslab = max(1, min(nslab, nz2))
#ifndef USE_HDF5
  if (.not. is_raw) call fatal('built without HDF5: use out_format=''raw'' or rebuild with -DUSE_HDF5')
#endif

  nplane = int(nx2,int64)*ny2
  recb   = 4_int64*nplane
  nsl    = (nz2 + nslab - 1)/nslab           ! slabs per variable
  ntask  = nsl*nvar
  allocate(buf(nx2, ny2, nslab, 2))

  write(*,'(a)')        '----- bin2hdf -----'
  write(*,'(a,3i6,a,i0)') ' grid: ', nx2, ny2, nz2, '   nvar=', nvar
  write(*,'(a,a,a,a)')  ' format: ', trim(out_format), '   filetype: ', trim(filetype)
  call units_report()
  do u = 1, nvar
    write(*,'(3x,a12,a,a,a,es12.5)') varname(u), ' [', trim(vlab(u)), ']  x', vfac(u)
  end do
  write(*,'(a,i0,a,f8.2,a)') ' slab: ', nslab, ' planes (', 8.0d0*recb*nslab/1024d0**3, ' GiB buffers)'

#ifdef USE_HDF5
  if (.not. is_raw) call h5open_f(herr)
#endif

  if (nfiles < 1) call fatal('no input file: nfiles must be >= 1')
  do ifile = 1, nfiles
    ntap1 = ntap + (ifile-1)*nstep
    fni = trim(make_fname(diri, prefix, ntap1, ndigit, suffix_in))
    fno = trim(make_fname(diro, prefix_out, ntap1, ndigit, suffix_out))
    fnh = trim(make_fname('',   prefix_out, ntap1, ndigit, suffix_out))
    fnx = trim(make_fname(diro, prefix_out, ntap1, ndigit, suffix_xmf))
    write(*,'(a)') ' input : '//fni
    write(*,'(a)') ' output: '//fno
    ta = wtime() ; tr = 0 ; tw = 0

    call seq_open_read(fin, fni, recb)
    call seq_check_size(fin, int(nvar,int64)*nz2)
    call open_output()

    ! ---- pipeline: read(t+1) || write(t) ----
    t0 = wtime()
    call read_task(1, 1)
    tr = tr + (wtime()-t0)
    !$omp parallel default(shared)
    !$omp single
    do t = 1, ntask
      cur = mod(t-1,2) + 1
      nxt = 3 - cur
      if (t < ntask) then
        if (overlap) then
          !$omp task default(shared) firstprivate(t, nxt)
          call read_task(t+1, nxt)
          !$omp end task
        end if
      end if
      t0 = wtime()
      call write_task(t, cur)
      tw = tw + (wtime()-t0)
      if (t < ntask .and. .not. overlap) then
        t0 = wtime()
        call read_task(t+1, nxt)
        tr = tr + (wtime()-t0)
      end if
      !$omp taskwait
    end do
    !$omp end single
    !$omp end parallel

    call close_output()
    call seq_close(fin)
    call write_xmf()
    call report_rate('  write (HDF5/raw, excl. overlapped read)', &
                     int(nvar,int64)*recb*nz2, tw)
    call report_rate('  total', 2_int64*nvar*recb*nz2, wtime()-ta)
  end do

#ifdef USE_HDF5
  if (.not. is_raw) call h5close_f(herr)
#endif
  deallocate(buf)

contains

  ! task t -> variable m and slab s
  pure subroutine task_index(t, m, kb, nk)
    integer, intent(in)  :: t
    integer, intent(out) :: m, kb, nk
    integer :: s
    m  = (t-1)/nsl + 1
    s  = mod(t-1, nsl) + 1
    kb = (s-1)*nslab + 1
    nk = min(nslab, nz2-kb+1)
  end subroutine task_index

  !----------------------------------------------------------------------------
  recursive subroutine read_task(t, ib)
    integer, intent(in) :: t, ib
    integer :: m, kb, nk, k
    call task_index(t, m, kb, nk)
    do k = 1, nk
      call seq_read_rec(fin, int(m-1,int64)*nz2 + kb + k - 1, buf(:,:,k,ib), nplane)
      if (vfac(m) /= 1.0_real64) buf(:,:,k,ib) = real(buf(:,:,k,ib)*vfac(m), real32)
    end do
  end subroutine read_task

  !----------------------------------------------------------------------------
  subroutine open_output()
    if (is_raw) then
      open(newunit=uraw, file=fno, access='stream', form='unformatted', &
           status='replace', action='write', iostat=ios, iomsg=msg)
      if (ios /= 0) call fatal('cannot create '//fno//' : '//trim(msg))
      return
    end if
#ifdef USE_HDF5
    call h5fcreate_f(fno, H5F_ACC_TRUNC_F, file_id, herr)
    if (herr /= 0) call fatal('h5fcreate failed: '//fno)
    dims = [int(nx2,HSIZE_T), int(ny2,HSIZE_T), int(nz2,HSIZE_T)]
    if (trim(filetype) == 'be') then
      ftype = H5T_IEEE_F32BE
    else
      ftype = H5T_NATIVE_REAL
    end if
#endif
  end subroutine open_output

  subroutine close_output()
    if (is_raw) then
      close(uraw) ; return
    end if
#ifdef USE_HDF5
    call h5fclose_f(file_id, herr)
#endif
  end subroutine close_output

  !----------------------------------------------------------------------------
  subroutine write_task(t, ib)
    integer, intent(in) :: t, ib
    integer :: m, kb, nk
    integer(int64) :: pos
    call task_index(t, m, kb, nk)
    if (is_raw) then
      pos = 1_int64 + 4_int64*nplane*(int(m-1,int64)*nz2 + kb - 1)
      write(uraw, pos=pos) buf(:,:,1:nk,ib)
      return
    end if
#ifdef USE_HDF5
    if (kb == 1) then     ! first slab of variable m -> create dataset
      call h5screate_simple_f(3, dims, fspace, herr)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, herr)
      if (compress > 0) then
        chunk = [int(nx2,HSIZE_T), int(ny2,HSIZE_T), 1_HSIZE_T]
        call h5pset_chunk_f(dcpl, 3, chunk, herr)
        call h5pset_deflate_f(dcpl, compress, herr)
      end if
      call h5dcreate_f(file_id, trim(varname(m)), ftype, fspace, dset_id, herr, dcpl)
      if (herr /= 0) call fatal('h5dcreate failed: '//trim(varname(m)))
      call h5pclose_f(dcpl, herr)
      call put_attrs(dset_id, m)
      write(*,'(a)') '   write '//trim(varname(m))
    end if
    cnt  = [int(nx2,HSIZE_T), int(ny2,HSIZE_T), int(nk,HSIZE_T)]
    offs = [0_HSIZE_T, 0_HSIZE_T, int(kb-1,HSIZE_T)]
    call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, offs, cnt, herr)
    call h5screate_simple_f(3, cnt, mspace, herr)
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, buf(:,:,1:nk,ib), cnt, herr, &
                    mem_space_id=mspace, file_space_id=fspace)
    if (herr /= 0) call fatal('h5dwrite failed')
    call h5sclose_f(mspace, herr)
    if (kb + nk - 1 == nz2) then   ! last slab
      call h5dclose_f(dset_id, herr)
      call h5sclose_f(fspace, herr)
    end if
#endif
  end subroutine write_task

#ifdef USE_HDF5
  ! unit attributes: units (string), scale_from_normalized (double)
  subroutine put_attrs(did, m)
    integer(HID_T), intent(in) :: did
    integer, intent(in) :: m
    integer(HID_T) :: asp, aty, aid
    integer(HSIZE_T) :: ad(1)
    character(len=32) :: us
    ad = 1
    us = vlab(m)
    if (len_trim(us) == 0) us = '-'
    call h5screate_f(H5S_SCALAR_F, asp, herr)
    call h5tcopy_f(H5T_FORTRAN_S1, aty, herr)
    call h5tset_size_f(aty, int(len_trim(us), SIZE_T), herr)
    call h5acreate_f(did, 'units', aty, asp, aid, herr)
    call h5awrite_f(aid, aty, us(1:len_trim(us)), ad, herr)
    call h5aclose_f(aid, herr)
    call h5tclose_f(aty, herr)
    call h5acreate_f(did, 'scale_from_normalized', H5T_NATIVE_DOUBLE, asp, aid, herr)
    call h5awrite_f(aid, H5T_NATIVE_DOUBLE, vfac(m), ad, herr)
    call h5aclose_f(aid, herr)
    call h5sclose_f(asp, herr)
  end subroutine put_attrs
#endif

  !----------------------------------------------------------------------------
  ! XDMF (same structure as the legacy code; variables = varname in the namelist)
  !----------------------------------------------------------------------------
  subroutine write_xmf()
    integer :: ux, m
    character(len=64) :: sdim, sorg, sspc, endian
    character(len=32) :: a, b, c
    integer(int64) :: seek
    write(sdim,'(i0,1x,i0,1x,i0)') nz2, ny2, nx2
    if (trim(origin) == 'planet') then
      ! coordinates with the planet center at the origin [planetary radii] (order z, y, x)
      write(a,'(f0.4)') (1.0-k0)*h ; write(b,'(f0.4)') (1.0-j0)*h ; write(c,'(f0.4)') (1.0-i0)*h
    else
      ! same values as the legacy code (grid indices used as Origin)
      write(a,'(f0.2)') k0 ; write(b,'(f0.2)') j0 ; write(c,'(f0.2)') i0
    end if
    sorg = trim(a)//' '//trim(b)//' '//trim(c)
    write(a,'(f0.4)') h
    if (a(1:1) == '.') a = '0'//trim(a)
    sspc = trim(a)//' '//trim(a)//' '//trim(a)
    endian = merge('Little', 'Big   ', little_endian())

    open(newunit=ux, file=fnx, form='formatted', status='replace')
    write(ux,'(a)') '<?xml version="1.0" ?>'
    write(ux,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(ux,'(a)') '<Xdmf Version="2.0">'
    write(ux,'(a)') '<Domain>'
    write(ux,'(a)') '<Grid Name="Structured Grid" GridType="Uniform">'
    write(ux,'(a)') '<Topology TopologyType="3DCORECTMesh" NumberOfElements="'//trim(sdim)//'"/>'
    write(ux,'(a)') '<Geometry GeometryType="ORIGIN_DXDYDZ">'
    write(ux,'(a)') '<DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="4" Format="XML">'
    write(ux,'(a)') trim(sorg)
    write(ux,'(a)') '</DataItem>'
    write(ux,'(a)') '<DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="4" Format="XML">'
    write(ux,'(a)') trim(sspc)
    write(ux,'(a)') '</DataItem>'
    write(ux,'(a)') '</Geometry>'
    do m = 1, nvar
      write(ux,'(a)') '<Attribute Name="'//trim(varname(m))//'" AttributeType="Scalar" Center="Node">'
      if (is_raw) then
        seek = 4_int64*nplane*nz2*(m-1)
        write(a,'(i0)') seek
        write(ux,'(a)') '<DataItem Dimensions="'//trim(sdim)//'" NumberType="Float" Precision="4"'// &
                        ' Format="Binary" Endian="'//trim(endian)//'" Seek="'//trim(a)//'">'
        write(ux,'(a)') fnh
      else
        write(ux,'(a)') '<DataItem Dimensions="'//trim(sdim)//'" NumberType="Float" Precision="4" Format="HDF">'
        write(ux,'(a)') fnh//':/'//trim(varname(m))
      end if
      write(ux,'(a)') '</DataItem>'
      if (len_trim(vlab(m)) > 0) &
        write(ux,'(a)') '<Information Name="units" Value="'//trim(vlab(m))//'"/>'
      write(ux,'(a)') '</Attribute>'
    end do
    write(ux,'(a)') '</Grid>'
    write(ux,'(a)') '</Domain>'
    write(ux,'(a)') '</Xdmf>'
    close(ux)
  end subroutine write_xmf

  logical function little_endian()
    integer(int32) :: one = 1
    character(len=4) :: c4
    c4 = transfer(one, c4)
    little_endian = (iachar(c4(1:1)) == 1)
  end function little_endian

end program bin2hdf_main
