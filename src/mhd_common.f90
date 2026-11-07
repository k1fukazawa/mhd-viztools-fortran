!==============================================================================
! mhd_common : common module of the viztools post-processing tools
!
!   - timing          : wtime()
!   - file names      : make_fname(), zpad()
!   - stream access to Fortran sequential unformatted files (type seqfile)
!       * record-marker size (4/8 bytes) and endianness detected from the first record
!       * fixed-length records read/written directly at any position (pos=)
!       * partial reads of a record
!     so that only the planes and variables needed are read, instead of reading all
!     records sequentially into huge arrays as the legacy codes did.
!
! Modernized in 2026 from the visualization tools of the late Prof. Tatsuki Ogino and the
! post-processing tools of K. FUKAZAWA (2007-)
!==============================================================================
module mhd_common
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, &
                                           error_unit, output_unit
  implicit none
  private

  public :: int32, int64, real32, real64
  public :: wtime, make_fname, zpad, fatal, bswap4, bswap_r4
  public :: seqfile, seq_open_read, seq_open_write, seq_close
  public :: seq_read_rec, seq_read_part, seq_write_rec, seq_check_size
  public :: report_rate

  !--- sequential unformatted file accessed as a stream ----------------------
  type :: seqfile
    integer        :: unit   = -1
    integer        :: marker = 4        ! record marker size [byte] (4 or 8)
    logical        :: swap   = .false.  ! .true. if file endianness /= native
    integer(int64) :: reclen = 0        ! data bytes of one (fixed) record
    character(len=:), allocatable :: name
  end type seqfile

contains

  !----------------------------------------------------------------------------
  ! wall clock time [s] (the legacy code truncated sub-second time by integer division)
  !----------------------------------------------------------------------------
  function wtime() result(t)
    real(real64) :: t
    integer(int64) :: c, r
    call system_clock(c, r)
    t = real(c, real64) / real(r, real64)
  end function wtime

  subroutine report_rate(label, nbytes, sec)
    character(len=*), intent(in) :: label
    integer(int64),   intent(in) :: nbytes
    real(real64),     intent(in) :: sec
    real(real64) :: gb
    gb = real(nbytes, real64) / 1024.0_real64**3
    if (sec > 0.0_real64) then
      write(output_unit,'(a,": ",f10.2," s, ",f10.2," GiB, ",f8.3," GiB/s")') &
           trim(label), sec, gb, gb/sec
    else
      write(output_unit,'(a,": ",f10.2," s, ",f10.2," GiB")') trim(label), sec, gb
    end if
    flush(output_unit)
  end subroutine report_rate

  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    write(error_unit,'(a)') 'ERROR: '//trim(msg)
    error stop 1
  end subroutine fatal

  !----------------------------------------------------------------------------
  ! zero padded integer:  zpad(7,3) -> '007'
  !----------------------------------------------------------------------------
  ! The result is fixed length (left aligned); trim() it at the caller.
  ! (avoids internal write and deferred-length results, which are not thread safe in gfortran)
  pure function zpad(n, nd) result(s)
    integer, intent(in) :: n, nd
    character(len=16) :: s
    character(len=16) :: buf
    integer :: v, p
    v = abs(n) ; p = 16 ; buf = ' '
    do
      buf(p:p) = achar(iachar('0') + mod(v,10))
      v = v / 10 ; p = p - 1
      if (v == 0 .and. 16-p >= max(nd,1)) exit
      if (p == 1) exit
    end do
    if (n < 0) then
      buf(p:p) = '-' ; p = p - 1
    end if
    s = buf(p+1:16)
  end function zpad

  !----------------------------------------------------------------------------
  ! dir // prefix // zpad(n) // suffix   (fixed length 1024; use with trim)
  !----------------------------------------------------------------------------
  pure function make_fname(dir, prefix, n, nd, suffix) result(s)
    character(len=*), intent(in) :: dir, prefix, suffix
    integer,          intent(in) :: n, nd
    character(len=1024) :: s
    s = trim(dir)//trim(prefix)//trim(zpad(n,nd))//trim(suffix)
  end function make_fname

  !----------------------------------------------------------------------------
  ! byte swap
  !----------------------------------------------------------------------------
  elemental function bswap4(i) result(o)
    integer(int32), intent(in) :: i
    integer(int32) :: o
    o = 0
    call mvbits(i,  0, 8, o, 24)
    call mvbits(i,  8, 8, o, 16)
    call mvbits(i, 16, 8, o,  8)
    call mvbits(i, 24, 8, o,  0)
  end function bswap4

  elemental function bswap8(i) result(o)
    integer(int64), intent(in) :: i
    integer(int64) :: o
    integer :: b
    o = 0
    do b = 0, 7
      call mvbits(i, 8*b, 8, o, 56-8*b)
    end do
  end function bswap8

  subroutine bswap_r4(a, n)
    integer(int64), intent(in)    :: n
    real(real32),   intent(inout) :: a(n)
    integer(int64) :: i
    integer(int32) :: t
    !$omp parallel do private(t) if(n > 1000000_int64)
    do i = 1, n
      t    = transfer(a(i), t)
      a(i) = transfer(bswap4(t), a(i))
    end do
  end subroutine bswap_r4

  !----------------------------------------------------------------------------
  ! Open an existing sequential-unformatted file as stream and detect the
  ! record-marker convention from the first record header.
  !   reclen_bytes : expected data length of every record [byte]
  !----------------------------------------------------------------------------
  subroutine seq_open_read(sf, fname, reclen_bytes)
    type(seqfile),    intent(out) :: sf
    character(len=*), intent(in)  :: fname
    integer(int64),   intent(in)  :: reclen_bytes
    integer :: ios
    integer(int32) :: h4
    character(len=256) :: msg

    sf%name   = trim(fname)
    sf%reclen = reclen_bytes
    open(newunit=sf%unit, file=sf%name, access='stream', form='unformatted', &
         status='old', action='read', iostat=ios, iomsg=msg)
    if (ios /= 0) call fatal('cannot open '//sf%name//' : '//trim(msg))

    ! accept the combination for which both the leading and trailing markers agree
    ! (the low 4 bytes of an 8-byte marker cannot be told from a 4-byte marker)
    read(sf%unit, pos=1, iostat=ios) h4
    if (ios /= 0) call fatal('cannot read header of '//sf%name)
    if (try4(.false.)) return
    if (try4(.true.))  return
    if (try8(.false.)) return
    if (try8(.true.))  return
    write(msg,'(a,i0,a,i0)') 'record length mismatch: expected ', reclen_bytes, &
                             ' bytes, header = ', h4
    call fatal(trim(msg)//' in '//sf%name//' (grid size setting wrong?)')
  contains
    logical function try4(sw)
      logical, intent(in) :: sw
      integer(int32) :: a, b
      integer :: io1, io2
      try4 = .false.
      read(sf%unit, pos=1, iostat=io1) a
      read(sf%unit, pos=5+reclen_bytes, iostat=io2) b
      if (io1 /= 0 .or. io2 /= 0) return
      if (sw) then
        a = bswap4(a) ; b = bswap4(b)
      end if
      if (int(a,int64) == reclen_bytes .and. int(b,int64) == reclen_bytes) then
        try4 = .true. ; sf%marker = 4 ; sf%swap = sw
      end if
    end function try4
    logical function try8(sw)
      logical, intent(in) :: sw
      integer(int64) :: a, b
      integer :: io1, io2
      try8 = .false.
      read(sf%unit, pos=1, iostat=io1) a
      read(sf%unit, pos=9+reclen_bytes, iostat=io2) b
      if (io1 /= 0 .or. io2 /= 0) return
      if (sw) then
        a = bswap8(a) ; b = bswap8(b)
      end if
      if (a == reclen_bytes .and. b == reclen_bytes) then
        try8 = .true. ; sf%marker = 8 ; sf%swap = sw
      end if
    end function try8
  end subroutine seq_open_read

  !----------------------------------------------------------------------------
  ! Create a sequential-unformatted-compatible file written through stream
  ! access (so that records can be written in any order).
  !----------------------------------------------------------------------------
  subroutine seq_open_write(sf, fname, reclen_bytes, marker, swap)
    type(seqfile),    intent(out) :: sf
    character(len=*), intent(in)  :: fname
    integer(int64),   intent(in)  :: reclen_bytes
    integer,          intent(in)  :: marker
    logical,          intent(in)  :: swap
    integer :: ios
    character(len=256) :: msg
    sf%name   = trim(fname)
    sf%reclen = reclen_bytes
    sf%marker = marker
    sf%swap   = swap
    open(newunit=sf%unit, file=sf%name, access='stream', form='unformatted', &
         status='replace', action='readwrite', iostat=ios, iomsg=msg)
    if (ios /= 0) call fatal('cannot create '//sf%name//' : '//trim(msg))
  end subroutine seq_open_write

  subroutine seq_close(sf)
    type(seqfile), intent(inout) :: sf
    if (sf%unit /= -1) close(sf%unit)
    sf%unit = -1
  end subroutine seq_close

  ! byte position (1-based) of the first data byte of record irec (1-based)
  pure function recpos(sf, irec) result(p)
    type(seqfile),  intent(in) :: sf
    integer(int64), intent(in) :: irec
    integer(int64) :: p
    p = (irec-1_int64)*(sf%reclen + 2_int64*sf%marker) + sf%marker + 1_int64
  end function recpos

  !----------------------------------------------------------------------------
  ! check that the file size equals nrec fixed-length records
  !----------------------------------------------------------------------------
  subroutine seq_check_size(sf, nrec)
    type(seqfile),  intent(in) :: sf
    integer(int64), intent(in) :: nrec
    integer(int64) :: fsize, expect
    character(len=200) :: msg
    inquire(unit=sf%unit, size=fsize)
    expect = nrec*(sf%reclen + 2_int64*sf%marker)
    if (fsize >= 0 .and. fsize < expect) then
      write(msg,'(a,i0,a,i0)') 'file too small: size=', fsize, ' expected>=', expect
      call fatal(trim(msg)//' : '//sf%name)
    end if
  end subroutine seq_check_size

  !----------------------------------------------------------------------------
  ! read whole record irec into buf(n)  (n*4 must equal reclen)
  !----------------------------------------------------------------------------
  subroutine seq_read_rec(sf, irec, buf, n)
    type(seqfile),  intent(in)  :: sf
    integer(int64), intent(in)  :: irec, n
    real(real32),   intent(out) :: buf(n)
    integer :: ios
    read(sf%unit, pos=recpos(sf,irec), iostat=ios) buf
    if (ios /= 0) call fatal('read error in '//sf%name)
    if (sf%swap) call bswap_r4(buf, n)
  end subroutine seq_read_rec

  !----------------------------------------------------------------------------
  ! read n words starting at word offset off (0-based) inside record irec
  !----------------------------------------------------------------------------
  subroutine seq_read_part(sf, irec, off, buf, n)
    type(seqfile),  intent(in)  :: sf
    integer(int64), intent(in)  :: irec, off, n
    real(real32),   intent(out) :: buf(n)
    integer :: ios
    read(sf%unit, pos=recpos(sf,irec)+4_int64*off, iostat=ios) buf
    if (ios /= 0) call fatal('read error (partial) in '//sf%name)
    if (sf%swap) call bswap_r4(buf, n)
  end subroutine seq_read_part

  !----------------------------------------------------------------------------
  ! write record irec (with leading/trailing markers)
  !----------------------------------------------------------------------------
  subroutine seq_write_rec(sf, irec, buf, n)
    type(seqfile),  intent(in)    :: sf
    integer(int64), intent(in)    :: irec, n
    real(real32),   intent(inout) :: buf(n)   ! swapped temporarily if needed
    integer(int64) :: p
    integer(int32) :: m4
    integer(int64) :: m8
    p = recpos(sf, irec) - sf%marker
    if (sf%marker == 4) then
      m4 = int(sf%reclen, int32)
      if (sf%swap) m4 = bswap4(m4)
      write(sf%unit, pos=p) m4
    else
      m8 = sf%reclen
      if (sf%swap) m8 = bswap8(m8)
      write(sf%unit, pos=p) m8
    end if
    if (sf%swap) call bswap_r4(buf, n)
    write(sf%unit) buf
    if (sf%swap) call bswap_r4(buf, n)
    if (sf%marker == 4) then
      write(sf%unit) m4
    else
      write(sf%unit) m8
    end if
  end subroutine seq_write_rec

end module mhd_common
