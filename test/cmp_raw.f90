! check that m.data (sequential) and the raw output of bin2hdf (stream) contain the same data
program cmp_raw
  implicit none
  integer :: nx2, ny2, nz2, nv, m, k, u1, u2, nbad
  real, allocatable :: a(:,:), b(:,:)
  character(len=256) :: f1, f2, s
  call get_command_argument(1, f1) ; call get_command_argument(2, f2)
  call get_command_argument(3, s) ; read(s,*) nx2
  call get_command_argument(4, s) ; read(s,*) ny2
  call get_command_argument(5, s) ; read(s,*) nz2
  call get_command_argument(6, s) ; read(s,*) nv
  allocate(a(nx2,ny2), b(nx2,ny2))
  open(newunit=u1, file=f1, form='unformatted', access='sequential', status='old')
  open(newunit=u2, file=f2, form='unformatted', access='stream', status='old')
  nbad = 0
  do m = 1, nv
    do k = 1, nz2
      read(u1) a ; read(u2) b
      nbad = nbad + count(a /= b)
    end do
  end do
  print '(a,i0)', '  mismatches: ', nbad
  if (nbad /= 0) error stop 1
end program cmp_raw
