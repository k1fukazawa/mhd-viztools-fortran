! compare interior points (excluding boundaries) of old/new .vort. Points where the old one is NaN are excluded.
! The new code computes in double precision, so agreement within a relative tol of the maximum of each variable is accepted
program cmp_vort
  implicit none
  integer :: nx2, ny2, nz2, nv, m, k, nbad, nnan, u1, u2
  real, allocatable :: a(:,:,:), b(:,:,:)
  real :: amax, dmax
  real, parameter :: tol = 1.0e-4
  character(len=256) :: f1, f2, s
  call get_command_argument(1, f1) ; call get_command_argument(2, f2)
  call get_command_argument(3, s) ; read(s,*) nx2
  call get_command_argument(4, s) ; read(s,*) ny2
  call get_command_argument(5, s) ; read(s,*) nz2
  nv = 2
  allocate(a(nx2,ny2,nz2), b(nx2,ny2,nz2))
  open(newunit=u1, file=f1, form='unformatted', access='sequential', status='old')
  open(newunit=u2, file=f2, form='unformatted', access='sequential', status='old')
  nbad = 0 ; nnan = 0
  do m = 1, nv
    do k = 1, nz2
      read(u1) a(:,:,k) ; read(u2) b(:,:,k)
    end do
    if (any(b /= b)) nbad = nbad + 1
    associate(ai => a(2:nx2-1,2:ny2-1,2:nz2-1), bi => b(2:nx2-1,2:ny2-1,2:nz2-1))
      nnan = nnan + count(ai /= ai)
      amax = maxval(abs(ai), mask=(ai == ai))
      dmax = maxval(abs(ai - bi), mask=(ai == ai))
      print '(a,i0,a,es10.2)', '  var ', m, '  max |old-new| / max|old| = ', dmax/max(amax, tiny(1.0))
      if (dmax > tol*amax) nbad = nbad + count(ai == ai .and. abs(ai - bi) > tol*amax)
    end associate
  end do
  print '(a,i0,a,i0)', '  NaN points in old (excluded): ', nnan, '   mismatches: ', nbad
  if (nbad /= 0) error stop 1
end program cmp_vort
