program genm
  implicit none
  integer, parameter :: nx2=40, ny2=30, nz2=20, nb=8
  real :: f(nx2,ny2)
  integer :: i,j,k,m,u
  open(newunit=u,file='ju02nt186m.data',form='unformatted',access='sequential',status='replace')
  do m=1,nb; do k=1,nz2
    do j=1,ny2; do i=1,nx2
      f(i,j) = sin(0.3*i*m+0.2*j)*cos(0.17*k+m) + 0.1*m
    end do; end do
    if (m==2 .and. k==5) f(10:15,10:15) = 1.0   ! uniform flow patch -> omega=0
    if (m>=6 .and. k==7) f(5:8,5:8) = 0.0         ! B=0 patch
    write(u) f
  end do; end do
  close(u)
end program
