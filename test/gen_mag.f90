! simple magnetosphere-like data for plot tests (302x202x102, h=1.5)
program gen_mag
  implicit none
  integer, parameter :: nx2=302, ny2=202, nz2=102, nxp=50
  real, parameter :: h=1.5
  real :: f(nx2,ny2,8), x, y, z, r, rc, rm, w, xm
  real :: i0, j0, k0
  integer :: i, j, k, m, u
  i0 = nx2/2-nxp+0.5 ; j0 = ny2/2+0.5 ; k0 = 1.5
  xm = 5.0e-2
  open(newunit=u, file='ju02nt171m.data', form='unformatted', access='stream', status='replace')
  do m = 1, 8
    do k = 1, nz2
      z = (k-k0)*h
      do j = 1, ny2
        y = (j-j0)*h
        do i = 1, nx2
          x = (i-i0)*h
          r = sqrt(x*x+y*y+z*z) ; rc = max(r, 7.0)
          w = exp(-(rc/60.0)**2)
          select case (m)
          case (1) ; f(i,j,m) = 1e-4*(10/rc)**2 + 1e-6
          case (2) ; f(i,j,m) = -1.38e-4*y*w + 3.27e-3*(1-w)
          case (3) ; f(i,j,m) =  1.38e-4*x*w + 1e-3*sin(y/30)*(1-w)
          case (4) ; f(i,j,m) =  1e-3*z/rc*w
          case (5) ; f(i,j,m) = (1e-4*(10/rc)**2 + 1e-6)*1.5e-6*exp(-rc/80)*(1+0.5*sin(x/20))
          case (6) ; f(i,j,m) = -3*xm*x*z/rc**5
          case (7) ; f(i,j,m) = -3*xm*y*z/rc**5
          case (8) ; f(i,j,m) = -xm*(3*z*z-rc*rc)/rc**5
          end select
        end do
      end do
      write(u) int(4*nx2*ny2), f(:,:,m), int(4*nx2*ny2)
    end do
  end do
  close(u)
end program gen_mag
