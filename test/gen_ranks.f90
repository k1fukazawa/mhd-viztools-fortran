program gen
  implicit none
  integer, parameter :: npex=2,npey=3,npez=4,nxx=10,nyy=8,nzz=6,nb=8
  real :: f(nxx,nyy,nzz)
  integer :: rx,ry,rz,irank,m,i,j,k,u
  character(len=64) :: fn
  do rz=0,npez-1; do ry=0,npey-1; do rx=0,npex-1
    irank = rx + npex*(ry+npey*rz)
    write(fn,'("in/ju03nt000.tdata.",i4.4)') irank
    open(newunit=u,file=fn,form='unformatted',access='sequential',status='replace')
    do m=1,nb
      do k=1,nzz; do j=1,nyy; do i=1,nxx
        f(i,j,k) = m*1.0e6 + (i+rx*nxx) + 100.0*(j+ry*nyy) + 1.0e4*(k+rz*nzz) + 0.25
      end do; end do; end do
      write(u) f
    end do
    close(u)
  end do; end do; end do
end program
