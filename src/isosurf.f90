!==============================================================================
! isosurf : isosurface extraction (marching tetrahedra)
!   Each cube is split into 6 tetrahedra and the surface is extracted as triangles in each.
!   Vertices are made unique per grid edge (shared with the neighbouring triangles), so the
!   file is small and the normals (gradient of the scalar) give smooth shading.
!   No large table (256 cases of marching cubes) is needed and the surface has no holes.
!==============================================================================
module isosurf
  use mhd_common, only: int64, real32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: iso_extract, sort_index

  type :: tri_store
    integer :: n = 0
    integer, allocatable :: e(:,:)          ! (2, 3*ntri): vertex = linear indices of both ends of a grid edge
  end type tri_store

contains

  !----------------------------------------------------------------------------
  ! s(n1,n2,n3): values at the grid points, org: position of grid point (1,1,1), st: spacing
  ! output: pv(3,nv) positions, nrm(3,nv) unit normals (toward smaller values), tri(3,nt) (1-based)
  !----------------------------------------------------------------------------
  subroutine iso_extract(s, org, st, level, nv, pv, nrm, nt, tri)
    real(real32), intent(in) :: s(:,:,:)
    real, intent(in) :: org(3), st, level
    integer, intent(out) :: nv, nt
    real, allocatable, intent(out) :: pv(:,:), nrm(:,:)
    integer, allocatable, intent(out) :: tri(:,:)
    integer, parameter :: co(3,0:7) = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], [3,8])
    integer, parameter :: tet(4,6) = reshape([0,5,1,6, 0,1,2,6, 0,2,3,6, 0,3,7,6, 0,7,4,6, 0,4,5,6], [4,6])
    type(tri_store), allocatable :: lay(:)
    integer :: n(3), i, j, k, t, c, ins(4), nin, m, nk, ci(3,4)
    integer(int64), allocatable :: key(:), ord(:), uk(:)
    integer, allocatable :: inv(:)
    integer(int64) :: nn
    real :: f(4)
    n = shape(s)
    nn = int(n(1),int64)*n(2)*n(3)
    nv = 0 ; nt = 0
    allocate(lay(max(n(3)-1,1)))
    !$omp parallel do schedule(dynamic,1) private(i,j,t,c,ins,nin,f,ci)
    do k = 1, n(3)-1
      allocate(lay(k)%e(2, 3*1024)) ; lay(k)%n = 0
      do j = 1, n(2)-1
        do i = 1, n(1)-1
          do t = 1, 6
            do c = 1, 4
              ci(:,c) = [i, j, k] + co(:, tet(c,t))
              f(c) = s(ci(1,c), ci(2,c), ci(3,c)) - level
            end do
            if (any(ieee_is_nan(f))) cycle
            nin = 0
            do c = 1, 4
              if (f(c) >= 0.0) then
                nin = nin + 1 ; ins(nin) = c
              end if
            end do
            if (nin == 0 .or. nin == 4) cycle
            call emit(lay(k), nin, ins, ci)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    nk = 0
    do k = 1, n(3)-1
      nk = nk + lay(k)%n
    end do
    if (nk == 0) then
      allocate(pv(3,0), nrm(3,0), tri(3,0)) ; return
    end if
    allocate(key(nk))
    m = 0
    do k = 1, n(3)-1
      do i = 1, lay(k)%n
        m = m + 1
        key(m) = int(min(lay(k)%e(1,i), lay(k)%e(2,i)),int64)*nn + max(lay(k)%e(1,i), lay(k)%e(2,i))
      end do
    end do
    deallocate(lay)
    allocate(ord(nk), inv(nk), uk(nk))
    call sort_index(key, ord)
    do i = 1, nk
      if (i == 1) then
        nv = 1 ; uk(1) = key(ord(1))
      else if (key(ord(i)) /= key(ord(i-1))) then
        nv = nv + 1 ; uk(nv) = key(ord(i))
      end if
      inv(ord(i)) = nv
    end do
    nt = nk/3
    allocate(tri(3,nt))
    tri = reshape(inv, [3, nt])
    allocate(pv(3,nv), nrm(3,nv))
    !$omp parallel do
    do i = 1, nv
      call vertex(int(uk(i)/nn), int(mod(uk(i), nn)), pv(:,i), nrm(:,i))
    end do
    !$omp end parallel do
  contains
    integer function lin(p)
      integer, intent(in) :: p(3)
      lin = (p(1)-1) + n(1)*((p(2)-1) + n(2)*(p(3)-1))
    end function lin
    subroutine unlin(l, p)
      integer, intent(in) :: l
      integer, intent(out) :: p(3)
      p(1) = mod(l, n(1)) + 1
      p(2) = mod(l/n(1), n(2)) + 1
      p(3) = l/(n(1)*n(2)) + 1
    end subroutine unlin
    subroutine emit(L, nin, ins, ci)
      type(tri_store), intent(inout) :: L
      integer, intent(in) :: nin, ins(4), ci(3,4)
      integer :: outs(4), no, c, x(4)
      logical :: isin(4)
      isin = .false. ; isin(ins(1:nin)) = .true.
      no = 0
      do c = 1, 4
        if (.not. isin(c)) then
          no = no + 1 ; outs(no) = c
        end if
      end do
      do c = 1, 4
        x(c) = lin(ci(:,c))
      end do
      if (nin == 1) then
        call put(L, x(ins(1)), x(outs(1)), x(ins(1)), x(outs(2)), x(ins(1)), x(outs(3)))
      else if (nin == 3) then
        call put(L, x(outs(1)), x(ins(1)), x(outs(1)), x(ins(2)), x(outs(1)), x(ins(3)))
      else
        ! quad (i1,o1)-(i1,o2)-(i2,o2)-(i2,o1)
        call put(L, x(ins(1)), x(outs(1)), x(ins(1)), x(outs(2)), x(ins(2)), x(outs(2)))
        call put(L, x(ins(1)), x(outs(1)), x(ins(2)), x(outs(2)), x(ins(2)), x(outs(1)))
      end if
    end subroutine emit
    subroutine put(L, a1, b1, a2, b2, a3, b3)
      type(tri_store), intent(inout) :: L
      integer, intent(in) :: a1, b1, a2, b2, a3, b3
      integer, allocatable :: tmp(:,:)
      if (L%n + 3 > size(L%e,2)) then
        allocate(tmp(2, 2*size(L%e,2))) ; tmp(:,1:L%n) = L%e(:,1:L%n) ; call move_alloc(tmp, L%e)
      end if
      L%e(:,L%n+1) = [a1,b1] ; L%e(:,L%n+2) = [a2,b2] ; L%e(:,L%n+3) = [a3,b3]
      L%n = L%n + 3
    end subroutine put
    subroutine vertex(la, lb, x, g)
      integer, intent(in) :: la, lb
      real, intent(out) :: x(3), g(3)
      integer :: pa(3), pb(3)
      real :: fa, fb, w, ga(3), gb(3)
      call unlin(la, pa) ; call unlin(lb, pb)
      fa = s(pa(1),pa(2),pa(3)) - level ; fb = s(pb(1),pb(2),pb(3)) - level
      w = 0.5
      if (fa /= fb) w = fa/(fa - fb)
      w = min(max(w, 0.0), 1.0)
      x = org + st*(real(pa-1) + w*real(pb-pa))
      call grad(pa, ga) ; call grad(pb, gb)
      g = (1.0-w)*ga + w*gb
      if (sum(g**2) > 0.0) then
        g = -g/sqrt(sum(g**2))
      else
        g = [0.0, 0.0, 1.0]
      end if
    end subroutine vertex
    subroutine grad(p, g)
      integer, intent(in) :: p(3)
      real, intent(out) :: g(3)
      integer :: d, a(3), b(3)
      do d = 1, 3
        a = p ; b = p
        a(d) = max(1, p(d)-1) ; b(d) = min(n(d), p(d)+1)
        g(d) = (s(b(1),b(2),b(3)) - s(a(1),a(2),a(3)))/real(max(b(d)-a(d),1))
        if (ieee_is_nan(g(d))) g(d) = 0.0
      end do
    end subroutine grad
  end subroutine iso_extract

  !----------------------------------------------------------------------------
  ! index sort of 64-bit keys (iterative merge sort, stable)
  !----------------------------------------------------------------------------
  subroutine sort_index(key, ord)
    integer(int64), intent(in) :: key(:)
    integer(int64), intent(out) :: ord(:)
    integer(int64), allocatable :: tmp(:)
    integer(int64) :: n, w, lo, mid, hi, i, j, k
    n = size(key)
    ord = [(i, i=1,n)]
    allocate(tmp(n))
    w = 1
    do while (w < n)
      lo = 1
      do while (lo <= n)
        mid = min(lo + w - 1, n) ; hi = min(lo + 2*w - 1, n)
        i = lo ; j = mid + 1 ; k = lo
        do while (i <= mid .and. j <= hi)
          if (key(ord(i)) <= key(ord(j))) then
            tmp(k) = ord(i) ; i = i + 1
          else
            tmp(k) = ord(j) ; j = j + 1
          end if
          k = k + 1
        end do
        do while (i <= mid)
          tmp(k) = ord(i) ; i = i + 1 ; k = k + 1
        end do
        do while (j <= hi)
          tmp(k) = ord(j) ; j = j + 1 ; k = k + 1
        end do
        lo = lo + 2*w
      end do
      ord = tmp
      w = 2*w
    end do
  end subroutine sort_index

end module isosurf
