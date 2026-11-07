!==============================================================================
! scene3d : build a 3-D scene (lines, triangle meshes, textured quads, text) and write it
!           as binary glTF 2.0 (.glb) and X3D (.x3d)
!
!   Successor of the legacy VRML subroutines (zvrsub*.f).
!     linebe/linep/lineen      -> add_polyline   (all field lines in one LINES primitive)
!     pixel1/pixebe/.../pixeen -> add_image_quad (PNG texture with transparency)
!     csuba (isosurface)       -> add_mesh       (triangles with normals and vertex colours)
!     sphere1, symblv2         -> add_sphere, add_label
!   Positions are given in display coordinates (X, Y, Z: Z north) and written to the Y-up
!   glTF / X3D as (X, Z, -Y) (a rotation keeping right-handedness; the legacy VRML (x,z,y) was a mirror image).
!   .glb is binary, about 1/5-1/10 of the size of VRML (ASCII), and loads fast.
!==============================================================================
module scene3d
  use mhd_common, only: int32, int64, real32, fatal
  use, intrinsic :: iso_fortran_env, only: int8
  implicit none
  private
  public :: scene_t, add_polyline, add_mesh_begin, add_triangles, add_image_quad, &
            add_sphere, add_label, write_glb, write_x3d, png_encode, u8, scene_reset, &
            add_image_file, add_textured_sphere, add_legend, set_extras, mesh_uv

  type :: lines_t
    character(len=64) :: name = ''
    integer :: nv = 0, ns = 0
    real(real32), allocatable :: p(:,:)        ! (3, nv)
    integer(int8), allocatable :: c(:,:)       ! (4, nv) RGBA
    integer(int32), allocatable :: e(:,:)      ! (2, ns) 0-based
  end type lines_t

  type :: mesh_t
    character(len=64) :: name = ''
    integer :: nv = 0, nt = 0
    real :: rgba(4) = 1.0
    logical :: lit = .true.
    real(real32), allocatable :: p(:,:), n(:,:)
    integer(int8), allocatable :: c(:,:)
    integer(int32), allocatable :: t(:,:)      ! (3, nt) 0-based
    integer :: tex = 0                          ! >0: image (scene images) mapped on the mesh
    real(real32), allocatable :: uv(:,:)       ! (2, nv) texture coordinates
    character(len=256) :: extras = ''           ! glTF node extras (JSON object)
  end type mesh_t

  type :: quad_t
    character(len=64) :: name = ''
    real :: corner(3,4) = 0.0                  ! lower left, lower right, upper right, upper left
    integer :: w = 0, h = 0
    integer(int8), allocatable :: img(:,:,:)   ! (4, w, h) row 1 = top
    character(len=256) :: extras = ''           ! glTF node extras (JSON object)
  end type quad_t

  ! image file (JPEG/PNG) embedded as is
  type :: image_t
    character(len=16) :: mime = ''
    character(len=256) :: path = ''
    integer(int8), allocatable :: d(:)
  end type image_t

  ! colour bar (shown by the viewer while a layer whose name starts with target is visible)
  type :: legend_t
    character(len=64) :: target = '', label = ''
    real :: vmin = 0.0, vmax = 1.0
    real :: rgb(3,33) = 0.0
  end type legend_t

  type :: label_t
    character(len=64) :: text = ''
    real :: p(3) = 0.0, size = 1.0, rgb(3) = 1.0
  end type label_t

  type :: scene_t
    character(len=128) :: title = ''
    real :: bg(3) = 0.0
    integer :: nl = 0, nm = 0, nq = 0, nlab = 0, nimf = 0, nleg = 0
    type(lines_t), allocatable :: l(:)
    type(mesh_t),  allocatable :: m(:)
    type(quad_t),  allocatable :: q(:)
    type(label_t), allocatable :: lab(:)
    type(image_t)  :: imf(8)
    type(legend_t) :: leg(64)
  end type scene_t

  ! glb builder state
  integer(int8), allocatable :: bin(:)
  integer(int64) :: nbin
  character(len=:), allocatable :: jacc, jbv, jmesh, jnode, jmat, jimg, jtex
  integer :: nacc, nbv, nmeshes, nmat, nimg
  real :: scene_radius = 100.0

contains

  elemental function u8(v) result(b)
    integer, intent(in) :: v
    integer(int8) :: b
    integer :: w
    w = min(max(v, 0), 255)
    if (w > 127) w = w - 256
    b = int(w, int8)
  end function u8

  pure integer function uval(b)
    integer(int8), intent(in) :: b
    uval = int(b)
    if (uval < 0) uval = uval + 256
  end function uval

  subroutine scene_reset(sc)
    type(scene_t), intent(inout) :: sc
    if (allocated(sc%l)) deallocate(sc%l)
    if (allocated(sc%m)) deallocate(sc%m)
    if (allocated(sc%q)) deallocate(sc%q)
    if (allocated(sc%lab)) deallocate(sc%lab)
    sc%nl = 0 ; sc%nm = 0 ; sc%nq = 0 ; sc%nlab = 0 ; sc%nleg = 0
    ! image files (planet images) are kept across time steps
  end subroutine scene_reset

  !----------------------------------------------------------------------------
  ! read an image file (JPEG / PNG) and return its index (0 if it cannot be read)
  !----------------------------------------------------------------------------
  integer function add_image_file(sc, path) result(k)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: path
    integer :: u, ios
    integer(int64) :: n
    k = 0
    do u = 1, sc%nimf
      if (trim(sc%imf(u)%path) == trim(path)) then
        k = u ; return
      end if
    end do
    if (sc%nimf >= size(sc%imf)) return
    open(newunit=u, file=path, access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    inquire(unit=u, size=n)
    if (n < 8) then
      close(u) ; return
    end if
    sc%nimf = sc%nimf + 1 ; k = sc%nimf
    allocate(sc%imf(k)%d(n))
    read(u) sc%imf(k)%d
    close(u)
    sc%imf(k)%path = path
    if (uval(sc%imf(k)%d(1)) == 255 .and. uval(sc%imf(k)%d(2)) == 216) then
      sc%imf(k)%mime = 'image/jpeg'
    else if (uval(sc%imf(k)%d(1)) == 137 .and. uval(sc%imf(k)%d(2)) == 80) then
      sc%imf(k)%mime = 'image/png'
    else
      deallocate(sc%imf(k)%d) ; sc%nimf = sc%nimf - 1 ; k = 0
    end if
  end function add_image_file

  ! textured sphere (longitude = azimuth, north pole at the top of the image). lon0 [deg] rotates it
  subroutine add_textured_sphere(sc, name, c, r, nlat, nlon, itex, lon0)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    real, intent(in) :: c(3), r, lon0
    integer, intent(in) :: nlat, nlon, itex
    real, allocatable :: p(:,:), nn(:,:), col(:,:), uv(:,:)
    integer, allocatable :: tri(:,:)
    integer :: i, j, k, im, nv, nt
    real :: th, ph
    nv = (nlat+1)*(nlon+1) ; nt = 2*nlat*nlon
    allocate(p(3,nv), nn(3,nv), col(4,nv), tri(3,nt), uv(2,nv))
    k = 0
    do i = 0, nlat
      th = 3.14159265*real(i)/real(nlat)
      do j = 0, nlon
        ph = 2.0*3.14159265*real(j)/real(nlon)
        k = k + 1
        nn(:,k) = [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)]
        p(:,k) = c + r*nn(:,k)
        col(:,k) = 1.0
        uv(:,k) = [real(j)/real(nlon) + 0.5 - lon0/360.0, real(i)/real(nlat)]
      end do
    end do
    k = 0
    do i = 0, nlat-1
      do j = 0, nlon-1
        k = k + 1 ; tri(:,k) = [i*(nlon+1)+j+1, (i+1)*(nlon+1)+j+1, (i+1)*(nlon+1)+j+2]
        k = k + 1 ; tri(:,k) = [i*(nlon+1)+j+1, (i+1)*(nlon+1)+j+2, i*(nlon+1)+j+2]
      end do
    end do
    im = add_mesh_begin(sc, name, [1.0, 1.0, 1.0, 1.0], .true.)
    call add_triangles(sc, im, nv, p, nn, col, nt, tri)
    call mesh_uv(sc, im, uv, itex)
  end subroutine add_textured_sphere

  ! attach texture coordinates and an image to mesh im (as many uv as vertices)
  subroutine mesh_uv(sc, im, uv, itex)
    type(scene_t), intent(inout) :: sc
    integer, intent(in) :: im, itex
    real, intent(in) :: uv(:,:)
    sc%m(im)%uv = real(uv, real32)
    sc%m(im)%tex = itex
  end subroutine mesh_uv

  ! attach glTF node extras (JSON) to the mesh / image quad called name
  subroutine set_extras(sc, name, js)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name, js
    integer :: k
    do k = 1, sc%nq
      if (trim(sc%q(k)%name) == trim(name)) sc%q(k)%extras = js
    end do
    do k = 1, sc%nm
      if (trim(sc%m(k)%name) == trim(name)) sc%m(k)%extras = js
    end do
  end subroutine set_extras

  ! colour bar (rgb: 33 colours in 0..1)
  subroutine add_legend(sc, target, label, vmin, vmax, rgb)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: target, label
    real, intent(in) :: vmin, vmax, rgb(3,33)
    integer :: k
    do k = 1, sc%nleg
      if (trim(sc%leg(k)%target) == trim(target)) return
    end do
    if (sc%nleg >= size(sc%leg)) return
    sc%nleg = sc%nleg + 1
    sc%leg(sc%nleg) = legend_t(target, label, vmin, vmax, rgb)
  end subroutine add_legend

  !============================================================================
  ! building
  !============================================================================
  integer function line_group(sc, name)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    type(lines_t), allocatable :: t(:)
    integer :: k
    do k = 1, sc%nl
      if (sc%l(k)%name == name) then
        line_group = k ; return
      end if
    end do
    if (.not. allocated(sc%l)) allocate(sc%l(8))
    if (sc%nl == size(sc%l)) then
      allocate(t(2*sc%nl)) ; t(1:sc%nl) = sc%l(1:sc%nl) ; call move_alloc(t, sc%l)
    end if
    sc%nl = sc%nl + 1
    sc%l(sc%nl)%name = name
    allocate(sc%l(sc%nl)%p(3,4096), sc%l(sc%nl)%c(4,4096), sc%l(sc%nl)%e(2,4096))
    line_group = sc%nl
  end function line_group

  ! one polyline (colour per vertex rgba(4,n) in 0..1, or one colour)
  subroutine add_polyline(sc, name, n, p, rgba)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real, intent(in) :: p(3,n), rgba(:,:)
    integer :: g, k, nv0
    real(real32), allocatable :: tp(:,:)
    integer(int8), allocatable :: tc(:,:)
    integer(int32), allocatable :: te(:,:)
    if (n < 2) return
    g = line_group(sc, name)
    associate(L => sc%l(g))
      do while (L%nv + n > size(L%p,2))
        allocate(tp(3,2*size(L%p,2)), tc(4,2*size(L%p,2)))
        tp(:,1:L%nv) = L%p(:,1:L%nv) ; tc(:,1:L%nv) = L%c(:,1:L%nv)
        call move_alloc(tp, L%p) ; call move_alloc(tc, L%c)
      end do
      do while (L%ns + n - 1 > size(L%e,2))
        allocate(te(2,2*size(L%e,2))) ; te(:,1:L%ns) = L%e(:,1:L%ns)
        call move_alloc(te, L%e)
      end do
      nv0 = L%nv
      L%p(:, nv0+1:nv0+n) = p
      do k = 1, n
        if (size(rgba,2) >= n) then
          L%c(:, nv0+k) = u8(nint(255.0*rgba(:,k)))
        else
          L%c(:, nv0+k) = u8(nint(255.0*rgba(:,1)))
        end if
      end do
      do k = 1, n-1
        L%e(:, L%ns+k) = [nv0+k-1, nv0+k]
      end do
      L%nv = nv0 + n ; L%ns = L%ns + n - 1
    end associate
  end subroutine add_polyline

  ! start a new mesh (add triangles with add_triangles)
  integer function add_mesh_begin(sc, name, rgba, lit)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    real, intent(in) :: rgba(4)
    logical, intent(in) :: lit
    type(mesh_t), allocatable :: t(:)
    if (.not. allocated(sc%m)) allocate(sc%m(8))
    if (sc%nm == size(sc%m)) then
      allocate(t(2*sc%nm)) ; t(1:sc%nm) = sc%m(1:sc%nm) ; call move_alloc(t, sc%m)
    end if
    sc%nm = sc%nm + 1
    sc%m(sc%nm)%name = name ; sc%m(sc%nm)%rgba = rgba ; sc%m(sc%nm)%lit = lit
    allocate(sc%m(sc%nm)%p(3,0), sc%m(sc%nm)%n(3,0), sc%m(sc%nm)%c(4,0), sc%m(sc%nm)%t(3,0))
    add_mesh_begin = sc%nm
  end function add_mesh_begin

  ! vertices p(3,nv), normals nrm(3,nv), vertex colours col(4,nv) (0..1), triangles tri(3,nt) (1-based)
  subroutine add_triangles(sc, im, nv, p, nrm, col, nt, tri)
    type(scene_t), intent(inout) :: sc
    integer, intent(in) :: im, nv, nt
    real, intent(in) :: p(3,nv), nrm(3,nv), col(4,nv)
    integer, intent(in) :: tri(3,nt)
    real(real32), allocatable :: tp(:,:), tn(:,:)
    integer(int8), allocatable :: tc(:,:)
    integer(int32), allocatable :: tt(:,:)
    integer :: n0, t0
    associate(M => sc%m(im))
      n0 = M%nv ; t0 = M%nt
      allocate(tp(3,n0+nv), tn(3,n0+nv), tc(4,n0+nv), tt(3,t0+nt))
      tp(:,1:n0) = M%p(:,1:n0) ; tn(:,1:n0) = M%n(:,1:n0) ; tc(:,1:n0) = M%c(:,1:n0)
      tt(:,1:t0) = M%t(:,1:t0)
      tp(:,n0+1:) = p ; tn(:,n0+1:) = nrm ; tc(:,n0+1:) = u8(nint(255.0*col))
      tt(:,t0+1:) = tri - 1 + n0
      call move_alloc(tp, M%p) ; call move_alloc(tn, M%n)
      call move_alloc(tc, M%c) ; call move_alloc(tt, M%t)
      M%nv = n0 + nv ; M%nt = t0 + nt
    end associate
  end subroutine add_triangles

  subroutine add_image_quad(sc, name, corner, w, h, rgba)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    real, intent(in) :: corner(3,4)
    integer, intent(in) :: w, h
    real, intent(in) :: rgba(4,w,h)        ! row 1 = bottom (natural for the caller)
    type(quad_t), allocatable :: t(:)
    integer :: j
    if (.not. allocated(sc%q)) allocate(sc%q(8))
    if (sc%nq == size(sc%q)) then
      allocate(t(2*sc%nq)) ; t(1:sc%nq) = sc%q(1:sc%nq) ; call move_alloc(t, sc%q)
    end if
    sc%nq = sc%nq + 1
    associate(QD => sc%q(sc%nq))
      QD%name = name ; QD%corner = corner ; QD%w = w ; QD%h = h
      allocate(QD%img(4, w, h))
      do j = 1, h                              ! PNG starts from the top row
        QD%img(:,:,j) = u8(nint(255.0*rgba(:,:,h-j+1)))
      end do
    end associate
  end subroutine add_image_quad

  subroutine add_sphere(sc, name, c, r, rgba, nlat, nlon)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: name
    real, intent(in) :: c(3), r, rgba(4)
    integer, intent(in) :: nlat, nlon
    real, allocatable :: p(:,:), nn(:,:), col(:,:)
    integer, allocatable :: tri(:,:)
    integer :: i, j, k, im, nv, nt
    real :: th, ph
    nv = (nlat+1)*(nlon+1) ; nt = 2*nlat*nlon
    allocate(p(3,nv), nn(3,nv), col(4,nv), tri(3,nt))
    k = 0
    do i = 0, nlat
      th = 3.14159265*real(i)/real(nlat)
      do j = 0, nlon
        ph = 2.0*3.14159265*real(j)/real(nlon)
        k = k + 1
        nn(:,k) = [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)]
        p(:,k) = c + r*nn(:,k)
        col(:,k) = rgba
        ! northern hemisphere slightly brighter (to show the orientation)
        if (cos(th) > 0.0) col(1:3,k) = min(1.0, rgba(1:3)*1.3)
      end do
    end do
    k = 0
    do i = 0, nlat-1
      do j = 0, nlon-1
        k = k + 1 ; tri(:,k) = [i*(nlon+1)+j+1, (i+1)*(nlon+1)+j+1, (i+1)*(nlon+1)+j+2]
        k = k + 1 ; tri(:,k) = [i*(nlon+1)+j+1, (i+1)*(nlon+1)+j+2, i*(nlon+1)+j+2]
      end do
    end do
    im = add_mesh_begin(sc, name, [1.0, 1.0, 1.0, rgba(4)], .true.)
    call add_triangles(sc, im, nv, p, nn, col, nt, tri)
  end subroutine add_sphere

  subroutine add_label(sc, text, p, hsz, rgb)
    type(scene_t), intent(inout) :: sc
    character(len=*), intent(in) :: text
    real, intent(in) :: p(3), hsz, rgb(3)
    type(label_t), allocatable :: t(:)
    if (.not. allocated(sc%lab)) allocate(sc%lab(16))
    if (sc%nlab == size(sc%lab)) then
      allocate(t(2*sc%nlab)) ; t(1:sc%nlab) = sc%lab(1:sc%nlab) ; call move_alloc(t, sc%lab)
    end if
    sc%nlab = sc%nlab + 1
    sc%lab(sc%nlab) = label_t(text, p, hsz, rgb)
  end subroutine add_label

  ! display (X,Y,Z north) -> Y-up (X, Z, -Y)
  pure function yup(p) result(q)
    real, intent(in) :: p(3)
    real :: q(3)
    q = [p(1), p(3), -p(2)]
  end function yup

  !============================================================================
  ! PNG (stored deflate: written anywhere without libraries)
  !============================================================================
  subroutine png_encode(w, h, img, out)
    integer, intent(in) :: w, h
    integer(int8), intent(in) :: img(4, w, h)
    integer(int8), allocatable, intent(out) :: out(:)
    integer(int8), allocatable :: raw(:), z(:)
    integer(int64) :: nraw, pos, nz, k, blk
    integer :: j
    integer(int64) :: a1, a2
    nraw = int(h,int64)*(1 + 4*int(w,int64))
    allocate(raw(nraw))
    pos = 0
    do j = 1, h
      raw(pos+1) = 0_int8 ; pos = pos + 1          ! filter: none
      raw(pos+1:pos+4*w) = reshape(img(:,:,j), [4*w]) ; pos = pos + 4*w
    end do
    ! zlib stream with stored blocks
    nz = 2 + nraw + 5*((nraw + 65534)/65535 + 1) + 4
    allocate(z(nz))
    z(1) = int(z'78', int8) ; z(2) = 1_int8
    pos = 2 ; k = 0
    do
      blk = min(65535_int64, nraw - k)
      pos = pos + 1
      z(pos) = merge(1_int8, 0_int8, k + blk >= nraw)
      z(pos+1) = u8(int(iand(blk, 255_int64)))  ; z(pos+2) = u8(int(ishft(blk, -8)))
      z(pos+3) = u8(int(iand(not(blk), 255_int64))) ; z(pos+4) = u8(int(iand(ishft(not(blk), -8), 255_int64)))
      pos = pos + 4
      z(pos+1:pos+blk) = raw(k+1:k+blk) ; pos = pos + blk
      k = k + blk
      if (k >= nraw) exit
    end do
    a1 = 1 ; a2 = 0
    do k = 1, nraw
      a1 = mod(a1 + uval(raw(k)), 65521_int64)
      a2 = mod(a2 + a1, 65521_int64)
    end do
    z(pos+1:pos+4) = be32(ior(ishft(a2, 16), a1)) ; pos = pos + 4
    nz = pos
    ! PNG file
    allocate(out(8 + 25 + (12 + nz) + 12))
    out(1:8) = [u8(137), u8(80), u8(78), u8(71), u8(13), u8(10), u8(26), u8(10)]
    pos = 8
    call chunk('IHDR', [be32(int(w,int64)), be32(int(h,int64)), u8(8), u8(6), u8(0), u8(0), u8(0)])
    call chunk('IDAT', z(1:nz))
    call chunk('IEND', [integer(int8) ::])
  contains
    function be32(v) result(b)
      integer(int64), intent(in) :: v
      integer(int8) :: b(4)
      b = [u8(int(iand(ishft(v,-24),255_int64))), u8(int(iand(ishft(v,-16),255_int64))), &
           u8(int(iand(ishft(v,-8),255_int64))),  u8(int(iand(v,255_int64)))]
    end function be32
    subroutine chunk(typ, d)
      character(len=4), intent(in) :: typ
      integer(int8), intent(in) :: d(:)
      integer(int8) :: tb(4)
      integer :: i
      tb = [(int(iachar(typ(i:i)), int8), i=1,4)]
      out(pos+1:pos+4) = be32(int(size(d),int64)) ; pos = pos + 4
      out(pos+1:pos+4) = tb ; pos = pos + 4
      out(pos+1:pos+size(d)) = d ; pos = pos + size(d)
      out(pos+1:pos+4) = be32(crc32([tb, d])) ; pos = pos + 4
    end subroutine chunk
  end subroutine png_encode

  function crc32(d) result(c)
    integer(int8), intent(in) :: d(:)
    integer(int64) :: c
    integer(int64), save :: tab(0:255)
    logical, save :: init = .false.
    integer(int64) :: x
    integer :: n, k
    if (.not. init) then
      do n = 0, 255
        x = n
        do k = 1, 8
          if (iand(x, 1_int64) /= 0) then
            x = ieor(int(z'EDB88320', int64), ishft(x, -1))
          else
            x = ishft(x, -1)
          end if
        end do
        tab(n) = x
      end do
      init = .true.
    end if
    c = int(z'FFFFFFFF', int64)
    do n = 1, size(d)
      c = ieor(tab(iand(ieor(c, int(uval(d(n)),int64)), 255_int64)), ishft(c, -8))
    end do
    c = ieor(c, int(z'FFFFFFFF', int64))
  end function crc32

  !============================================================================
  ! glTF 2.0 binary
  !============================================================================
  subroutine bin_add(d)
    integer(int8), intent(in) :: d(:)
    integer(int8), allocatable :: t(:)
    integer(int64) :: need
    need = nbin + size(d) + 3
    if (need > size(bin)) then
      allocate(t(max(need, 2*size(bin,kind=int64))))
      t(1:nbin) = bin(1:nbin) ; call move_alloc(t, bin)
    end if
    bin(nbin+1:nbin+size(d)) = d ; nbin = nbin + size(d)
    do while (mod(nbin, 4_int64) /= 0)
      nbin = nbin + 1 ; bin(nbin) = 0_int8
    end do
  end subroutine bin_add

  ! add a bufferView and return its index
  integer function add_bv(d, target)
    integer(int8), intent(in) :: d(:)
    integer, intent(in) :: target
    integer(int64) :: off
    off = nbin
    call bin_add(d)
    if (nbv > 0) jbv = jbv//','
    if (target > 0) then
      jbv = jbv//'{"buffer":0,"byteOffset":'//i2s(off)//',"byteLength":'//i2s(int(size(d),int64))// &
            ',"target":'//i2s(int(target,int64))//'}'
    else
      jbv = jbv//'{"buffer":0,"byteOffset":'//i2s(off)//',"byteLength":'//i2s(int(size(d),int64))//'}'
    end if
    add_bv = nbv ; nbv = nbv + 1
  end function add_bv

  integer function add_acc(bv, ctype, count, typ, extra)
    integer, intent(in) :: bv, ctype, count
    character(len=*), intent(in) :: typ, extra
    if (nacc > 0) jacc = jacc//','
    jacc = jacc//'{"bufferView":'//i2s(int(bv,int64))//',"componentType":'//i2s(int(ctype,int64))// &
           ',"count":'//i2s(int(count,int64))//',"type":"'//typ//'"'//extra//'}'
    add_acc = nacc ; nacc = nacc + 1
  end function add_acc

  integer function add_mat(s)
    character(len=*), intent(in) :: s
    if (nmat > 0) jmat = jmat//','
    jmat = jmat//s
    add_mat = nmat ; nmat = nmat + 1
  end function add_mat

  subroutine add_node(name, prim, extras)
    character(len=*), intent(in) :: name, prim
    character(len=*), intent(in), optional :: extras
    if (nmeshes > 0) then
      jmesh = jmesh//',' ; jnode = jnode//','
    end if
    jmesh = jmesh//'{"name":"'//trim(name)//'","primitives":['//prim//']}'
    if (present(extras)) then
      if (len_trim(extras) > 0) then
        jnode = jnode//'{"name":"'//trim(name)//'","mesh":'//i2s(int(nmeshes,int64))// &
                ',"extras":'//trim(extras)//'}'
        nmeshes = nmeshes + 1
        return
      end if
    end if
    jnode = jnode//'{"name":"'//trim(name)//'","mesh":'//i2s(int(nmeshes,int64))//'}'
    nmeshes = nmeshes + 1
  end subroutine add_node

  function f2b(a) result(b)
    real(real32), intent(in) :: a(:,:)
    integer(int8), allocatable :: b(:)
    b = transfer(a, [0_int8])
  end function f2b

  function i2b(a) result(b)
    integer(int32), intent(in) :: a(:,:)
    integer(int8), allocatable :: b(:)
    b = transfer(a, [0_int8])
  end function i2b

  function minmax(p) result(s)
    real(real32), intent(in) :: p(:,:)
    character(len=:), allocatable :: s
    s = ',"min":['//r2s(minval(p(1,:)))//','//r2s(minval(p(2,:)))//','//r2s(minval(p(3,:)))// &
        '],"max":['//r2s(maxval(p(1,:)))//','//r2s(maxval(p(2,:)))//','//r2s(maxval(p(3,:)))//']'
  end function minmax

  function yup_all(p) result(q)
    real(real32), intent(in) :: p(:,:)
    real(real32) :: q(3, size(p,2))
    q(1,:) = p(1,:) ; q(2,:) = p(3,:) ; q(3,:) = -p(2,:)
  end function yup_all

  subroutine write_glb(sc, fname)
    type(scene_t), intent(in) :: sc
    character(len=*), intent(in) :: fname
    integer :: k, ap, ac, ai, an, at, mt, im, u, ios, bvi
    integer(int64) :: jl, bl
    character(len=:), allocatable :: js, ext, texs
    logical :: wrapuse
    real(real32), allocatable :: q(:,:)
    integer(int8), allocatable :: png(:)
    integer(int32) :: head(3), ch(2)
    integer(int32), parameter :: tri2(3,2) = reshape([0,1,2, 0,2,3], [3,2])
    real(real32), parameter :: uv(2,4) = reshape([0.0,1.0, 1.0,1.0, 1.0,0.0, 0.0,0.0], [2,4])

    allocate(bin(1048576)) ; nbin = 0
    jacc = '' ; jbv = '' ; jmesh = '' ; jnode = '' ; jmat = '' ; jimg = '' ; jtex = ''
    nacc = 0 ; nbv = 0 ; nmeshes = 0 ; nmat = 0 ; nimg = 0 ; wrapuse = .false. ; at = -1

    ! ---- line groups (one LINES primitive per group)
    do k = 1, sc%nl
      associate(L => sc%l(k))
        if (L%ns == 0) cycle
        q = yup_all(L%p(:,1:L%nv))
        ap = add_acc(add_bv(f2b(q), 34962), 5126, L%nv, 'VEC3', minmax(q))
        ac = add_acc(add_bv(reshape(L%c(:,1:L%nv),[4*L%nv]), 34962), 5121, L%nv, 'VEC4', ',"normalized":true')
        ai = add_acc(add_bv(i2b(L%e(:,1:L%ns)), 34963), 5125, 2*L%ns, 'SCALAR', '')
        mt = add_mat('{"name":"'//trim(L%name)//'","pbrMetallicRoughness":{"baseColorFactor":[1,1,1,1]},'// &
                     '"extensions":{"KHR_materials_unlit":{}}}')
        call add_node(L%name, '{"attributes":{"POSITION":'//i2s(int(ap,int64))//',"COLOR_0":'// &
             i2s(int(ac,int64))//'},"indices":'//i2s(int(ai,int64))//',"mode":1,"material":'//i2s(int(mt,int64))//'}')
      end associate
    end do
    ! ---- meshes
    do k = 1, sc%nm
      associate(M => sc%m(k))
        if (M%nt == 0) cycle
        q = yup_all(M%p)
        ap = add_acc(add_bv(f2b(q), 34962), 5126, M%nv, 'VEC3', minmax(q))
        q = yup_all(M%n)
        an = add_acc(add_bv(f2b(q), 34962), 5126, M%nv, 'VEC3', '')
        ac = add_acc(add_bv(reshape(M%c,[4*M%nv]), 34962), 5121, M%nv, 'VEC4', ',"normalized":true')
        ai = add_acc(add_bv(i2b(M%t), 34963), 5125, 3*M%nt, 'SCALAR', '')
        texs = ''
        if (M%tex > 0 .and. allocated(M%uv)) then
          at = add_acc(add_bv(f2b(M%uv), 34962), 5126, M%nv, 'VEC2', '')
          bvi = add_bv(sc%imf(M%tex)%d, 0)
          if (nimg > 0) then
            jimg = jimg//',' ; jtex = jtex//','
          end if
          jimg = jimg//'{"bufferView":'//i2s(int(bvi,int64))//',"mimeType":"'//trim(sc%imf(M%tex)%mime)//'"}'
          jtex = jtex//'{"source":'//i2s(int(nimg,int64))//',"sampler":1}'
          texs = ',"baseColorTexture":{"index":'//i2s(int(nimg,int64))//'}'
          nimg = nimg + 1 ; wrapuse = .true.
        end if
        ext = '"pbrMetallicRoughness":{"baseColorFactor":['//r2s(M%rgba(1))//','//r2s(M%rgba(2))//','// &
              r2s(M%rgba(3))//','//r2s(M%rgba(4))//']'//texs//',"metallicFactor":0,"roughnessFactor":0.8},'// &
              '"doubleSided":true'
        if (M%rgba(4) < 0.999 .or. any(M%c(4,:) /= -1_int8)) ext = ext//',"alphaMode":"BLEND"'
        if (.not. M%lit) ext = ext//',"extensions":{"KHR_materials_unlit":{}}'
        mt = add_mat('{"name":"'//trim(M%name)//'",'//ext//'}')
        if (len_trim(texs) > 0) then
          call add_node(M%name, '{"attributes":{"POSITION":'//i2s(int(ap,int64))//',"NORMAL":'// &
               i2s(int(an,int64))//',"TEXCOORD_0":'//i2s(int(at,int64))//'},"indices":'// &
               i2s(int(ai,int64))//',"mode":4,"material":'//i2s(int(mt,int64))//'}', M%extras)
        else
          call add_node(M%name, '{"attributes":{"POSITION":'//i2s(int(ap,int64))//',"NORMAL":'// &
               i2s(int(an,int64))//',"COLOR_0":'//i2s(int(ac,int64))//'},"indices":'//i2s(int(ai,int64))// &
               ',"mode":4,"material":'//i2s(int(mt,int64))//'}', M%extras)
        end if
      end associate
    end do
    ! ---- image quads
    do k = 1, sc%nq
      associate(QD => sc%q(k))
        q = yup_all(real(QD%corner, real32))
        ap = add_acc(add_bv(f2b(q), 34962), 5126, 4, 'VEC3', minmax(q))
        at = add_acc(add_bv(f2b(uv), 34962), 5126, 4, 'VEC2', '')
        ai = add_acc(add_bv(i2b(tri2), 34963), 5125, 6, 'SCALAR', '')
        call png_encode(QD%w, QD%h, QD%img, png)
        bvi = add_bv(png, 0)
        if (nimg > 0) then
          jimg = jimg//',' ; jtex = jtex//','
        end if
        jimg = jimg//'{"bufferView":'//i2s(int(bvi,int64))//',"mimeType":"image/png"}'
        jtex = jtex//'{"source":'//i2s(int(nimg,int64))//',"sampler":0}'
        im = nimg ; nimg = nimg + 1
        mt = add_mat('{"name":"'//trim(QD%name)//'","pbrMetallicRoughness":{"baseColorTexture":{"index":'// &
             i2s(int(im,int64))//'}},"alphaMode":"BLEND","doubleSided":true,'// &
             '"extensions":{"KHR_materials_unlit":{}}}')
        call add_node(QD%name, '{"attributes":{"POSITION":'//i2s(int(ap,int64))//',"TEXCOORD_0":'// &
             i2s(int(at,int64))//'},"indices":'//i2s(int(ai,int64))//',"mode":4,"material":'//i2s(int(mt,int64))// &
             '}', QD%extras)
      end associate
    end do

    ! ---- JSON
    js = '{"asset":{"version":"2.0","generator":"jupiter7 viztools mag3d"},'// &
         '"extensionsUsed":["KHR_materials_unlit"],"scene":0,'// &
         '"scenes":[{"name":"'//esc(sc%title)//'","nodes":['
    do k = 0, nmeshes-1
      if (k > 0) js = js//','
      js = js//i2s(int(k,int64))
    end do
    js = js//'],"extras":{"title":"'//esc(sc%title)//'","background":['//r2s(sc%bg(1))//','// &
         r2s(sc%bg(2))//','//r2s(sc%bg(3))//'],"labels":['
    do k = 1, sc%nlab
      if (k > 1) js = js//','
      q = yup_all(reshape(real(sc%lab(k)%p, real32), [3,1]))
      js = js//'{"text":"'//esc(sc%lab(k)%text)//'","p":['//r2s(q(1,1))//','//r2s(q(2,1))//','// &
           r2s(q(3,1))//'],"size":'//r2s(sc%lab(k)%size)//',"color":['//r2s(sc%lab(k)%rgb(1))//','// &
           r2s(sc%lab(k)%rgb(2))//','//r2s(sc%lab(k)%rgb(3))//']}'
    end do
    js = js//'],"legends":['
    do k = 1, sc%nleg
      if (k > 1) js = js//','
      js = js//'{"target":"'//esc(sc%leg(k)%target)//'","label":"'//esc(sc%leg(k)%label)// &
           '","vmin":'//r2s(sc%leg(k)%vmin)//',"vmax":'//r2s(sc%leg(k)%vmax)//',"colors":['
      do im = 1, 33
        if (im > 1) js = js//','
        js = js//'['//i2s(int(nint(255.0*sc%leg(k)%rgb(1,im)),int64))//','// &
             i2s(int(nint(255.0*sc%leg(k)%rgb(2,im)),int64))//','// &
             i2s(int(nint(255.0*sc%leg(k)%rgb(3,im)),int64))//']'
      end do
      js = js//']}'
    end do
    js = js//']}}],"nodes":['//jnode//'],"meshes":['//jmesh//'],"materials":['//jmat// &
         '],"accessors":['//jacc//'],"bufferViews":['//jbv//'],'// &
         '"buffers":[{"byteLength":'//i2s(nbin)//'}]'
    if (nimg > 0) js = js//',"images":['//jimg//'],"textures":['//jtex// &
         '],"samplers":[{"magFilter":9729,"minFilter":9729,"wrapS":33071,"wrapT":33071},'// &
         '{"magFilter":9729,"minFilter":9987,"wrapS":10497,"wrapT":33071}]'
    js = js//'}'
    do while (mod(len(js), 4) /= 0)
      js = js//' '
    end do
    jl = len(js) ; bl = nbin
    open(newunit=u, file=fname, access='stream', form='unformatted', status='replace', iostat=ios)
    if (ios /= 0) call fatal('cannot create '//fname)
    head = [int(z'46546C67', int32), 2_int32, int(12 + 8 + jl + 8 + bl, int32)]
    write(u) head
    ch = [int(jl, int32), int(z'4E4F534A', int32)]
    write(u) ch, js
    ch = [int(bl, int32), int(z'004E4942', int32)]
    write(u) ch, bin(1:nbin)
    close(u)
    deallocate(bin)
  end subroutine write_glb

  !============================================================================
  ! X3D (XML). Images are written to <base>_texN.png and referenced by ImageTexture
  !============================================================================
  subroutine write_x3d(sc, fname)
    type(scene_t), intent(in) :: sc
    character(len=*), intent(in) :: fname
    integer :: u, k, i, up
    real :: q(3)
    character(len=:), allocatable :: base, tex
    integer(int8), allocatable :: png(:)
    base = fname ; tex = ''
    if (index(base, '.', back=.true.) > index(base, '/', back=.true.)) base = base(1:index(base,'.',back=.true.)-1)
    scene_radius = 50.0
    do k = 1, sc%nl
      if (sc%l(k)%nv > 0) scene_radius = max(scene_radius, maxval(abs(sc%l(k)%p(:,1:sc%l(k)%nv))))
    end do
    do k = 1, sc%nq
      scene_radius = max(scene_radius, maxval(abs(sc%q(k)%corner)))
    end do
    open(newunit=u, file=fname, status='replace', form='formatted')
    write(u,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(u,'(a)') '<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.3//EN" "http://www.web3d.org/specifications/x3d-3.3.dtd">'
    write(u,'(a)') '<X3D profile="Immersive" version="3.3">'
    write(u,'(a)') '<head><meta name="generator" content="jupiter7 viztools mag3d"/>'
    write(u,'(a)') '<meta name="title" content="'//xesc(sc%title)//'"/></head>'
    write(u,'(a)') '<Scene>'
    write(u,'(a,3(1x,f5.3),a)') '<Background skyColor="', sc%bg, '"/>'
    write(u,'(a)') '<NavigationInfo type=''"EXAMINE" "ANY"''/>'
    call viewpoint(u, 'oblique', yup([0.55, 0.6, 0.45]))
    call viewpoint(u, 'from the Sun', yup([1.0, 0.0, 0.0]))
    call viewpoint(u, 'from dusk', yup([0.0, 1.0, 0.0]))
    call viewpoint(u, 'from north', yup([-0.01, 0.0, 1.0]))
    do k = 1, sc%nl
      associate(L => sc%l(k))
        if (L%ns == 0) cycle
        write(u,'(a)') '<Shape DEF="'//trim(L%name)//'"><IndexedLineSet colorPerVertex="true" coordIndex="'
        ! polylines: consecutive segments share vertices -> rebuild strips
        do i = 1, L%ns
          if (i == 1) then
            write(u,'(i0,1x,i0)', advance='no') L%e(1,i), L%e(2,i)
          else if (L%e(1,i) /= L%e(2,i-1)) then
            write(u,'(a,i0,1x,i0)', advance='no') ' -1 ', L%e(1,i), L%e(2,i)
          else
            write(u,'(1x,i0)', advance='no') L%e(2,i)
          end if
          if (mod(i, 40) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') ' -1">'
        write(u,'(a)', advance='no') '<Coordinate point="'
        do i = 1, L%nv
          q = yup(real(L%p(:,i)))
          write(u,'(3(f0.3,1x))', advance='no') q
          if (mod(i, 8) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') '"/>'
        write(u,'(a)', advance='no') '<ColorRGBA color="'
        do i = 1, L%nv
          write(u,'(4(f0.3,1x))', advance='no') real([(uval(L%c(up,i)), up=1,4)])/255.0
          if (mod(i, 8) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') '"/></IndexedLineSet></Shape>'
      end associate
    end do
    do k = 1, sc%nm
      associate(M => sc%m(k))
        if (M%nt == 0) cycle
        if (M%tex > 0 .and. allocated(M%uv)) then
          tex = base//'_img'//i2s(int(M%tex,int64))//merge('.jpg', '.png', sc%imf(M%tex)%mime == 'image/jpeg')
          call write_bytes(tex, sc%imf(M%tex)%d)
          write(u,'(a)') '<Shape DEF="'//trim(M%name)//'"><Appearance><Material diffuseColor="1 1 1"/>'// &
                         '<ImageTexture url=''"'//basename(tex)//'"'' repeatS="true" repeatT="false"/></Appearance>'
        else
          write(u,'(a)') '<Shape DEF="'//trim(M%name)//'"><Appearance><Material diffuseColor="1 1 1" transparency="'// &
                         r2s(1.0-M%rgba(4))//'"/></Appearance>'
        end if
        write(u,'(a)', advance='no') '<IndexedTriangleSet solid="false" colorPerVertex="true" normalPerVertex="true" index="'
        do i = 1, M%nt
          write(u,'(3(i0,1x))', advance='no') M%t(:,i)
          if (mod(i, 16) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') '">'
        write(u,'(a)', advance='no') '<Coordinate point="'
        do i = 1, M%nv
          write(u,'(3(f0.3,1x))', advance='no') yup(real(M%p(:,i)))
          if (mod(i, 8) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') '"/>'
        write(u,'(a)', advance='no') '<Normal vector="'
        do i = 1, M%nv
          write(u,'(3(f0.3,1x))', advance='no') yup(real(M%n(:,i)))
          if (mod(i, 8) == 0) write(u,'(a)') ''
        end do
        write(u,'(a)') '"/>'
        if (M%tex > 0 .and. allocated(M%uv)) then
          write(u,'(a)', advance='no') '<TextureCoordinate point="'
          do i = 1, M%nv
            write(u,'(2(f0.4,1x))', advance='no') M%uv(1,i), 1.0 - M%uv(2,i)
            if (mod(i, 8) == 0) write(u,'(a)') ''
          end do
          write(u,'(a)') '"/></IndexedTriangleSet></Shape>'
        else
          write(u,'(a)', advance='no') '<ColorRGBA color="'
          do i = 1, M%nv
            write(u,'(4(f0.3,1x))', advance='no') real([(uval(M%c(up,i)), up=1,4)])/255.0
            if (mod(i, 8) == 0) write(u,'(a)') ''
          end do
          write(u,'(a)') '"/></IndexedTriangleSet></Shape>'
        end if
      end associate
    end do
    do k = 1, sc%nq
      associate(QD => sc%q(k))
        tex = base//'_tex'//i2s(int(k,int64))//'.png'
        call png_encode(QD%w, QD%h, QD%img, png)
        call write_bytes(tex, png)
        write(u,'(a)') '<Shape DEF="'//trim(QD%name)//'"><Appearance><Material emissiveColor="1 1 1" diffuseColor="0 0 0"/>'
        write(u,'(a)') '<ImageTexture url=''"'//basename(tex)//'"''><TextureProperties boundaryModeS="CLAMP_TO_EDGE" '// &
                       'boundaryModeT="CLAMP_TO_EDGE"/></ImageTexture></Appearance>'
        write(u,'(a)') '<IndexedFaceSet solid="false" coordIndex="0 1 2 3 -1" texCoordIndex="0 1 2 3 -1">'
        write(u,'(a)', advance='no') '<Coordinate point="'
        do i = 1, 4
          write(u,'(3(f0.3,1x))', advance='no') yup(QD%corner(:,i))
        end do
        write(u,'(a)') '"/><TextureCoordinate point="0 0 1 0 1 1 0 1"/></IndexedFaceSet></Shape>'
      end associate
    end do
    do k = 1, sc%nlab
      write(u,'(a,3(f0.3,1x),a)') '<Transform translation="', yup(sc%lab(k)%p), '"><Billboard axisOfRotation="0 0 0">'
      write(u,'(a,3(f5.3,1x),a)') '<Shape><Appearance><Material emissiveColor="', sc%lab(k)%rgb, &
                                  '" diffuseColor="0 0 0"/></Appearance>'
      write(u,'(a)') '<Text string=''"'//xesc(sc%lab(k)%text)//'"''><FontStyle size="'//r2s(sc%lab(k)%size)// &
                     '" justify=''"MIDDLE" "MIDDLE"''/></Text></Shape></Billboard></Transform>'
    end do
    write(u,'(a)') '</Scene>'
    write(u,'(a)') '</X3D>'
    close(u)
  contains
    function basename(s) result(o)
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: o
      o = s(index(s,'/',back=.true.)+1:)
    end function basename
  end subroutine write_x3d

  ! Viewpoint looking at the origin (dir: view direction, Y-up). Distance from the scene size
  subroutine viewpoint(u, desc, dir)
    integer, intent(in) :: u
    character(len=*), intent(in) :: desc
    real, intent(in) :: dir(3)
    real :: bk(3), rt(3), up(3), r(3,3), ang, ax(3), d, tr, pos(3)
    d = 2.2*scene_radius
    bk = dir/sqrt(sum(dir**2))
    rt = [bk(3), 0.0, -bk(1)]                      ! cross((0,1,0), back)
    if (sum(rt**2) < 1.0e-8) rt = [1.0, 0.0, 0.0]
    rt = rt/sqrt(sum(rt**2))
    up = [bk(2)*rt(3)-bk(3)*rt(2), bk(3)*rt(1)-bk(1)*rt(3), bk(1)*rt(2)-bk(2)*rt(1)]
    r(:,1) = rt ; r(:,2) = up ; r(:,3) = bk
    tr = r(1,1) + r(2,2) + r(3,3)
    ang = acos(min(1.0, max(-1.0, 0.5*(tr - 1.0))))
    if (abs(sin(ang)) > 1.0e-4) then
      ax = [r(3,2)-r(2,3), r(1,3)-r(3,1), r(2,1)-r(1,2)]/(2.0*sin(ang))
    else if (ang < 1.0) then
      ax = [0.0, 1.0, 0.0] ; ang = 0.0
    else                                           ! 180 deg
      ax = sqrt(max(0.0, ([r(1,1), r(2,2), r(3,3)] + 1.0)/2.0))
      if (r(1,2) < 0.0) ax(2) = -ax(2)
      if (r(1,3) < 0.0) ax(3) = -ax(3)
    end if
    pos = d*bk
    write(u,'(a,3(f0.2,1x),a,4(f0.5,1x),a)') '<Viewpoint description="'//desc//'" position="', pos, &
         '" orientation="', ax, ang, '" centerOfRotation="0 0 0"/>'
  end subroutine viewpoint

  subroutine write_bytes(fname, d)
    character(len=*), intent(in) :: fname
    integer(int8), intent(in) :: d(:)
    integer :: u
    open(newunit=u, file=fname, access='stream', form='unformatted', status='replace')
    write(u) d
    close(u)
  end subroutine write_bytes

  !============================================================================
  function i2s(i) result(s)
    integer(int64), intent(in) :: i
    character(len=:), allocatable :: s
    character(len=24) :: b
    write(b,'(i0)') i
    s = trim(b)
  end function i2s

  function r2s(x) result(s)
    real, intent(in) :: x
    character(len=:), allocatable :: s
    character(len=24) :: b
    if (x == 0.0) then
      s = '0' ; return
    end if
    write(b,'(es14.7)') x
    s = trim(adjustl(b))
  end function r2s

  ! for JSON: avoid " and \
  function esc(s) result(o)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: o
    integer :: i
    o = ''
    do i = 1, len_trim(s)
      select case (s(i:i))
      case ('"')  ; o = o//"'"
      case ('\')  ; o = o//'/'
      case default; o = o//s(i:i)
      end select
    end do
  end function esc

  ! for XML attributes
  function xesc(s) result(o)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: o
    integer :: i
    o = ''
    do i = 1, len_trim(s)
      select case (s(i:i))
      case ('"')  ; o = o//'&quot;'
      case ("'")  ; o = o//'&apos;'
      case ('<')  ; o = o//'&lt;'
      case ('>')  ; o = o//'&gt;'
      case ('&')  ; o = o//'&amp;'
      case default; o = o//s(i:i)
      end select
    end do
  end function xesc

end module scene3d
