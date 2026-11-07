!==============================================================================
! psplot : PostScript plotting library (F90 version of the needed parts of the legacy plotting library)
!
!   correspondence with the legacy library
!     xyopen  -> ps_open        plots  -> ps_page     plote -> ps_endpage
!     plot    -> ps_plot        newpen -> ps_pen      newrgb -> ps_rgb
!     symbol  -> ps_text        symblc -> ps_text(font='sym')
!     imagc3  -> ps_image       grap3a -> ps_vectors   garohd -> ps_arrowhead
!     diff2   -> smooth2d
!   PostScript operators, units (idci=2: cm) and colour maps (ico=1..8) are the same as
!   the legacy library. Hex strings of images are assembled in a character buffer, which is
!   faster than the legacy imagc3 (formatted write 60z1).
!   Output goes to a file (given to ps_open) instead of the standard output.
!==============================================================================
module psplot
  use mhd_common, only: real32, real64, int64, fatal
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: ps_open, ps_close, ps_page, ps_endpage, ps_plot, ps_stroke, ps_pen, &
            ps_rgb, ps_text, ps_image, ps_vectors, ps_arrowhead, ps_box,       &
            ps_line, ps_colorbar, palette_rgb, smooth2d, ps_raw,         &
            ps_contour, ps_circle

  integer :: ups  = -1      ! output unit
  integer :: ldci = 2       ! 1: dot, 2: cm, 3: inch
  character(len=16), parameter :: hexd = '0123456789ABCDEF'

contains

  !----------------------------------------------------------------------------
  subroutine ps_open(fname)
    character(len=*), intent(in) :: fname
    integer :: ios
    character(len=256) :: msg
    open(newunit=ups, file=fname, form='formatted', status='replace', &
         action='write', iostat=ios, iomsg=msg)
    if (ios /= 0) call fatal('cannot create '//trim(fname)//' : '//trim(msg))
    write(ups,'(a)') '%!PS-Adobe-2.0'
    write(ups,'(a)') '%%Creator: jupiter7 viztools slice_plot'
    write(ups,'(a)') '/m { moveto } def'
    write(ups,'(a)') '/l { lineto } def'
    write(ups,'(a)') '/s { show } def'
    write(ups,'(a)') '/aws { awidthshow } def'
    write(ups,'(a)') '/sl { setlinewidth } def'
    write(ups,'(a)') '/sn { stroke newpath } def'
    write(ups,'(a)') '/tr { /Times-Roman findfont } def'
    write(ups,'(a)') '/sb { /Symbol findfont } def'
    write(ups,'(a)') '/sf { scalefont } def'
    write(ups,'(a)') '/se { setfont } def'
    write(ups,'(a)') '/ro { rotate } def'
    write(ups,'(a)') '/sd { setdash } def'
    write(ups,'(a)') '/st { stroke } def'
    write(ups,'(a)') '/sp { stroke showpage } def'
    write(ups,'(a)') '/tl { translate } def'
    write(ups,'(a)') '/sc { scale } def'
    write(ups,'(a)') '/gs { gsave } def'
    write(ups,'(a)') '/gr { grestore } def'
    write(ups,'(a)') '/grc { grestore 0.0 0.0 0.0 setrgbcolor } def'
    write(ups,'(a)') '/crg { setrgbcolor } def'
    write(ups,'(a)') '/f { fill } def'
    write(ups,'(a)') '/cup { currentpoint } def'
    ! right / center aligned text
    write(ups,'(a)') '/rshow { dup stringwidth pop neg 0 rmoveto show } def'
    write(ups,'(a)') '/cshow { dup stringwidth pop -2 div 0 rmoveto show } def'
  end subroutine ps_open

  subroutine ps_close()
    if (ups /= -1) close(ups)
    ups = -1
  end subroutine ps_close

  subroutine ps_raw(s)
    character(len=*), intent(in) :: s
    write(ups,'(a)') s
  end subroutine ps_raw

  !----------------------------------------------------------------------------
  ! legacy plots(idci,ity,isz,ino)
  !   idci=1 dot, 2 cm, 3 inch ; ity=1 portrait, 2 landscape ; isz=1 A4, 2 B4
  !----------------------------------------------------------------------------
  subroutine ps_page(idci, ity, isz, ino)
    integer, intent(in) :: idci, ity, isz, ino
    ldci = idci
    write(ups,'(a,i0,1x,i0)') '%%Page: ', ino, ino
    select case (idci)
    case (1)
      if (ity == 1) then
        write(ups,'(a)') ' 1.00 1.00 sc' ; write(ups,'(a)') ' 28.0 33.0 tl'
      else if (isz == 1) then
        write(ups,'(a)') ' 1.00 1.00 sc' ; write(ups,'(a)') ' 90.0 ro'
        write(ups,'(a)') ' 35.0 -565.0 tl'
      else
        write(ups,'(a)') ' 1.00 1.00 sc' ; write(ups,'(a)') ' 90.0 ro'
        write(ups,'(a)') ' 35.0 -696.0 tl'
      end if
    case (2)
      if (ity == 1) then
        write(ups,'(a)') ' 28.38  28.17 sc' ; write(ups,'(a)') '  1.00  1.20 tl'
      else if (isz == 1) then
        write(ups,'(a)') ' 90.0 ro' ; write(ups,'(a)') ' 28.17  28.38 sc'
        write(ups,'(a)') ' 1.15  -19.90 tl'
      else
        write(ups,'(a)') ' 90.0 ro' ; write(ups,'(a)') ' 28.17  28.38 sc'
        write(ups,'(a)') ' 1.15  -24.60 tl'
      end if
    case default
      if (ity == 1) then
        write(ups,'(a)') ' 72.08 71.55 sc' ; write(ups,'(a)') ' 0.39 0.47 tl'
      else if (isz == 1) then
        write(ups,'(a)') ' 90.0 ro' ; write(ups,'(a)') ' 71.55 72.08 sc'
        write(ups,'(a)') ' 0.49 -7.88 tl'
      else
        write(ups,'(a)') ' 90.0 ro' ; write(ups,'(a)') ' 71.55 72.08 sc'
        write(ups,'(a)') ' 0.49 -9.66 tl'
      end if
    end select
  end subroutine ps_page

  subroutine ps_endpage()
    write(ups,'(a)') 'sp'
  end subroutine ps_endpage

  subroutine ps_stroke()
    write(ups,'(a)') 'st'
  end subroutine ps_stroke

  !----------------------------------------------------------------------------
  ! legacy plot(x,y,im) : im=3 moveto, im=2 lineto
  !----------------------------------------------------------------------------
  subroutine ps_plot(x, y, im)
    real, intent(in) :: x, y
    integer, intent(in) :: im
    character(len=2) :: c
    c = merge(' m', ' l', im == 3)
    if (ldci == 1) then
      write(ups,'(f7.1,1x,f7.1,a)') x, y, c
    else
      write(ups,'(f7.3,1x,f7.3,a)') x, y, c
    end if
  end subroutine ps_plot

  subroutine ps_line(x1, y1, x2, y2)
    real, intent(in) :: x1, y1, x2, y2
    call ps_plot(x1, y1, 3)
    call ps_plot(x2, y2, 2)
  end subroutine ps_line

  subroutine ps_box(x0, y0, w, h)
    real, intent(in) :: x0, y0, w, h
    call ps_plot(x0,   y0,   3)
    call ps_plot(x0+w, y0,   2)
    call ps_plot(x0+w, y0+h, 2)
    call ps_plot(x0,   y0+h, 2)
    call ps_plot(x0,   y0,   2)
  end subroutine ps_box

  !----------------------------------------------------------------------------
  ! legacy newpen(ipn) (+ pen 21 of newpe1 = magenta)
  !   ipn = 10*(dash-1) + width  (1..10: solid, 11..: dashed)
  !----------------------------------------------------------------------------
  subroutine ps_pen(ipn)
    integer, intent(in) :: ipn
    real :: pn, fp
    integer :: ip, is
    real, parameter :: d2(10) = [1.00,0.05,0.10,0.20,0.30,0.40,0.20,0.40,0.60,0.80]
    real, parameter :: g2(10) = [0.00,0.05,0.10,0.20,0.30,0.40,0.10,0.10,0.10,0.10]
    pn = 0.012
    if (ldci == 1) pn = 0.350
    if (ldci == 3) pn = 0.005
    ip = max(mod(ipn-1,10)+1, 1)
    is = mod((ipn/10)-1,10)+1
    fp = real(ip)*pn
    write(ups,'(a,f6.3,a)') 'sn ', fp, ' sl'
    if (ipn > 10) then
      write(ups,'(a,f5.2,1x,f5.2,a)') '[', d2(is)*cfac(), g2(is)*cfac(), '] 0 sd'
    else
      write(ups,'(a)') '[] 0 sd'
    end if
    if (ipn == 21) call ps_rgb(1.0, 0.0, 1.0)
  contains
    real function cfac()
      cfac = 1.0
      if (ldci == 1) cfac = 30.0
      if (ldci == 3) cfac = 0.4
    end function cfac
  end subroutine ps_pen

  subroutine ps_rgb(r, g, b)
    real, intent(in) :: r, g, b
    write(ups,'(3(1x,f4.2),a)') min(max(r,0.0),1.0), min(max(g,0.0),1.0), &
                                min(max(b,0.0),1.0), ' crg'
  end subroutine ps_rgb

  !----------------------------------------------------------------------------
  ! legacy symbol / symblc
  !   align: 'l' (default), 'c', 'r' ; font: 'tr' Times-Roman, 'sym' Symbol
  !----------------------------------------------------------------------------
  subroutine ps_text(x, y, h, str, ang, align, font)
    real, intent(in) :: x, y, h
    character(len=*), intent(in) :: str
    real, intent(in), optional :: ang
    character(len=*), intent(in), optional :: align, font
    character(len=8) :: f
    character(len=:), allocatable :: s
    real :: a
    integer :: i
    a = 0.0 ; if (present(ang)) a = ang
    f = 'tr'  ; if (present(font)) f = font
    if (trim(f) == 'sym') f = 'sb'
    ! escape ( ) and backslash except "\ooo" octal codes
    s = ''
    do i = 1, len_trim(str)
      if (str(i:i) == '(' .or. str(i:i) == ')') then
        s = s//'\'//str(i:i)
      else
        s = s//str(i:i)
      end if
    end do
    write(ups,'(a)') 'gs'
    write(ups,'(a,1x,f6.3,a)') trim(f), h*1.5, ' sf se'
    write(ups,'(f7.3,1x,f7.3,a)') x, y, ' m'
    write(ups,'(f6.1,a)') a, ' ro'
    if (present(align)) then
      if (align == 'c') then
        write(ups,'(a)') '('//s//') cshow' ; write(ups,'(a)') 'gr' ; return
      else if (align == 'r') then
        write(ups,'(a)') '('//s//') rshow' ; write(ups,'(a)') 'gr' ; return
      end if
    end if
    write(ups,'(f6.3,a,f6.3,a)') h*0.06, ' 0.0 8#040 ', h*0.02, ' 0.0'
    write(ups,'(a)') '('//s//') aws'
    write(ups,'(a)') 'gr'
  end subroutine ps_text

  !----------------------------------------------------------------------------
  ! colour maps of the legacy imagc3 (x: 0..1)
  !----------------------------------------------------------------------------
  pure subroutine palette_rgb(x0, ico, r, g, b)
    real, intent(in) :: x0
    integer, intent(in) :: ico
    real, intent(out) :: r, g, b
    real :: x1, y1, a0
    x1 = x0 ; r = 0.0 ; g = 0.0 ; b = 0.0
    select case (ico)
    case (2)
      r = x1 ; g = x1 ; b = x1
    case (3)
      x1 = 7.0*x1
      b = x1
      if (x1 >= 2.0) b = 3.0-x1
      if (x1 >= 5.0) b = 0.5*(x1-5.0)
      g = x1-1.0
      if (x1 >= 4.0) g = 5.0-x1
      if (x1 >= 5.0) g = 0.5*(x1-5.0)
      r = x1-3.0
    case (4)
      x1 = 10.0*x1
      b = x1
      if (x1 >= 3.0) b = 4.0-x1
      if (x1 >= 6.0) b = 0.25*(x1-6.0)
      g = x1-2.0
      if (x1 >= 5.0) g = 6.0-x1
      if (x1 >= 6.0) g = 0.25*(x1-6.0)
      r = x1
      if (x1 >= 1.0) r = 2.0-x1
      if (x1 >= 4.0) r = x1-4.0
    case (5)
      a0 = 2.0
      y1 = x1*(a0+4.0)
      b = a0-y1
      if (y1 >= a0) b = y1-(a0+1.0)
      r = a0+3.0-y1
      g = a0+1.0-y1
      if (y1 >= a0+2.0) g = y1-(a0+3.0)
    case (6)
      a0 = 2.0
      y1 = x1*(a0+4.0)
      r = a0-y1
      if (y1 >= a0) r = y1-(a0+1.0)
      b = a0+3.0-y1
      g = a0+1.0-y1
      if (y1 >= a0+2.0) g = y1-(a0+3.0)
    case (8)                ! blue-white-red (diverging)
      if (x1 <= 0.495) then
        r = 4.0*(x1-0.245) ; b = 1.0
      end if
      if (x1 > 0.495 .and. x1 < 0.505) then
        r = 1.0 ; g = 1.0 ; b = 1.0
      end if
      if (x1 >= 0.505) then
        r = 1.0 ; b = -4.0*(x1-0.755)
      end if
      if (x1 <= 0.495 .and. x1 >= 0.245) then
        g = 4.0*(x1-0.245)
      else if (x1 <= 0.245) then
        g = -4.0*(x1-0.245)
      end if
      if (x1 >= 0.505 .and. x1 <= 0.755) then
        g = -4.0*(x1-0.755)
      else if (x1 >= 0.755) then
        g = 4.0*(x1-0.755)
      end if
    case default            ! ico = 1
      x1 = 2.0*x1
      r = x1 ; g = 2.0-x1 ; b = 0.0
    end select
    r = min(max(r,0.0),1.0) ; g = min(max(g,0.0),1.0) ; b = min(max(b,0.0),1.0)
  end subroutine palette_rgb

  !----------------------------------------------------------------------------
  ! legacy imagc3 : draw u(nx,ny) as a colour image in the area xl x yl [cm] at (xb,yb)
  !   written in iny strips (limit of PostScript string length, as before)
  !   4 bits per colour (as before)
  !----------------------------------------------------------------------------
  subroutine ps_image(nx, ny, xb, yb, xl, yl, u, vmin, vmax, ico, iny, bits)
    integer, intent(in) :: nx, ny, ico, iny
    real, intent(in) :: xb, yb, xl, yl, vmin, vmax
    real(real32), intent(in) :: u(nx, ny)
    integer, intent(in), optional :: bits        ! 4 (as before, default) or 8
    integer, parameter :: nc = 60
    integer(kind=2), allocatable :: cq(:,:,:)   ! (3, nx, ny) 0..15 or 0..255
    character(len=nc) :: line
    integer :: i, j, k, jj, iny1, iny2, ny1, ngy, nstrip, np, ic, nb, cmax
    real :: x, r, g, b, scale
    real :: yb1, yl1

    nb = 4
    if (present(bits)) nb = bits
    if (nb /= 8) nb = 4
    cmax = 2**nb - 1
    ! quantize (parallel)
    allocate(cq(3, nx, ny))
    scale = vmax - vmin
    if (scale == 0.0) scale = 1.0
    !$omp parallel do private(i,x,r,g,b)
    do j = 1, ny
      do i = 1, nx
        x = (u(i,j)-vmin)/scale
        if (ieee_is_nan(x)) then            ! NaN (outside the data) is white
          r = 1.0 ; g = 1.0 ; b = 1.0
        else
          call palette_rgb(x, ico, r, g, b)
        end if
        cq(1,i,j) = int(r*real(cmax), 2)
        cq(2,i,j) = int(g*real(cmax), 2)
        cq(3,i,j) = int(b*real(cmax), 2)
      end do
    end do
    !$omp end parallel do

    nstrip = max(1, min(iny, ny))
    iny1 = ny/nstrip
    iny2 = ny - iny1*(nstrip-1)
    do jj = 1, nstrip
      if (jj /= 1) write(ups,'(a)') 'gr'
      write(ups,'(a)') 'gs'
      ny1 = min(iny1 + 2, ny - iny1*(jj-1))      ! strips overlap by 2 rows (as before)
      if (jj == nstrip) ny1 = iny2
      ngy = iny1*(jj-1)
      yb1 = yb + yl*real(ngy)/real(ny)
      yl1 = yl*real(ny1)/real(ny)
      write(ups,'(f9.4,1x,f9.4,a)') xb, yb1, ' tl'
      write(ups,'(f9.5,1x,f9.5,a)') xl, yl1, ' sc'
      do k = 1, 3
        if (k == 1) write(ups,'(a)') '/Red <'
        if (k == 2) write(ups,'(a)') '/Green <'
        if (k == 3) write(ups,'(a)') '/Blue <'
        np = 0
        do j = ngy+1, ngy+ny1
          do i = 1, nx
            ic = cq(k,i,j)
            if (nb == 8) then
              call put(hexd(ic/16+1:ic/16+1)) ; call put(hexd(mod(ic,16)+1:mod(ic,16)+1))
            else
              call put(hexd(ic+1:ic+1))
            end if
          end do
          ! PostScript image rows start at byte boundaries: pad one digit for 4 bits and odd width
          if (nb == 4 .and. mod(nx,2) == 1) call put('0')
        end do
        if (np > 0) write(ups,'(a)') line(1:np)
        write(ups,'(a)') '> def'
      end do
      write(ups,'(i0,1x,i0,1x,i0)') nx, ny1, nb
      write(ups,'(a,i0,a,i0,a)') '[', nx, ' 0 0 ', ny1, ' 0 0]'
      write(ups,'(a)') '{Red} {Green} {Blue}'
      write(ups,'(a)') 'true 3 colorimage'
    end do
    write(ups,'(a)') 'grc'
    deallocate(cq)
  contains
    subroutine put(c)
      character(len=1), intent(in) :: c
      np = np + 1
      line(np:np) = c
      if (np == nc) then
        write(ups,'(a)') line ; np = 0
      end if
    end subroutine put
  end subroutine ps_image

  !----------------------------------------------------------------------------
  ! horizontal colour bar
  !----------------------------------------------------------------------------
  subroutine ps_colorbar(xb, yb, xl, yl, ico)
    real, intent(in) :: xb, yb, xl, yl
    integer, intent(in) :: ico
    integer, parameter :: n = 256
    real(real32) :: u(n,1)
    integer :: i
    do i = 1, n
      u(i,1) = real(i-1)/real(n-1)
    end do
    call ps_image(n, 1, xb, yb, xl, yl, u, 0.0, 1.0, ico, 1)
    call ps_pen(1)
    call ps_rgb(0.0, 0.0, 0.0)
    call ps_box(xb, yb, xl, yl)
    call ps_stroke()
  end subroutine ps_colorbar

  !----------------------------------------------------------------------------
  ! legacy garohd : arrow head at the end of the segment (x0,y0)->(x1,y1)
  !----------------------------------------------------------------------------
  subroutine ps_arrowhead(x0, y0, x1, y1, b, cr)
    real, intent(in) :: x0, y0, x1, y1, b, cr
    real :: c, aa, aa2, a1, a2, xh, yh, x2, x3, y2, y3
    real, parameter :: ep0 = 1.0e-10, ep1 = 1.0e-3, ep2 = 1.0e3
    c = cr*b
    if (abs(x1-x0) < ep0) then
      call vert()
    else if (abs(y1-y0) < ep0) then
      call horiz()
    else
      aa = (y1-y0)/(x1-x0)
      if (abs(aa) > ep2) then
        call vert()
      else if (abs(aa) < ep1) then
        call horiz()
      else
        aa2 = aa*aa
        a1 = b/sqrt(1.0+aa2)
        a2 = c/sqrt(1.0+1.0/aa2)
        xh = x1-a1
        if (x1 < x0) xh = x1+a1
        x2 = xh+a2
        x3 = xh-a2
        yh = y1+aa*(xh-x1)
        y2 = yh-(x2-xh)/aa
        y3 = yh-(x3-xh)/aa
      end if
    end if
    call ps_plot(x0, y0, 3)
    call ps_plot(x1, y1, 2)
    call ps_plot(x2, y2, 3)
    call ps_plot(x1, y1, 2)
    call ps_plot(x3, y3, 2)
  contains
    subroutine vert()
      x2 = x1-c ; x3 = x1+c
      y2 = y1-b
      if (y1 < y0) y2 = y1+b
      y3 = y2
    end subroutine vert
    subroutine horiz()
      y2 = y1-c ; y3 = y1+c
      x2 = x1-b
      if (x1 < x0) x2 = x1+b
      x3 = x2
    end subroutine horiz
  end subroutine ps_arrowhead

  !----------------------------------------------------------------------------
  ! legacy grap3a : draw a 2-component vector field (ux,uy) as short streamline arrows
  !   n1 x n2 grid (points 1 and n1 at the left and right edges) mapped to xb..xb+xl, yb..yb+yl
  !   traced las1 steps from nsx x nsy seed points, with an arrow head at the end
  !   umax : normalization (legacy vmax). urmin : minimum drawn (relative to umax)
  !   iar=1 : arrow head size proportional to the mean strength (as before)
  !   (fixes the coefficient (1-Z)*R -> (1-R)*Z of the uy interpolation in the legacy grap3a)
  !----------------------------------------------------------------------------
  subroutine ps_vectors(n1, n2, ux, uy, nsx, nsy, xb, yb, xl, yl, umax0, &
                        urmin, las1, iar, ipen)
    integer, intent(in) :: n1, n2, nsx, nsy, las1, iar, ipen
    real(real32), intent(in) :: ux(n1,n2), uy(n1,n2)
    real, intent(in) :: xb, yb, xl, yl, umax0, urmin
    real :: hx, hy, dx, dy, ad, ar, cr, umax, urmi
    real :: r1, z1, x1, y1, x2, y2, rr, r, z, vx, vy, vm, r2, z2, xx1
    integer :: i, j, iii, ir1, iz1
    logical :: drawn

    umax = umax0
    if (umax <= 0.0) umax = 1.0
    urmi = urmin*umax
    hx = real(n1-1)/real(nsx+1)
    hy = real(n2-1)/real(nsy+1)
    dx = xl/real(n1-1)
    dy = yl/real(n2-1)
    ad = 0.1*hx
    ar = dx*hx*0.60
    cr = 0.40

    do j = 1, nsy
      call ps_pen(ipen)
      do i = 1, nsx
        r1 = real(i)*hx
        z1 = real(j)*hy
        x2 = xb + dx*r1
        y2 = yb + dy*z1
        x1 = x2 ; y1 = y2
        rr = 0.0
        drawn = .false.
        do iii = 1, las1
          x1 = x2 ; y1 = y2
          ir1 = int(r1) + 1
          iz1 = int(z1) + 1
          if (ir1 < 1 .or. ir1 > n1-1) cycle
          if (iz1 < 1 .or. iz1 > n2-1) cycle
          r = r1 - int(r1)
          z = z1 - int(z1)
          vx = r*z*ux(ir1+1,iz1+1) + (1.0-r)*z*ux(ir1,iz1+1) &
             + r*(1.0-z)*ux(ir1+1,iz1) + (1.0-r)*(1.0-z)*ux(ir1,iz1)
          vy = r*z*uy(ir1+1,iz1+1) + (1.0-r)*z*uy(ir1,iz1+1) &
             + r*(1.0-z)*uy(ir1+1,iz1) + (1.0-r)*(1.0-z)*uy(ir1,iz1)
          vm = sqrt(vx*vx + vy*vy)
          if (vm < 0.01*urmi) then      ! legacy: GO TO 12 (no arrow head for this seed)
            drawn = .false. ; exit
          end if
          if (vm > umax) then
            vx = vx*umax/vm ; vy = vy*umax/vm
          end if
          vx = ad*vx/umax
          vy = ad*vy/umax
          rr = rr + min(vm/umax, 1.0)
          r2 = r1 + vx
          z2 = z1 + vy
          x2 = xb + dx*r2
          y2 = yb + dy*z2
          call ps_plot(x1, y1, 3)
          call ps_plot(x2, y2, 2)
          drawn = .true.
          r1 = r2 ; z1 = z2
          if (vm < urmi) exit
          if (r1 <= 0.0 .or. r1 >= real(n1-1)) exit
          if (z1 <= 0.0 .or. z1 >= real(n2-1)) exit
        end do
        if (.not. drawn) cycle
        xx1 = rr/real(las1)
        if (iar == 1) then
          call ps_arrowhead(x1, y1, x2, y2, ar*xx1*1.5,      cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*xx1*1.5*0.35, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*xx1*1.5*0.55, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*xx1*1.5*0.65, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*xx1*1.5*0.85, cr)
        else
          call ps_arrowhead(x1, y1, x2, y2, ar,      cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*0.35, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*0.55, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*0.65, cr)
          call ps_arrowhead(x1, y1, x2, y2, ar*0.85, cr)
        end if
      end do
      call ps_stroke()
    end do
    call ps_pen(1)
    call ps_rgb(0.0, 0.0, 0.0)
  end subroutine ps_vectors

  !----------------------------------------------------------------------------
  ! legacy diff2 : 5-point smoothing (coefficient 1/8) m times; edges mirrored.
  !   (the legacy behaviour of skipping the middle row was removed)
  !----------------------------------------------------------------------------
  subroutine smooth2d(m, nx, ny, u)
    integer, intent(in) :: m, nx, ny
    real(real32), intent(inout) :: u(nx, ny)
    real(real32), allocatable :: v(:,:)
    real(real32), parameter :: alp = 1.0/8.0, c0 = 1.0 - 4.0*alp
    integer :: ii, i, j
    if (m <= 0 .or. nx < 3 .or. ny < 3) return
    allocate(v(nx, ny))
    do ii = 1, m
      !$omp parallel do
      do j = 1, ny
        v(:,j) = u(:,j)
      end do
      !$omp end parallel do
      !$omp parallel do private(i)
      do j = 2, ny-1
        do i = 2, nx-1
          u(i,j) = c0*v(i,j) + alp*(v(i-1,j)+v(i+1,j)+v(i,j-1)+v(i,j+1))
        end do
        u(1,j)  = c0*v(1,j)  + alp*(2.0*v(2,j)    + v(1,j-1)  + v(1,j+1))
        u(nx,j) = c0*v(nx,j) + alp*(2.0*v(nx-1,j) + v(nx,j-1) + v(nx,j+1))
      end do
      !$omp end parallel do
      do i = 2, nx-1
        u(i,1)  = c0*v(i,1)  + alp*(v(i-1,1)  + v(i+1,1)  + 2.0*v(i,2))
        u(i,ny) = c0*v(i,ny) + alp*(v(i-1,ny) + v(i+1,ny) + 2.0*v(i,ny-1))
      end do
    end do
    deallocate(v)
  end subroutine smooth2d

  !----------------------------------------------------------------------------
  ! contour lines (marching squares). u(nx,ny) mapped to xb..xb+xl, yb..yb+yl;
  ! the contour of value level is drawn with the current pen. Cells with NaN are skipped.
  !----------------------------------------------------------------------------
  subroutine ps_contour(nx, ny, u, xb, yb, xl, yl, level)
    integer, intent(in) :: nx, ny
    real(real32), intent(in) :: u(nx, ny)
    real, intent(in) :: xb, yb, xl, yl, level
    integer :: i, j, c, np
    real :: v(4), px(4), py(4), sx(4), sy(4), dx, dy, t
    integer, parameter :: e1(4) = [1, 2, 3, 4], e2(4) = [2, 3, 4, 1]
    dx = xl/real(nx-1) ; dy = yl/real(ny-1)
    do j = 1, ny-1
      do i = 1, nx-1
        v = [u(i,j), u(i+1,j), u(i+1,j+1), u(i,j+1)]
        if (any(ieee_is_nan(v))) cycle
        if (all(v >= level) .or. all(v < level)) cycle
        px = [xb+dx*(i-1), xb+dx*i, xb+dx*i,     xb+dx*(i-1)]
        py = [yb+dy*(j-1), yb+dy*(j-1), yb+dy*j, yb+dy*j]
        np = 0
        do c = 1, 4
          if ((v(e1(c)) >= level) .neqv. (v(e2(c)) >= level)) then
            t = (level - v(e1(c)))/(v(e2(c)) - v(e1(c)))
            np = np + 1
            sx(np) = px(e1(c)) + t*(px(e2(c))-px(e1(c)))
            sy(np) = py(e1(c)) + t*(py(e2(c))-py(e1(c)))
          end if
        end do
        if (np >= 2) call ps_line(sx(1), sy(1), sx(2), sy(2))
        if (np == 4) call ps_line(sx(3), sy(3), sx(4), sy(4))
      end do
    end do
    call ps_stroke()
  end subroutine ps_contour

  subroutine ps_circle(x, y, r)
    real, intent(in) :: x, y, r
    write(ups,'(f8.3,1x,f8.3,a)') x+r, y, ' m'
    write(ups,'(3(1x,f8.3),a)') x, y, r, ' 0 360 arc'
    write(ups,'(a)') 'st'
  end subroutine ps_circle

end module psplot
