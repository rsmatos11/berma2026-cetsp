! =====================================================================
!  module blockcoordinatedescent: BCD for the visit points of a fixed
!  tour (BCD) and the single-disk subproblem (mindisk)
! =====================================================================

module blockcoordinatedescent

  use instance, only: instance_type,counter_type,evalf

  implicit none

  logical :: debug = .false.

  ! Choice of the point in Case 1 (segment [p,q] meets the disk, so every
  ! point of I = [p,q] cap B(C,R) is optimal): 'projection' (point of the
  ! segment closest to C), 'midpoint' (midpoint of I) or 'maximin' (point
  ! of I farthest from both p and q). 'maximin' keeps consecutive visit
  ! points apart.
  character(len=16) :: tie_rule = 'maximin'

  contains

  subroutine BCD(inst,count,tour,x,cost)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! Block coordinate descent: optimizes the visit points x of a fixed
    ! tour, one point at a time, in cyclic sweeps. On entry, x must be
    ! feasible and cost must be the length of (tour,x). Stops when a sweep
    ! reduces the length by a relative amount below sqrt(machine epsilon),
    ! or after count%maxbcdcic sweeps.

    ! PARAMETERS
    real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )

    ! LOCAL SCALARS
    integer :: k,kprev,knext,cic
    real(kind=8) :: gmin,oldloc,curcost,prevcost

    ! LOCAL ARRAYS
    real(kind=8) :: p(2),q(2)

    count%bcdmcalls = count%bcdmcalls + 1
    prevcost = cost
    curcost = cost

    cic = 0

    do

      cic = cic + 1

      ! Sweep over all tour positions (the anchor's visit point is also
      ! optimized; for a depot, a radius-zero disk, it does not move).
      do k = 1, inst%ntot
        if ( k .eq. 1 ) then
          kprev = inst%ntot           ! close the tour at the anchor
        else
          kprev = k - 1
        end if
        if ( k .eq. inst%ntot ) then
          knext = 1
        else
          knext = k + 1
        end if

        p(1:2) = x(1:2,kprev)
        q(1:2) = x(1:2,knext)

        ! Local cost of the two edges adjacent to x(1:2,k) before mindisk
        oldloc = norm2( x(1:2,k) - p(1:2) ) + norm2( x(1:2,k) - q(1:2) )

        call mindisk(inst,tour(k),p,q,x(1:2,k),gmin)

        ! Only the edges (kprev,k) and (k,knext) changed
        curcost = curcost + ( gmin - oldloc )
      end do

      ! Stopping criteria

      if ( curcost .ge. prevcost * ( 1.0d0 - macheps12 ) ) then
        count%stopbcdf = count%stopbcdf + 1
        count%bcdmcic = count%bcdmcic + cic
        cost = curcost
        return
      end if

      if ( cic .ge. count%maxbcdcic ) then
        count%stopbcdcic = count%stopbcdcic + 1
        count%bcdmcic = count%bcdmcic + cic
        cost = curcost
        return
      end if

      prevcost = curcost
    end do

  end subroutine BCD

  ! ------------------------------------------------------------

  subroutine mindisk(inst,pol,p,q,xmin,dmin)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst

    ! SCALAR ARGUMENTS
    integer, intent(in) :: pol
    real(kind=8), intent(out) :: dmin

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: p(2),q(2)
    real(kind=8), intent(inout) :: xmin(2)

    ! Single-disk subproblem: xmin = argmin_{x in B(C,R)} |x-p| + |x-q|.
    !   Case 1: segment [p,q] meets the disk; dmin = |p-q| and xmin is the
    !           point of the segment chosen by tie_rule.
    !   Case 2: the optimum lies on the boundary; found by 1D Newton on the
    !           angle theta of x(theta) = C + R (cos theta, sin theta).
    ! If R <= 0, xmin is the center.

    ! LOCAL SCALARS
    real(kind=8) :: Cx,Cy,R,a1,b1,a2,b2,dseg,theta_best,dcen

    if ( debug ) then
      ! Valid node range: 1..ntot (the depot included).
      if ( pol .lt. 1 .or. pol .gt. inst%ntot ) then
        write(*,*) 'There is something wrong!'
        write(*,*) 'pol = ',pol,' must be between 1 and ntot = ',inst%ntot
        stop
      end if
    end if

    Cx = inst%Cx(pol)
    Cy = inst%Cy(pol)
    R  = inst%Rad(pol)

    if ( R .le. 0.0d0 ) then
      xmin(1) = Cx
      xmin(2) = Cy
      dmin = norm2( xmin(1:2) - p(1:2) ) + norm2( xmin(1:2) - q(1:2) )
      return
    end if

    ! Shift to the disk-centered frame.
    a1 = p(1) - Cx;  b1 = p(2) - Cy
    a2 = q(1) - Cx;  b2 = q(2) - Cy

    dseg = seg_dist_origin(a1,b1,a2,b2)

    ! Tolerance on the Case 1 / Case 2 boundary (avoids sending a neighbor
    ! lying on the boundary to Newton, which would divide by zero).
    if ( dseg .le. R + 1.0d-12 * max( 1.0d0, hypot(a1,b1), hypot(a2,b2) ) ) then
      ! ---- Case 1: closed form ----
      select case ( trim(tie_rule) )
      case ('projection')
        call closest_point_on_segment(p(1),p(2),q(1),q(2),Cx,Cy,xmin(1),xmin(2))
      case ('midpoint')
        call segment_disk_midpoint(p(1),p(2),q(1),q(2),Cx,Cy,R,xmin(1),xmin(2))
      case ('maximin')
        call segment_disk_maximin(p(1),p(2),q(1),q(2),Cx,Cy,R,xmin(1),xmin(2))
      case default
        write(*,*) 'Unknown tie-break rule: ',trim(tie_rule),' - mindisk.'
        write(*,*) 'Valid options: projection, midpoint, maximin.'
        stop
      end select
      dmin = norm2( p(1:2) - q(1:2) )
      ! If the segment passes slightly outside the disk (within the
      ! tolerance), the point is pulled onto the boundary.
      dcen = hypot( xmin(1) - Cx, xmin(2) - Cy )
      if ( dcen .gt. R ) then
        xmin(1) = Cx + R * ( xmin(1) - Cx ) / dcen
        xmin(2) = Cy + R * ( xmin(2) - Cy ) / dcen
        dmin = norm2( xmin(1:2) - p(1:2) ) + norm2( xmin(1:2) - q(1:2) )
      end if
    else
      ! ---- Case 2: boundary optimum by 1D Newton ----
      call solve_newton(a1,b1,a2,b2,R,theta_best)
      xmin(1) = Cx + R * cos(theta_best)
      xmin(2) = Cy + R * sin(theta_best)
      dmin = norm2( xmin(1:2) - p(1:2) ) + norm2( xmin(1:2) - q(1:2) )
    end if

  end subroutine mindisk

  ! ------------------------------------------------------------
  ! Single-disk solver (used by mindisk)
  ! ------------------------------------------------------------

  real(kind=8) function seg_dist_origin(a1,b1,a2,b2)
    ! Distance from the origin to the segment [(a1,b1),(a2,b2)].
    implicit none
    real(kind=8), intent(in) :: a1,b1,a2,b2
    real(kind=8) :: v1,v2,vv,t,px,py
    v1 = a2 - a1
    v2 = b2 - b1
    vv = v1 * v1 + v2 * v2
    if ( vv .lt. 1.0d-14 ) then
      ! Degenerate segment: distance to the nearer endpoint.
      seg_dist_origin = min( hypot(a1,b1), hypot(a2,b2) )
      return
    else
      t = -(a1 * v1 + b1 * v2) / vv
      t = max( 0.0d0, min( 1.0d0, t ) )
    end if
    px = a1 + t * v1
    py = b1 + t * v2
    seg_dist_origin = hypot(px,py)
  end function seg_dist_origin

  subroutine closest_point_on_segment(P1x,P1y,P2x,P2y,Cx,Cy,projx,projy)
    ! Point of segment [P1,P2] closest to C.
    implicit none
    real(kind=8), intent(in)  :: P1x,P1y,P2x,P2y,Cx,Cy
    real(kind=8), intent(out) :: projx,projy
    real(kind=8) :: vx,vy,wx,wy,vv,t
    vx = P2x - P1x;  vy = P2y - P1y
    wx = Cx  - P1x;  wy = Cy  - P1y
    vv = vx * vx + vy * vy
    if ( vv .lt. 1.0d-14 ) then
      t = 0.0d0
    else
      t = (wx * vx + wy * vy) / vv
      t = max( 0.0d0, min( 1.0d0, t ) )
    end if
    projx = P1x + t * vx
    projy = P1y + t * vy
  end subroutine closest_point_on_segment

  real(kind=8) function f_theta(theta,a1,b1,a2,b2,R)
    ! f(theta) = |x(theta)-P1| + |x(theta)-P2|,
    ! x(theta) = C + R(cos theta, sin theta), disk-centered frame.
    implicit none
    real(kind=8), intent(in) :: theta,a1,b1,a2,b2,R
    real(kind=8) :: g1,g2
    g1 = hypot(R * cos(theta) - a1,R * sin(theta) - b1)
    g2 = hypot(R * cos(theta) - a2,R * sin(theta) - b2)
    f_theta = g1 + g2
  end function f_theta

  real(kind=8) function fprime_theta(theta,a1,b1,a2,b2,R)
    implicit none
    real(kind=8), intent(in) :: theta,a1,b1,a2,b2,R
    real(kind=8) :: g1,g2,N1,N2
    g1 = hypot(R * cos(theta) - a1,R * sin(theta) - b1)
    g2 = hypot(R * cos(theta) - a2,R * sin(theta) - b2)
    N1 = a1 * sin(theta) - b1 * cos(theta)
    N2 = a2 * sin(theta) - b2 * cos(theta)
    fprime_theta = R * ( N1 / g1 + N2 / g2 )
  end function fprime_theta

  real(kind=8) function fsecond_theta(theta,a1,b1,a2,b2,R)
    implicit none
    real(kind=8), intent(in) :: theta,a1,b1,a2,b2,R
    real(kind=8) :: g1,g2,N1,N2,N1p,N2p,g1p,g2p,g1pp,g2pp
    g1 = hypot(R * cos(theta) - a1,R * sin(theta) - b1)
    g2 = hypot(R * cos(theta) - a2,R * sin(theta) - b2)
    N1 = a1 * sin(theta) - b1 * cos(theta)
    N2 = a2 * sin(theta) - b2 * cos(theta)
    g1p = R * N1 / g1
    g2p = R * N2 / g2
    N1p = a1 * cos(theta) + b1 * sin(theta)
    N2p = a2 * cos(theta) + b2 * sin(theta)
    g1pp = ( R * N1p - g1p ** 2 ) / g1
    g2pp = ( R * N2p - g2p ** 2 ) / g2
    fsecond_theta = g1pp + g2pp
  end function fsecond_theta

  subroutine newton_one(theta0,a1,b1,a2,b2,R,theta_star)
    ! Newton's method in one variable (tol = 1e-13, at most 50 iterations).
    implicit none
    real(kind=8), intent(in)  :: theta0,a1,b1,a2,b2,R
    real(kind=8), intent(out) :: theta_star
    real(kind=8), parameter :: tol = 1.0d-13
    real(kind=8), parameter :: pi = 3.14159265358979323846d0
    real(kind=8) :: th,step
    integer :: k
    th = theta0
    do k = 1, 50
      step = fprime_theta(th,a1,b1,a2,b2,R) / fsecond_theta(th,a1,b1,a2,b2,R)
      th = th - step
      if ( abs(step) .lt. tol ) exit
    end do
    theta_star = modulo( th, 2.0d0 * pi )
  end subroutine newton_one

  subroutine solve_newton(a1,b1,a2,b2,R,theta_best)
    ! Adaptive multi-start Newton. Anchors: theta1 = direction of P1,
    ! theta2 = direction of P2, theta3 = direction of uhat1 + uhat2
    ! (bisector). Newton starts from whichever of theta1, theta2 has the
    ! smaller f; it also starts from the other one if the two values of f
    ! are close (relative gap < F_AMBIG_THRESHOLD), and from theta3 if f is
    ! nearly flat (rho = R |uhat1+uhat2| / (|P1|+|P2|) < RHO_THRESHOLD).
    ! The best result is kept. Here |P1|, |P2| > R > 0 (Case 2).
    implicit none
    real(kind=8), intent(in)  :: a1,b1,a2,b2,R
    real(kind=8), intent(out) :: theta_best
    real(kind=8), parameter :: RHO_THRESHOLD = 5.0d-2
    real(kind=8), parameter :: F_AMBIG_THRESHOLD = 1.0d-4
    real(kind=8) :: theta1,theta2,f1,f2,th,fval,fbest
    real(kind=8) :: nrm1,nrm2,sx,sy,rho

    theta1 = atan2(b1,a1)
    theta2 = atan2(b2,a2)
    f1 = f_theta(theta1,a1,b1,a2,b2,R)
    f2 = f_theta(theta2,a1,b1,a2,b2,R)

    ! primary anchor: smaller f
    if ( f1 .le. f2 ) then
      call newton_one(theta1,a1,b1,a2,b2,R,th)
    else
      call newton_one(theta2,a1,b1,a2,b2,R,th)
    end if
    fbest = f_theta(th,a1,b1,a2,b2,R)
    theta_best = th

    ! second anchor if the two values of f are close
    if ( abs(f1 - f2) .lt. F_AMBIG_THRESHOLD * min(f1,f2) ) then
      if ( f1 .le. f2 ) then
        call newton_one(theta2,a1,b1,a2,b2,R,th)
      else
        call newton_one(theta1,a1,b1,a2,b2,R,th)
      end if
      fval = f_theta(th,a1,b1,a2,b2,R)
      if ( fval .lt. fbest ) then
        fbest = fval
        theta_best = th
      end if
    end if

    ! bisector anchor if f is nearly flat
    nrm1 = hypot(a1,b1)
    nrm2 = hypot(a2,b2)
    sx = a1 / nrm1 + a2 / nrm2
    sy = b1 / nrm1 + b2 / nrm2
    rho = R * hypot(sx,sy) / (nrm1 + nrm2)
    if ( rho .lt. RHO_THRESHOLD ) then
      call newton_one(atan2(sy,sx),a1,b1,a2,b2,R,th)
      fval = f_theta(th,a1,b1,a2,b2,R)
      if ( fval .lt. fbest ) then
        fbest = fval
        theta_best = th
      end if
    end if
  end subroutine solve_newton

  ! ------------------------------------------------------------

  subroutine segment_disk_roots(P1x,P1y,P2x,P2y,Cx,Cy,R,tlo,thi,ok)
    ! Interval [tlo,thi] of t such that P1 + t (P2-P1), t in [0,1], lies in
    ! the disk. ok = .false. if the segment is degenerate or the interval is
    ! empty.
    implicit none
    real(kind=8), intent(in)  :: P1x,P1y,P2x,P2y,Cx,Cy,R
    real(kind=8), intent(out) :: tlo,thi
    logical, intent(out) :: ok
    real(kind=8) :: vx,vy,fx,fy,a,b,c,disc,sq,q,t1,t2
    real(kind=8) :: t0,footx,footy,dperp

    ok = .false.
    tlo = 0.0d0
    thi = 0.0d0

    vx = P2x - P1x;  vy = P2y - P1y
    fx = P1x - Cx;   fy = P1y - Cy
    a = vx * vx + vy * vy

    if ( a .lt. 1.0d-14 ) return

    b = 2.0d0 * (fx * vx + fy * vy)
    c = fx * fx + fy * fy - R * R

    ! Discriminant computed from the distance of C to the line (stable).
    t0 = -(fx * vx + fy * vy) / a
    footx = P1x + t0 * vx
    footy = P1y + t0 * vy
    dperp = hypot(footx - Cx,footy - Cy)
    disc = 4.0d0 * a * (R - dperp) * (R + dperp)

    if ( disc .lt. 0.0d0 ) return

    ! Stable quadratic formula.
    sq = sqrt(disc)
    if ( b .ge. 0.0d0 ) then
      q = -0.5d0 * (b + sq)
    else
      q = -0.5d0 * (b - sq)
    end if
    if ( abs(q) .gt. 1.0d-14 ) then
      t1 = q / a
      t2 = c / q
    else
      t1 = (-b - sq) / (2.0d0 * a)
      t2 = (-b + sq) / (2.0d0 * a)
    end if

    tlo = max( 0.0d0, min(t1,t2) )
    thi = min( 1.0d0, max(t1,t2) )

    if ( tlo .gt. thi ) return

    ok = .true.
  end subroutine segment_disk_roots

  subroutine segment_disk_midpoint(P1x,P1y,P2x,P2y,Cx,Cy,R,outx,outy)
    ! Midpoint of I = [P1,P2] cap B(C,R) (closest point if I is empty).
    implicit none
    real(kind=8), intent(in)  :: P1x,P1y,P2x,P2y,Cx,Cy,R
    real(kind=8), intent(out) :: outx,outy
    real(kind=8) :: tlo,thi,tmid
    logical :: ok

    call segment_disk_roots(P1x,P1y,P2x,P2y,Cx,Cy,R,tlo,thi,ok)

    if ( .not. ok ) then
      call closest_point_on_segment(P1x,P1y,P2x,P2y,Cx,Cy,outx,outy)
      return
    end if

    tmid = 0.5d0 * (tlo + thi)
    outx = P1x + tmid * (P2x - P1x)
    outy = P1y + tmid * (P2y - P1y)
  end subroutine segment_disk_midpoint

  subroutine segment_disk_maximin(P1x,P1y,P2x,P2y,Cx,Cy,R,outx,outy)
    ! Point of I maximizing min(|z-P1|,|z-P2|): the segment midpoint
    ! clamped to [tlo,thi].
    implicit none
    real(kind=8), intent(in)  :: P1x,P1y,P2x,P2y,Cx,Cy,R
    real(kind=8), intent(out) :: outx,outy
    real(kind=8) :: tlo,thi,tmid
    logical :: ok

    call segment_disk_roots(P1x,P1y,P2x,P2y,Cx,Cy,R,tlo,thi,ok)

    if ( .not. ok ) then
      call closest_point_on_segment(P1x,P1y,P2x,P2y,Cx,Cy,outx,outy)
      return
    end if

    tmid = max( tlo, min(0.5d0,thi) )
    outx = P1x + tmid * (P2x - P1x)
    outy = P1y + tmid * (P2y - P1y)
  end subroutine segment_disk_maximin

end module blockcoordinatedescent
