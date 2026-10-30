! =====================================================================
!  module instance: data types, input (inputdata), MetaPost figure
!  (drawsol), tour length (evalf) and random numbers (drand)
! =====================================================================

module instance

  implicit none

  type :: instance_type
    integer :: nballs,ntot,idepot,ianchor,vsize,lwsize,lenfilename,ncand
    ! nballs  : number of disks in the instance file
    ! ntot    : number of nodes in the tour (nballs, plus the depot if any)
    ! idepot  : index of the depot (0 if there is none)
    ! ianchor : node fixed at tour position 1 (the depot, or disk 1); its
    !           position is fixed, its visit point is optimized by BCD
    logical :: hasdepot
    real(kind=8) :: xymax
    character(len=1024) :: filename
    character(len=256)  :: basename
    ! basename: file name without directory and extension (used to name
    ! the output files)
    integer, allocatable :: cand(:,:) ! Arrays holding the candidate (neighbor) lists
    real(kind=8), allocatable :: Cx(:),Cy(:),Rad(:),canddist(:,:)
  end type instance_type

  type :: parameter_type
    integer :: jmin,jmax,ell_np_max,kappamax,ncool_steps,strlen,reverse,cardZ,&
    numdbmoves,bsk,vnd_type,seedrun,ncand_max,bestfirst,ncool_probe,maxbcdcic,refine
    ! ncool_steps : number of threshold reductions (kappa in the paper)
    ! strlen      : maximum string length (L); reverse: 1 = reversed strings
    ! cardZ       : positions reassigned by the Gutin neighborhood (|Z|)
    ! numdbmoves  : double-bridge moves per perturbation (ell)
    ! bsk         : Balas-Simonetti parameter (k), 2 or 3
    ! vnd_type    : 0 = fixed neighborhood order; 1 = random order
    ! bestfirst   : 0 = first improvement; 1 = best improvement (2-opt, swap,
    !               insertion and 3-opt)
    ! ncand_max   : size of the candidate lists (K)
    ! maxbcdcic   : maximum BCD sweeps per evaluation in VND and after each
    !               perturbation (s_max in the paper)
    ! refine      : 1 = complete BCD (up to 200 sweeps) after each VND; 0 = no
    logical :: cool_marginal
    real(kind=8) :: seed,eps,h,cool_probe_frac
    ! usecplex: .false. = visit points computed by BCD; .true. = by CPLEX.
    ! Read from the environment variable ILSUSECPLEX.
    logical :: usecplex
  end type parameter_type

  type :: counter_type
    integer :: bcdmcalls,bcdmcic,stopbcdf,stopbcdcic,vndcal,nscal,&
    hcalls,ilsiter,twoptcall,treeoptcall,swapcall,insercall,dbcall,&
    relocall,bscall,gutincall,maxbcdcic,incnimp,cplexcalls
    ! incnimp : how many improvements have been reported so far
    real(kind=8) :: inccost,inct0
    ! inccost : best cost found so far, in NORMALIZED units
  end type counter_type

  contains

  subroutine inputdata(inst,count,param)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(inout) :: inst
    type(parameter_type), intent(inout) :: param
    type(counter_type), intent(inout) :: count

    ! LOCAL SCALARS
    integer :: i,j,k,m,iworst,allocerr,depotflag,ibeg,iend
    real(kind=8) :: dum,worst,dmin
    character(len=256) :: dumcar
    
    ! LOCAL ARRAYS
    real(kind=8), allocatable :: rowd(:)

    ! ---- parameters read from standard input ----
    read(*,*) inst%filename ! instance file path
    if ( len_trim(inst%filename) .ge. len(inst%filename) ) then
      write(*,*) 'Instance file path too long (maximum ',len(inst%filename)-1, &
      ' characters) - inputdata.'
      stop 1
    end if
    read(*,*) depotflag ! 1 = the instance file carries a depot line, 0 = it does not
    read(*,*) param%seed, param%seedrun ! seed and repetition index
    read(*,*) param%ncool_steps, param%eps, param%h, param%strlen, param%reverse, &
    param%cardZ, param%numdbmoves, param%bsk, param%vnd_type, param%bestfirst, &
    param%ncand_max, param%maxbcdcic, param%refine

    ! Validation of the parameters
    if ( param%bsk .ne. 2 .and. param%bsk .ne. 3 ) then
      write(*,*) 'bsk = ',param%bsk,' - only 2 and 3 are implemented in ', &
      'BalaSimonetti - inputdata.'
      stop 1
    end if

    ! 0 = fixed neighborhood order, 1 = order drawn at random per VND call.
    if ( param%vnd_type .ne. 0 .and. param%vnd_type .ne. 1 ) then
      write(*,*) 'vnd_type = ',param%vnd_type,' - use 0 for the fixed ', &
      'neighborhood order or 1 for the random one.'
      stop 1
    end if

    if ( param%bestfirst .ne. 0 .and. param%bestfirst .ne. 1 ) then
      write(*,*) 'bestfirst = ',param%bestfirst,' - it is a flag: 0 for ', &
      'first improvement, 1 for best improvement - inputdata.'
      stop 1
    end if

    if ( param%reverse .ne. 0 .and. param%reverse .ne. 1 ) then
      write(*,*) 'reverse = ',param%reverse,' - it is a flag, use 0 or 1; ', &
      'StringRelocN tests it against 1, so any other value silently ', &
      'means 0 - inputdata.'
      stop 1
    end if

    if ( param%cardZ .lt. 0 .or. param%strlen .lt. 0 .or. &
         param%numdbmoves .lt. 0 ) then
      write(*,*) 'cardZ, strlen and numdbmoves must all be non-negative', &
      ' - inputdata.'
      stop 1
    end if

    ! ncand_max is the size K of the candidate lists (capped at ntot-1).
    if ( param%ncand_max .lt. 1 ) then
      write(*,*) 'ncand_max = ',param%ncand_max,' must be at least 1 ', &
      '(it is K, the size of the candidate lists) - inputdata.'
      stop 1
    end if

    if ( param%maxbcdcic .lt. 1 ) then
      write(*,*) 'maxbcdcic = ',param%maxbcdcic,' must be at least 1 - inputdata.'
      stop 1
    end if

    if ( param%refine .ne. 0 .and. param%refine .ne. 1 ) then
      write(*,*) 'refine = ',param%refine,' - it is a flag, use 0 or 1 - inputdata.'
      stop 1
    end if

    if ( param%ncool_steps .lt. 1 ) then
      write(*,*) 'ncool_steps = ',param%ncool_steps,' must be at least 1 ', &
      '(it is the number of coolings of eps over the whole run) - inputdata.'
      stop 1
    end if

    if ( param%eps .lt. 0.0d0 ) then
      write(*,*) 'eps = ',param%eps,' must be non-negative - inputdata.'
      stop 1
    end if

    ! h must be in (0,1].
    if ( param%h .le. 0.0d0 .or. param%h .gt. 1.0d0 ) then
      write(*,*) 'h = ',param%h,' must lie in (0,1] - inputdata.'
      stop 1
    end if

    ! Stagnation limit: the ILS stops after ell_np_max consecutive
    ! iterations without improvement of the best solution (applied only
    ! when ntot <= 200; larger instances stop by the time limit).
    param%ell_np_max = 5000 ! stagnation limit (used only if ntot <= 200)

    ! Time-based schedule of the threshold reductions: during a short
    ! probe (ncool_probe iterations or a fraction cool_probe_frac of the
    ! time limit, whichever comes first) the time per iteration is
    ! measured; the reductions are then spread over the iterations that
    ! still fit in the time limit.
    param%ncool_probe     = 50
    param%cool_probe_frac = 0.05d0
    param%cool_marginal   = .true.

    ! Reduction period used when no time limit is given (ilstlim unset).
    param%kappamax = max( 1, param%ell_np_max / param%ncool_steps )
    param%jmin       = 1 ! first VND neighborhood index
    param%jmax       = 8 ! last VND neighborhood index

    ! Engine for the visit points (environment variable ILSUSECPLEX):
    ! unset or not starting with 1/T/t -> BCD.
    param%usecplex = .false.
    call get_environment_variable('ILSUSECPLEX',dumcar)
    if ( len_trim(dumcar) .gt. 0 ) then
      if ( dumcar(1:1) .eq. '1' .or. dumcar(1:1) .eq. 'T' .or. dumcar(1:1) .eq. 't' ) then
        param%usecplex = .true.
      end if
    end if

    inst%filename    = adjustl(inst%filename)
    inst%lenfilename = len_trim(inst%filename)

    ! Strip the directory path and the extension.
    ibeg = index( inst%filename(1:inst%lenfilename),'/',back=.true. ) + 1
    iend = index( inst%filename(1:inst%lenfilename),'.',back=.true. ) - 1
    if ( iend .lt. ibeg ) iend = inst%lenfilename
    inst%basename = inst%filename(ibeg:iend)
    count%maxbcdcic  = 200

    param%seed = param%seed + 1000.0d0 * dble(param%seedrun)
    do k = 1, param%seedrun
      dum = drand(param%seed)
    end do

    ! ---- instance file ----
    ! With a depot (depotflag /= 0):
    !   line 1: number of disks; line 2: depot "x y"; then "x y r" per disk.
    !   The depot is stored as the last node (radius 0) and is the anchor.
    ! Without a depot (depotflag == 0):
    !   line 1: number of disks; then "x y r" per disk. Disk 1 is the anchor.
    inst%hasdepot = ( depotflag .ne. 0 )

    open(unit=10,file=trim(adjustl(inst%filename)),status='old')
    read(10,*) inst%nballs

    if ( inst%hasdepot ) then
      inst%ntot    = inst%nballs + 1     ! total number of nodes (customers + depot)
      inst%idepot  = inst%nballs + 1     ! depot stored as the last node, radius 0
      inst%ianchor = inst%idepot
    else
      inst%ntot    = inst%nballs
      inst%idepot  = 0                   ! no depot in this instance
      inst%ianchor = 1
    end if

    if ( inst%ntot .lt. 3 ) then
      write(*,*) 'Instance has ',inst%ntot,' nodes; at least 3 are needed', &
      ' to form a tour - inputdata.'
      stop
    end if

    inst%ncand = min( param%ncand_max, inst%ntot - 1 )

    allocate(inst%Cx(inst%ntot),inst%Cy(inst%ntot),inst%Rad(inst%ntot),&
    inst%cand(inst%ncand,inst%ntot),inst%canddist(inst%ncand,inst%ntot),&
    rowd(inst%ntot),stat=allocerr)

    if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error - inputdata.'
      stop
    end if

    if ( inst%hasdepot ) then
      read(10,*) inst%Cx(inst%idepot), inst%Cy(inst%idepot) ! depot (line 2)
      inst%Rad(inst%idepot) = 0.0d0
    end if

    do i = 1,inst%nballs
      read(10,*) inst%Cx(i), inst%Cy(i), inst%Rad(i)
      if ( inst%Rad(i) .lt. 0.0d0 ) then
        write(*,*) 'Negative radius at disk ',i,' - inputdata.'
        stop
      end if
    end do
    close(10)

    ! ---- normalization ----
    ! Coordinates and radii are divided by xymax, the largest
    ! |coordinate| + radius, so that all visit points lie in the unit box
    ! (the tolerances of the single-disk solver assume coordinates of
    ! order 1). Costs are reported in the original units (cost * xymax).
    inst%xymax = - huge( 1.0d0 )
    do k = 1,inst%ntot
      inst%xymax = max( inst%xymax, abs(inst%Cx(k)) + inst%Rad(k), &
      abs(inst%Cy(k)) + inst%Rad(k) )
    end do

    ! Degenerate instance (every disk a point at the origin).
    if ( inst%xymax .le. 0.0d0 ) inst%xymax = 1.0d0

    do k = 1,inst%ntot
      inst%Cx(k)  = inst%Cx(k)  / inst%xymax
      inst%Cy(k)  = inst%Cy(k)  / inst%xymax
      inst%Rad(k) = inst%Rad(k) / inst%xymax
    end do

    inst%vsize  = (inst%ntot + 1) * 8
    inst%lwsize = inst%vsize * 3

    ! Candidate lists: cand(1:ncand,i) holds the ncand disks closest to
    ! disk i, by increasing canddist = max(0, |c_i - c_j| - r_i - r_j), a
    ! lower bound on the length of an edge between them.
    do i = 1, inst%ntot
      ! Distance (lower bound) from disk i to every other disk j
      do j = 1, inst%ntot
        if ( j .eq. i ) then
          rowd(j) = huge( 1.0d0 )
        else
          ! max(0, distance between centers - r_i - r_j)
          dmin = max(0.0d0, norm2( (/ inst%Cx(i) - inst%Cx(j), &
                                      inst%Cy(i) - inst%Cy(j) /) ) &
                            - inst%Rad(i) - inst%Rad(j))
          rowd(j) = dmin
        end if
      end do

      do m = 1, inst%ncand
        iworst = 1
        worst = rowd(1)
        do j = 2, inst%ntot
          if ( rowd(j) .lt. worst ) then
            worst = rowd(j)
            iworst = j
          end if
        end do
        inst%cand(m, i) = iworst
        inst%canddist(m, i) = worst
        ! Invalidates the current position to find the next smallest
        rowd(iworst) = huge( 1.0d0 )
      end do
    end do
    deallocate(rowd)

    ! ---- initialize counters ----
    count%bcdmcalls = 0 ; count%bcdmcic = 0
    count%stopbcdf  = 0 ; count%stopbcdcic = 0
    count%vndcal    = 0 ; count%nscal = 0 ; count%hcalls = 0
    count%ilsiter   = 0
    count%twoptcall = 0 ; count%treeoptcall = 0 ; count%swapcall = 0
    count%insercall = 0 ; count%dbcall = 0 ; count%relocall = 0
    count%bscall    = 0 ; count%gutincall = 0
    count%cplexcalls = 0

    ! ---- incumbent record: nothing found yet ----
    count%inccost = huge( 1.0d0 )
    count%incnimp = 0
    call cpu_time( count%inct0 )

  end subroutine inputdata

  ! ------------------------------------------------------------

  subroutine drawsol(inst,xsol,plotfn,ball)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ball

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: xsol(2,inst%ntot)
    character(len=512) :: plotfn

    ! Writes a MetaPost figure with the disks, their centers, the visit
    ! points and the tour. If ball is a valid node index, its visit point
    ! is highlighted (0 = no highlight). Disks are filled from the largest
    ! to the smallest, so that smaller disks remain visible.

    ! LOCAL SCALARS
    integer :: i,ip1,k,kmax,tmp
    real(kind=8) :: xyamax
    character(len=20) :: col

    ! LOCAL ARRAYS
    integer :: order(inst%ntot)

    xyamax = - huge( 1.0d0 )
    do i = 1,inst%ntot
      xyamax = max( xyamax, 1.0d0, abs( inst%Cx(i) ), abs( inst%Cy(i) ) )
    end do

    ! Permutation of 1..ntot by non-increasing radius (selection sort)
    do i = 1,inst%ntot
      order(i) = i
    end do

    do i = 1,inst%ntot - 1
      kmax = i
      do k = i + 1,inst%ntot
        if ( inst%Rad(order(k)) .gt. inst%Rad(order(kmax)) ) kmax = k
      end do
      if ( kmax .ne. i ) then
        tmp = order(i)
        order(i) = order(kmax)
        order(kmax) = tmp
      end if
    end do

    open(unit=1030,file=plotfn)

    ! Prologue and variable declarations
    write(1030,10) inst%ntot

    ! Disk centers (w[])
    do i = 1,inst%ntot
      write(1030,20) i,inst%Cx(i) / xyamax,inst%Cy(i) / xyamax
    end do

    ! Disk radii (r[]); fullcircle has diameter 1, hence the factor 2
    do i = 1,inst%ntot
      write(1030,30) i,( inst%Rad(i) / xyamax ) * 2.0d0
    end do

    ! Visiting points (v[])
    do i = 1,inst%ntot
      write(1030,40) i,xsol(1,i) / xyamax,xsol(2,i) / xyamax
    end do

    ! ==================
    ! Draw the disks
    ! ==================

    ! All fills first, largest disk to smallest
    col = 'lavanda'
    do k = 1,inst%ntot
      i = order(k)
      write(1030,50) i,i,trim(col)
    end do

    ! All outlines only afterwards, same order
    col = 'black'
    do k = 1,inst%ntot
      i = order(k)
      write(1030,60) i,i,trim(col)
    end do

    ! ==================
    ! Draw the route
    ! ==================

    col = 'blue'
    do i = 1,inst%ntot
      ip1 = i + 1
      if ( ip1 .gt. inst%ntot ) ip1 = 1
      write(1030,70) i,ip1,trim(col)
    end do

    ! =====================
    ! Draw the disk centers
    ! =====================

    col = 'black'
    do i = 1,inst%ntot
      write(1030,80) i,trim(col)
    end do

    ! =====================
    ! Draw the visit points
    ! =====================
    
    do i = 1,inst%ntot
       if ( i .eq. ball .or. ( inst%hasdepot .and. i .eq. 1 ) ) then
          col = 'orange'
       else
          col = 'red'
       end if
       write(1030,90) i,trim(col)
    end do

    write(1030,95)
    close(1030)

    ! NON-EXECUTABLE STATEMENTS
    10 format('prologues := 3;',/, &
    'outputtemplate := "%j.mps";',/, &
    'beginfig(1);',/, &
    'u = 10.0cm;',/, &
    'pair w[];',/, &
    'pair v[];',/, &
    'numeric r[];',/, &
    'number_of_balls = ',I0,';',/, &
    'color lavanda;  lavanda = (235/255,200/255,231/255);',/, &
    'color orange;  orange = (255/255,165/255,0/255);',/, &
    'color vermelho_carmesim; vermelho_carmesim = (139/255, 0/255, 0/255);',/, &
    'color roxo_ameixa; roxo_ameixa = (93/255, 58/255, 155/255);')
    20 format('w',I0,' = ( ',f20.17,'*u,',f20.17,'*u);')
    30 format('r',I0,' = ',f20.17,'*u;')
    40 format('v',I0,' = ( ',f20.17,'*u,',f20.17,'*u);')
    50 format('fill fullcircle scaled r',I0,' shifted w',I0,' withcolor ',a,';')
    60 format('draw fullcircle scaled r',I0,' shifted w',I0, &
    ' withpen pencircle scaled 0.5pt withcolor ',a,';')
    70 format('draw v',I0,'--v',I0, &
    ' withpen pencircle scaled 0.5pt withcolor ',a,';')
    80 format('drawdot w',I0,' withpen pencircle scaled 0.75pt withcolor ',a,';')
    90 format('drawdot v',I0,' withpen pencircle scaled 1.5pt withcolor ',a,';')
    95 format('endfig;',/,'end;')

  end subroutine drawsol

  ! ------------------------------------------------------------

  subroutine dswapr(a,b)

    implicit none

    ! SCALAR ARGUMENTS
    real(kind=8), intent(inout) :: a,b

    ! On output, the values of a and b have been interchanged.
    ! Swap between REAL numbers.

    ! LOCAL SCALARS
    real(kind=8) :: temp

    temp = a
    a = b
    b = temp
  end subroutine dswapr

  ! ------------------------------------------------------------

  subroutine dswapi(a,b)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: a,b

    ! On output, the values of a and b have been interchanged.
    ! Swap between INTEGER numbers.

    ! LOCAL SCALARS
    integer :: temp

    temp = a
    a = b
    b = temp
  end subroutine dswapi

  ! ------------------------------------------------------------

  subroutine evalf(n,x,f)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(2,n)

    ! Length of the closed tour through the visit points x.

    ! LOCAL SCALARS
    integer :: k
    real(8) :: dx,dy

    f = 0.0d0

    do k = 1, n - 1
      dx = x(1,k) - x(1,k + 1)
      dy = x(2,k) - x(2,k + 1)
      f = f + hypot(dx,dy)
    end do

    dx = x(1,n) - x(1,1)
    dy = x(2,n) - x(2,1)
    f = f + hypot(dx,dy)
  end subroutine evalf

  ! ------------------------------------------------------------

  real(kind=8) function drand(ix)
    implicit none
    ! This is the random number generator of Schrage:
    !
    ! L. Schrage, A more portable Fortran random number generator, ACM
    ! Transactions on Mathematical Software 5 (1979), 132-138.
    ! SCALAR ARGUMENT
    real(kind=8), intent(inout) :: ix

    ! PARAMETERS
    real(kind=8), parameter :: a = 16807.0d0, b15 = 32768.0d0, b16 = 65536.0d0, p = 2147483647.0d0

    ! LOCAL SCALARS
    real(kind=8) :: xhi,xalo,leftlo,fhi,k

    xhi = ix / b16
    xhi = xhi - dmod(xhi, 1.0d0)
    xalo = (ix - xhi * b16) * a
    leftlo = xalo / b16
    leftlo = leftlo - dmod(leftlo, 1.0d0)
    fhi = xhi * a + leftlo
    k = fhi / b15
    k = k - dmod(k, 1.0d0)
    ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k
    if (ix .lt. 0) ix = ix + p
    drand = ix * 4.65661030875d-10

    return

  end function drand

end module instance
