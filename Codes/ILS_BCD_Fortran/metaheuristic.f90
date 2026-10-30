! =====================================================================
!  module metaheuristic: nearest neighbor (nn), VND neighborhoods, VND,
!  double-bridge perturbation and ILS
! =====================================================================

module metaheuristic

  use instance, only: instance_type,parameter_type,counter_type,evalf,dswapr,dswapi,drand,drawsol
  use blockcoordinatedescent, only: BCD,mindisk
  use cetsp_fixedtour_interface, only: SolveTourCPLEX,GetCPLEXTotalIter
  use iso_fortran_env, only: output_unit

  implicit none

  ! Number of Gutin moves discarded because of an invalid assignment
  ! (expected to stay 0).
  integer(kind=8) :: ngutindesc = 0

  contains

  subroutine WriteBestRecord(inst,count,cost,tour,x,usecplex,incnum,telapsed, &
  ilsiter_val)

    ! Writes the current best solution (overwritten at each improvement):
    !   figures/<instance>-<engine>.mp          figure
    !   last_incumbent/<instance>-<engine>.csv  one-line CSV record
    !   solution_report/<instance>.txt          tour and visit points
    ! BCD and CPLEX counters are cumulative since the start of the program.

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(in) :: count

    ! SCALAR ARGUMENTS
    real(kind=8), intent(in) :: cost,telapsed
    logical, intent(in) :: usecplex
    integer, intent(in) :: incnum,ilsiter_val

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(inst%ntot)
    real(kind=8), intent(in) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    character(len=512) :: plotfn,csvfn,enginedir,reportfn
    character(len=8000) :: line
    character(len=32) :: buf
    character(len=48) :: fbuf
    integer :: bcdc,bcds,cplexc,cplexi,k,linelen,buflen

    if ( usecplex ) then
      enginedir = 'cplex'
    else
      enginedir = 'bcd'
    end if

    bcdc   = count%bcdmcalls
    bcds   = count%bcdmcic
    cplexc = count%cplexcalls
    cplexi = GetCPLEXTotalIter()

    ! ---- figure ----
    plotfn = 'figures/'//trim(inst%basename)//'-'//trim(enginedir)//'.mp'
    call drawsol(inst,x,plotfn,0)

    ! ---- one-line CSV record ----
    line = trim(inst%basename)
    write(buf,'(I0)')    inst%nballs     ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')    incnum          ; line = trim(line)//','//trim(buf)
    write(buf,'(G0)') cost*inst%xymax    ; line = trim(line)//','//trim(buf)
    write(buf,'(G0)') telapsed           ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   ilsiter_val      ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%vndcal     ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%nscal      ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%twoptcall  ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%treeoptcall; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%swapcall   ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%insercall  ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%relocall   ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%dbcall     ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%gutincall  ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%bscall     ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   count%hcalls     ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   bcdc             ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   bcds             ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   cplexc           ; line = trim(line)//','//trim(buf)
    write(buf,'(I0)')   cplexi           ; line = trim(line)//','//trim(buf)//',"'

    ! The length of line is tracked explicitly (trim() would remove the
    ! separating blank).
    linelen = len_trim(line)
    do k = 1, inst%ntot
      write(buf,'(I0)') tour(k)
      buflen = len_trim(buf)
      if ( k .gt. 1 ) then
        line(linelen+1:linelen+1) = ' '
        linelen = linelen + 1
      end if
      line(linelen+1:linelen+buflen) = buf(1:buflen)
      linelen = linelen + buflen
    end do
    line = line(1:linelen)//'"'

    csvfn = 'last_incumbent/'//trim(inst%basename)//'-'//trim(enginedir)//'.csv'
    open(unit=1032,file=trim(csvfn),status='replace',action='write')
    write(1032,'(A)') 'Instance,NDiscs,IncumbentN,TourCost,CPUTime,ILSiter,'// &
    'VNDcalls,NSearchs,N2opt,N3opt,NSwap,NInsertion,NStringReloc,'// &
    'NDoubleBridge,NGutin,NBalasSimonetti,HCcalls,BCDcalls,BCDsweeps,'// &
    'CPLEXcalls,CPLEXiter,Tour'
    write(1032,'(A)') trim(line)
    close(1032)

    ! ---- solution report: tour and visit points, in the original units,
    ! full precision (F40.17). Line k of the visit points corresponds to
    ! node tour(k); with a depot, it is the first node. ----
    reportfn = 'solution_report/'//trim(inst%basename)//'.txt'
    open(unit=1033,file=trim(reportfn),status='replace',action='write')

    write(1033,'(A,I0)') 'Number of disks = ', inst%nballs
    if ( inst%hasdepot ) then
      write(1033,'(A,I0,A)') 'Depot = yes (node ', inst%idepot, ')'
    else
      write(1033,'(A)') 'Depot = no'
    end if
    write(fbuf,'(F40.17)') cost * inst%xymax
    write(1033,'(A,A)') 'Tour length = ', trim(adjustl(fbuf))

    write(1033,'(A)',advance='no') 'Tour sequence ='
    do k = 1, inst%ntot
      write(1033,'(1X,I0)',advance='no') tour(k)
    end do
    write(1033,*)

    write(1033,'(A)') 'Tour visit points:'
    do k = 1, inst%ntot
      write(fbuf,'(F40.17)') x(1,k) * inst%xymax
      write(1033,'(2X,A)',advance='no') trim(adjustl(fbuf))
      write(fbuf,'(F40.17)') x(2,k) * inst%xymax
      write(1033,'(1X,A)') trim(adjustl(fbuf))
    end do

    close(1033)

  end subroutine WriteBestRecord

  ! ------------------------------------------------------------

  subroutine ReportIncumbent(inst,count,cost,tour,x,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    real(kind=8), intent(in) :: cost
    logical, intent(in) :: usecplex

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(inst%ntot)
    real(kind=8), intent(in) :: x(2,inst%ntot)

    ! Records a new incumbent, prints it immediately (with flush, so that a
    ! run stopped by the time limit keeps it) and calls WriteBestRecord.
    ! Costs are printed in the original units.

    ! PARAMETERS
    real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )

    ! LOCAL SCALARS
    real(kind=8) :: tnow

    if ( cost .ge. count%inccost - macheps12 * max( 1.0d0, abs(count%inccost) ) ) return

    count%inccost = cost
    count%incnimp = count%incnimp + 1

    call cpu_time( tnow )

    ! Six decimals (the improvement threshold is relative, ~1e-8).
    write(*,'(A,I0,A,F0.6,A,F0.6,A,I0)') ' [incumbent ',count%incnimp, &
    '] cost = ',cost * inst%xymax,'   t = ',tnow - count%inct0, &
    ' s   ILS iter = ',count%ilsiter
    if ( usecplex ) then
      write(*,'(A,I0,A,I0)') '   CPLEX calls = ',count%cplexcalls, &
      '   CPLEX iters = ',GetCPLEXTotalIter()
    else
      write(*,'(A,I0,A,I0)') '   BCD calls = ',count%bcdmcalls, &
      '   BCD sweeps = ',count%bcdmcic
    end if
    flush( output_unit )

    ! CPU time since the instance was read (includes the initial solution).
    call WriteBestRecord(inst,count,cost,tour,x,usecplex,count%incnimp, &
    tnow - count%inct0,count%ilsiter)

  end subroutine ReportIncumbent

  ! ------------------------------------------------------------

  subroutine nn(inst,count,inipol,tour,x,cost,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: usecplex
    integer, intent(in) :: inipol
    real(kind=8), intent(out) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(out) :: tour(inst%ntot)
    real(kind=8), intent(out) :: x(2,inst%ntot)

    ! Nearest neighbor starting at node inipol (between disk centers); the
    ! cycle is then rotated so that the anchor is at position 1, and BCD
    ! computes the visit points.

    ! LOCAL VARIABLES
    logical :: visited(inst%ntot)
    integer :: i,j,k,ipos,currentpol,nextpol,cplexstatus
    real(kind=8) :: mindist,tempdist

    ! LOCAL ARRAYS
    integer :: rottour(inst%ntot)
    real(kind=8) :: rotx(2,inst%ntot)

    count%hcalls = count%hcalls + 1

    ! Starting node (out-of-range values fall back to the anchor).
    currentpol = inipol
    if ( currentpol .lt. 1 .or. currentpol .gt. inst%ntot ) then
      currentpol = inst%ianchor
    end if

    visited(1:inst%ntot) = .false.

    tour(1)  = currentpol
    x(1:2,1) = (/ inst%Cx(currentpol), inst%Cy(currentpol) /)
    visited(currentpol) = .true.

    do i = 2, inst%ntot
      nextpol = 0
      mindist = huge( 1.0d0 )
      do j = 1, inst%ntot
        if ( .not. visited(j) ) then
          tempdist = norm2( (/ inst%Cx(currentpol), inst%Cy(currentpol) /) - &
          (/ inst%Cx(j), inst%Cy(j) /) )
          if ( tempdist .lt. mindist ) then
            mindist = tempdist
            nextpol = j
          end if
        end if
      end do

      ! Defensive check (cannot happen).
      if ( nextpol .eq. 0 ) then
        write(*,*) 'No unvisited ball left at position ',i,' - nn.'
        stop
      end if

      tour(i) = nextpol
      x(1:2,i) = (/ inst%Cx(nextpol), inst%Cy(nextpol) /)
      visited(nextpol) = .true.
      currentpol = nextpol
    end do

    ! ---- rotate the cycle so that the anchor sits at position 1 ----
    ipos = 0
    do k = 1, inst%ntot
      if ( tour(k) .eq. inst%ianchor ) then
        ipos = k
        exit
      end if
    end do

    if ( ipos .gt. 1 ) then
      do k = 1, inst%ntot
        j = mod( ipos - 2 + k, inst%ntot ) + 1
        rottour(k) = tour(j)
        rotx(1:2,k) = x(1:2,j)
      end do
      tour = rottour
      x = rotx
    end if

    call evalf(inst%ntot,x,cost)

    if ( usecplex ) then
      call SolveTourCPLEX(inst%ntot,tour,x,cost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else
      call BCD(inst,count,tour,x,cost)
    end if

  end subroutine nn
  
  ! ------------------------------------------------------------

  subroutine SwapN(inst,count,tour,x,cost,alcallbcd,bestfirst,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    integer, intent(in) :: bestfirst ! 0 = first improvement, 1 = best improvement
    logical, intent(in) :: alcallbcd,usecplex
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    logical :: stopswap
    integer :: i,j,m,p,a,b,pm,pp,scan
    real(kind=8) :: swapcost,bestcost,rad

    ! LOCAL ARRAYS
    integer :: swaptour(inst%ntot),bestour(inst%ntot),pos(inst%ntot)
    real(kind=8) :: swaptsol(2,inst%ntot),bestsol(2,inst%ntot)

    count%swapcall = count%swapcall + 1

    bestour = tour
    bestsol = x
    bestcost = cost

    do i = 1, inst%ntot
      pos(tour(i)) = i
    end do

    stopswap = .false.

    do p = 1, inst%ntot
      if ( stopswap ) exit

      pm = merge(inst%ntot, p - 1, p .eq. 1)
      pp = merge(1, p + 1, p .eq. inst%ntot)

      do scan = 1, 2
        if ( stopswap ) exit

        if ( scan .eq. 1 ) then
          j = pm
          rad = norm2( x(1:2,pm) - x(1:2,p) )
        else
          j = pp
          rad = norm2( x(1:2,p) - x(1:2,pp) )
        end if

        do m = 1, inst%ncand
          if ( inst%canddist(m,tour(j)) .ge. rad ) exit

          b = pos(inst%cand(m,tour(j)))
          a = min(p,b)
          b = max(p,b)

          ! Position 1 holds the anchor and must never be swapped.
          if ( a .lt. 2 ) cycle
          if ( b .lt. a + 2 ) cycle

          call SwapMove(inst%ntot,tour,swaptour,x,swaptsol,cost,swapcost,a,b)
          if ( alcallbcd ) then
            call OptimizePoints(swaptour,swaptsol,swapcost)
          end if
          if ( swapcost .lt. bestcost  ) then
            bestcost = swapcost
            bestour = swaptour
            bestsol = swaptsol
            ! First improvement stops at the first improving move; best
            ! improvement scans the whole (pruned) neighborhood.
            if ( bestfirst .eq. 0 ) then
              stopswap = .true.
              exit
            end if
          end if
        end do
      end do
    end do

    if ( .not.(alcallbcd) ) then
      call OptimizePoints(bestour,bestsol,bestcost)
    end if

    tour = bestour
    x = bestsol
    cost = bestcost

  contains

    ! Visit points of a fixed tour: BCD or CPLEX (usecplex).
    subroutine OptimizePoints(tvec,xvec,cvec)

      implicit none

      integer, intent(in) :: tvec(inst%ntot)
      real(kind=8), intent(inout) :: xvec(2,inst%ntot),cvec

      integer :: cplexstatus

      if ( usecplex ) then
        call SolveTourCPLEX(inst%ntot,tvec,xvec,cvec,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else
        call BCD(inst,count,tvec,xvec,cvec)
      end if

    end subroutine OptimizePoints

  end subroutine SwapN

  ! ------------------------------------------------------------

  subroutine SwapMove(nballs,tour,newtour,x,newx,cost,newcost,i,j)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i,j,nballs
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: im,ip,jm,jp
    real(kind=8) :: brokencost,addedcost

    newtour = tour
    newtour(i) = tour(j)
    newtour(j) = tour(i)

    newx = x
    newx(1:2,i) = x(1:2,j)
    newx(1:2,j) = x(1:2,i)

    ! Cyclic neighbors of positions i and j (i <= nballs-2 and j >= i+2 >= 3)
    im = merge(nballs, i - 1, i .eq. 1)
    ip = i + 1
    jm = j - 1
    jp = merge(1, j + 1, j .eq. nballs)

    if ( i .eq. 1 .and. j .eq. nballs ) then
      brokencost = edge(x, jm,j) + edge(x, i,ip)
      addedcost  = edge(x, jm,i) + edge(x, j,ip)
    else
      brokencost = edge(x, im,i) + edge(x, i,ip) + edge(x, jm,j) + edge(x, j,jp)
      addedcost  = edge(x, im,j) + edge(x, j,ip) + edge(x, jm,i) + edge(x, i,jp)
    end if

    newcost = cost + ( addedcost - brokencost )

  end subroutine SwapMove

  ! ------------------------------------------------------------

  subroutine DoubleBridgeN(inst,count,tour,x,cost,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: usecplex
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    integer :: i,j,i1,j1,i2,j2,typpe,cplexstatus
    real(kind=8), parameter :: infty = huge( 1.0d0 )
    real(kind=8) :: best_t,costc,costd,newcost

    ! LOCAL ARRAYS
    integer :: pred_t1(2),pred_t2(2),newtour(inst%ntot),pred_reach_d(inst%ntot),&
    pred_reach_c(inst%ntot),pred_cross_d(inst%ntot,2),pred_cross_c(inst%ntot,2)
    real(kind=8) :: newx(2,inst%ntot),best_reach_d(inst%ntot),best_reach_c(inst%ntot),&
    best_cross_d(inst%ntot),best_cross_c(inst%ntot),edgelen(inst%ntot)

    count%dbcall = count%dbcall + 1

    do i = 1, inst%ntot - 1
      edgelen(i) = norm2( x(1:2,i) - x(1:2,i + 1) )
    end do
    edgelen(inst%ntot) = norm2( x(1:2,inst%ntot) - x(1:2,1) )

    pred_t1(1:2) = -1
    pred_t2(1:2) = -1

    best_t = infty
    do j = 3, inst%ntot - 1
      best_reach_d(j) = cost_d(inst%ntot,1,j,x,edgelen)
      pred_reach_d(j) = 1

      best_reach_c(j) = cost_c(inst%ntot,1,j,x,edgelen)
      pred_reach_c(j) = 1

      costc = best_reach_c(j)
      if ( costc .lt. best_t ) then
        best_t = costc
        pred_t2 = (/ 1, j /)
        typpe = 0
      end if
    end do

    do i = 2, inst%ntot - 2
      best_cross_d(i + 1) = best_reach_d(i + 1)
      pred_cross_d(i + 1,1:2) = (/ pred_reach_d(i + 1), i + 1 /)

      best_cross_c(i + 1) = best_reach_c(i + 1)
      pred_cross_c(i + 1,1:2) = (/ pred_reach_c(i + 1), i + 1 /)

      do j = i + 2, inst%ntot
        costd = cost_d(inst%ntot,i,j,x,edgelen)
        costc = cost_c(inst%ntot,i,j,x,edgelen)

        if ( j .le. inst%ntot - 1 ) then
          if ( best_reach_d(j) .lt. best_cross_d(j - 1) ) then
            best_cross_d(j) = best_reach_d(j)
            pred_cross_d(j,1:2) = (/ pred_reach_d(j), j /)
          else
            best_cross_d(j) = best_cross_d(j - 1)
            pred_cross_d(j,1:2) = pred_cross_d(j - 1,1:2)
          end if

          if ( best_reach_c(j) .lt. best_cross_c(j - 1) ) then
            best_cross_c(j) = best_reach_c(j)
            pred_cross_c(j,1:2) = (/ pred_reach_c(j), j /)
          else
            best_cross_c(j) = best_cross_c(j - 1)
            pred_cross_c(j,1:2) = pred_cross_c(j - 1,1:2)
          end if

          if ( costd .lt. best_reach_d(j) ) then
            best_reach_d(j) = costd
            pred_reach_d(j) = i
          end if

          if ( costc .lt. best_reach_c(j) ) then
            best_reach_c(j) = costc
            pred_reach_c(j) = i
          end if
        end if

        if ( costd + best_cross_d(j - 1) .lt. best_t ) then
          best_t = costd + best_cross_d(j - 1)
          pred_t1 = pred_cross_d(j - 1,1:2)
          pred_t2 = (/ i, j /)
          typpe = 1
        end if

        if ( costd + best_cross_c(j - 1) .lt. best_t ) then
          best_t = costd + best_cross_c(j - 1)
          pred_t1 = pred_cross_c(j - 1,1:2)
          pred_t2 = (/ i, j /)
          typpe = 2
        end if

        if ( costc + best_cross_d(j - 1) .lt. best_t ) then
          best_t = costc + best_cross_d(j - 1)
          pred_t1 = pred_cross_d(j - 1,1:2)
          pred_t2 = (/ i, j /)
          typpe = 3
        end if

        if ( costc .lt. best_t ) then
          best_t = costc
          pred_t2 = (/ i, j /)
          typpe = 0
        end if
      end do
    end do

    if ( best_t .ge. -1.0d-8 ) return

    if ( typpe .ne. 0 ) then
      i1 = pred_t1(1)
      j1 = pred_t1(2)
    end if
    i2 = pred_t2(1)
    j2 = pred_t2(2)

    if ( typpe .eq. 0 ) then
      call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,i2,j2)
    else if ( typpe .eq. 1 ) then
      call DoubleBridgeMove(inst%ntot,tour,newtour,x,newx,cost,newcost,i1,i2,j1,j2)
    else if ( typpe .eq. 2 .or. typpe .eq. 3 ) then
      call Type2Move(inst%ntot,tour,newtour,x,newx,cost,newcost,i1,i2,j1,j2,typpe)
    end if

    if ( usecplex ) then
      call SolveTourCPLEX(inst%ntot,newtour,newx,newcost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else
      call BCD(inst,count,newtour,newx,newcost)
    end if

    if ( newcost .lt. cost ) then
      tour = newtour
      x = newx
      cost = newcost
    end if
  end subroutine DoubleBridgeN

  ! ------------------------------------------------------------

  pure function edge(x,a,b) result(dab)
    real(kind=8), intent(in) :: x(:,:)
    integer, intent(in) :: a,b
    real(kind=8) :: dab,dx,dy
    dx = x(1,a) - x(1,b)
    dy = x(2,a) - x(2,b)
    dab = hypot(dx,dy)
  end function edge

  ! ------------------------------------------------------------

  subroutine DoubleBridgeMove(nballs,tour,newtour,x,newx,cost,newcost,i,j,k,ell)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i,j,k,ell,nballs
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: ellp
    real(kind=8) :: brokencost,addedcost

    newtour(1:i) = tour(1:i)
    newtour(i + 1:i + ell - k) = tour(k + 1:ell)
    newtour(i + ell - k + 1:i + ell - j) = tour(j + 1:k)
    newtour(i + ell - j + 1:ell) = tour(i + 1:j)
    newtour(ell + 1:nballs) = tour(ell + 1:nballs)

    newx(1:2,1:i) = x(1:2,1:i)
    newx(1:2,i + 1:i + ell - k) = x(1:2,k + 1:ell)
    newx(1:2,i + ell - k + 1:i + ell - j) = x(1:2,j + 1:k)
    newx(1:2,i + ell - j + 1:ell) = x(1:2,i + 1:j)
    newx(1:2,ell + 1:nballs) = x(1:2,ell + 1:nballs)

    ellp = merge(1, ell + 1, ell .eq. nballs)

    brokencost = edge(x, i,i + 1) + edge(x, j,j + 1) + edge(x, k,k + 1) + edge(x, ell,ellp)
    addedcost = edge(x, i,k + 1) + edge(x, ell,j + 1) + edge(x, k,i + 1) + edge(x, j,ellp)

    newcost = cost + ( addedcost - brokencost )

  end subroutine DoubleBridgeMove

  ! ------------------------------------------------------------

  real(kind=8)  function cost_d(nballs,i,j,x,edgelen)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nballs,i,j

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(2,nballs),edgelen(nballs)

    ! LOCAL SCALAR
    integer :: inext,jnext

    inext = merge(1, i + 1, i == nballs)
    jnext = merge(1, j + 1, j == nballs)

    cost_d = ( norm2( x(1:2,i) - x(1:2,jnext) ) + norm2( x(1:2,inext) - x(1:2,j) ) ) &
    - ( edgelen(i) + edgelen(j) )

    return

  end function cost_d

  ! ------------------------------------------------------------

  real(kind=8)  function cost_c(nballs,i,j,x,edgelen)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nballs,i,j

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(2,nballs),edgelen(nballs)

    ! LOCAL SCALAR
    integer :: inext,jnext

    inext = merge(1, i + 1, i == nballs)
    jnext = merge(1, j + 1, j == nballs)

    cost_c = ( norm2( x(1:2,i) - x(1:2,j) ) + norm2( x(1:2,inext) - x(1:2,jnext) ) ) &
    - ( edgelen(i) + edgelen(j) )

    return

  end function cost_c

  ! ------------------------------------------------------------

  subroutine Type2Move(nballs,tour,newtour,x,newx,cost,newcost,i,j,k,ell,typpe)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i,j,k,ell,nballs,typpe
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: ellp
    real(kind=8) :: brokencost,addedcost

    ellp = merge(1, ell + 1, ell .eq. nballs)

    brokencost = edge(x, i,i + 1) + edge(x, j,j + 1) + edge(x, k,k + 1) + edge(x, ell,ellp)
    addedcost = brokencost

    if ( typpe .eq. 2 ) then
      newtour(1:i) = tour(1:i)
      newtour(i + 1:i + k - j) = tour(k:j + 1:-1) ! C reversed
      newtour(i + k - j + 1:i + ell - j) = tour(ell:k + 1:-1) ! D reversed
      newtour(i + ell - j + 1:ell) = tour(i + 1:j) ! B forward
      newtour(ell + 1:nballs) = tour(ell + 1:nballs)

      newx(1:2,1:i) = x(1:2,1:i)
      newx(1:2,i + 1:i + k - j) = x(1:2,k:j + 1:-1)
      newx(1:2,i + k - j + 1:i + ell - j) = x(1:2,ell:k + 1:-1)
      newx(1:2,i + ell - j + 1:ell) = x(1:2,i + 1:j)
      newx(1:2,ell + 1:nballs) = x(1:2,ell + 1:nballs)

      addedcost = edge(x, i,k) + edge(x, j + 1,ell) + edge(x, k + 1,i + 1) + edge(x, j,ellp)

    else if (typpe .eq. 3) then
      newtour(1:i) = tour(1:i)
      newtour(i + 1:i + ell - k) = tour(k + 1:ell) ! D forward
      newtour(i + ell - k + 1:i + ell - k + j - i) = tour(j:i + 1:-1) ! B reversed
      newtour(i + ell - k + j - i + 1:ell) = tour(k:j + 1:-1) ! C reversed
      newtour(ell + 1:nballs) = tour(ell + 1:nballs)

      newx(1:2,1:i) = x(1:2,1:i)
      newx(1:2,i + 1:i + ell - k) = x(1:2,k + 1:ell)
      newx(1:2,i + ell - k + 1:i + ell - k + j - i) = x(1:2,j:i + 1:-1)
      newx(1:2,i + ell - k + j - i + 1:ell) = x(1:2,k:j + 1:-1)
      newx(1:2,ell + 1:nballs) = x(1:2,ell + 1:nballs)

      addedcost = edge(x, i,k + 1) + edge(x, ell,j) + edge(x, i + 1,k) + edge(x, j + 1,ellp)
    end if

    newcost = cost + ( addedcost - brokencost )

  end subroutine Type2Move

  ! ------------------------------------------------------------

  subroutine G2opt(inst,count,tour,x,cost,alcallbcd,bestfirst,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    integer, intent(in) :: bestfirst ! 0 = first improvement, 1 = best improvement
    logical, intent(in) :: alcallbcd,usecplex
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    logical :: stop2opt
    integer :: i,j,m,a,b,ii,inext,iprev,jprev,polii
    real(kind=8) :: twoptcost,bestcost,ledge

    ! LOCAL ARRAYS
    integer :: twoptour(inst%ntot),bestour(inst%ntot),pos(inst%ntot)
    real(kind=8) :: twoptsol(2,inst%ntot),bestsol(2,inst%ntot)

    count%twoptcall = count%twoptcall + 1

    bestour = tour
    bestsol = x
    bestcost = cost

    do i = 1, inst%ntot
      pos(tour(i)) = i
    end do

    stop2opt = .false.
    ii = 1
    do while ( ii .le. inst%ntot .and. .not.(stop2opt) )

      polii = tour(ii)

      ! ---------------------------------------------------------
      ! Scan A: the new edge leaves position ii
      ! ---------------------------------------------------------
      inext = merge(1, ii + 1, ii .eq. inst%ntot)
      ledge = norm2( x(1:2,ii) - x(1:2,inext) )

      m = 1
      do while ( m .le. inst%ncand .and. .not.(stop2opt) )
        if ( inst%canddist(m,polii) .ge. ledge ) exit

        j = pos(inst%cand(m,polii))
        a = min(ii,j)
        b = max(ii,j)

        ! With a = 1 the 2-opt reverses positions 2..b, so the anchor at
        ! position 1 is never moved.
        if ( b .ge. a + 2 .and. .not.( a .eq. 1 .and. b .eq. inst%ntot ) ) then
          call TwoOptMove(inst%ntot,tour,twoptour,x,twoptsol,cost,twoptcost,a,b)
          if ( alcallbcd ) then
            call OptimizePoints(twoptour,twoptsol,twoptcost)
          end if
          if ( twoptcost .lt. bestcost  ) then
            bestcost = twoptcost
            bestour = twoptour
            bestsol = twoptsol
            ! First improvement stops at the first improving move; best
            ! improvement scans the whole (pruned) neighborhood.
            if ( bestfirst .eq. 0 ) stop2opt = .true.
          end if
        end if
        m = m + 1
      end do

      ! ---------------------------------------------------------
      ! Scan B: the new edge arrives at position ii
      ! ---------------------------------------------------------
      if ( .not.(stop2opt) ) then
        iprev = merge(inst%ntot, ii - 1, ii .eq. 1)
        ledge = norm2( x(1:2,iprev) - x(1:2,ii) )

        m = 1
        do while ( m .le. inst%ncand .and. .not.(stop2opt) )
          if ( inst%canddist(m,polii) .ge. ledge ) exit

          j = pos(inst%cand(m,polii))
          jprev = merge(inst%ntot, j - 1, j .eq. 1)
          a = min(iprev,jprev)
          b = max(iprev,jprev)

          if ( b .ge. a + 2 .and. .not.( a .eq. 1 .and. b .eq. inst%ntot ) ) then
            call TwoOptMove(inst%ntot,tour,twoptour,x,twoptsol,cost,twoptcost,a,b)
            if ( alcallbcd ) then
              call OptimizePoints(twoptour,twoptsol,twoptcost)
            end if
            if ( twoptcost .lt. bestcost  ) then
              bestcost = twoptcost
              bestour = twoptour
              bestsol = twoptsol
              ! First improvement stops at the first improving move; best
              ! improvement scans the whole (pruned) neighborhood.
              if ( bestfirst .eq. 0 ) stop2opt = .true.
            end if
          end if
          m = m + 1
        end do
      end if

      ii = ii + 1
    end do

    if ( .not.(alcallbcd) ) then
      call OptimizePoints(bestour,bestsol,bestcost)
    end if

    tour = bestour
    x = bestsol
    cost = bestcost

  contains

    ! Visit points of a fixed tour: BCD or CPLEX (usecplex).
    subroutine OptimizePoints(tvec,xvec,cvec)

      implicit none

      integer, intent(in) :: tvec(inst%ntot)
      real(kind=8), intent(inout) :: xvec(2,inst%ntot),cvec

      integer :: cplexstatus

      if ( usecplex ) then
        call SolveTourCPLEX(inst%ntot,tvec,xvec,cvec,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else
        call BCD(inst,count,tvec,xvec,cvec)
      end if

    end subroutine OptimizePoints

  end subroutine G2opt

  ! ------------------------------------------------------------

  subroutine TwoOptMove(nballs,tour,newtour,x,newx,cost,newcost,i,j)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i,j,nballs
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: ip,jp
    real(kind=8) :: brokencost,addedcost

    newtour(1:i) = tour(1:i)
    newtour(i + 1:j) = tour(j:i + 1:-1)
    newtour(j + 1:nballs) = tour(j + 1:nballs)

    newx(1:2,1:i) = x(1:2,1:i)
    newx(1:2,i + 1:j) = x(1:2,j:i + 1:-1)
    newx(1:2,j + 1:nballs) = x(1:2,j + 1:nballs)

    ! Cyclic wrap-around (i <= nballs-2)
    ip = i + 1
    jp = merge(1, j + 1, j .eq. nballs)

    brokencost = edge(x, i,ip) + edge(x, j,jp)
    addedcost  = edge(x, i,j)  + edge(x, ip,jp)

    newcost = cost + ( addedcost - brokencost )

  end subroutine TwoOptMove

  ! ------------------------------------------------------------

  subroutine G3opt(inst,count,tour,x,cost,alcallbcd,bestfirst,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    integer, intent(in) :: bestfirst ! 0 = first improvement, 1 = best improvement
    logical, intent(in) :: alcallbcd,usecplex
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! PARAMETERS
    real(kind=8), parameter :: tol = 1.0d-12

    ! LOCAL SCALARS
    logical :: stop3opt
    integer :: i,j,k,ip,jp,kp,m1,m2,mm,mtype,itmp
    integer :: t1,t2,t3,t4,t5,t6,dd1,dd2,dd3,dir1,dir2,dir3
    integer :: b1,b2,b3
    real(kind=8) :: bestcost,hatcost,cost1,cost2,cost3,cost4,newcost
    real(kind=8) :: dd,g1,g2,gain

    ! LOCAL ARRAYS
    integer :: newtour(inst%ntot),bestour(inst%ntot),hatour(inst%ntot), &
    tour1(inst%ntot),tour2(inst%ntot),tour3(inst%ntot),tour4(inst%ntot), &
    pos(inst%ntot),addp(2,3),cnd(2,3)
    real(kind=8) :: bestsol(2,inst%ntot),hatsol(2,inst%ntot),sol1(2,inst%ntot), &
    sol2(2,inst%ntot),sol3(2,inst%ntot),sol4(2,inst%ntot),newx(2,inst%ntot)

    count%treeoptcall = count%treeoptcall + 1

    bestour = tour
    bestsol = x
    bestcost = cost
    stop3opt = .false.

    do i = 1, inst%ntot
      pos(tour(i)) = i
    end do

    ! =====================================================================
    ! DEGENERATE CASES (j == i+1 or k == j+1)
    ! =====================================================================
    deg_outer1: do i = 1, inst%ntot - 4
      j = i + 1
      do k = j + 1, inst%ntot - 1
        call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,j,k)
        call TwoOptMove(inst%ntot,newtour,hatour,newx,hatsol,newcost,hatcost,i,k)
        if ( alcallbcd ) call OptimizePoints(hatour,hatsol,hatcost)
        if ( hatcost .lt. bestcost ) then
          bestcost = hatcost
          bestour = hatour
          bestsol = hatsol
          ! First improvement stops at the first improving move; best
          ! improvement scans the whole (pruned) neighborhood.
          if ( bestfirst .eq. 0 ) then
            stop3opt = .true.
            exit deg_outer1
          end if
        end if
      end do
    end do deg_outer1

    if ( .not. stop3opt ) then
      deg_outer2: do i = 1, inst%ntot - 4
        do j = i + 2, inst%ntot - 2
          k = j + 1
          call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,i,j)
          call TwoOptMove(inst%ntot,newtour,hatour,newx,hatsol,newcost,hatcost,i,k)
          if ( alcallbcd ) call OptimizePoints(hatour,hatsol,hatcost)
          if ( hatcost .lt. bestcost ) then
            bestcost = hatcost
            bestour = hatour
            bestsol = hatsol
            ! First improvement stops at the first improving move; best
            ! improvement scans the whole (pruned) neighborhood.
            if ( bestfirst .eq. 0 ) then
              stop3opt = .true.
              exit deg_outer2
            end if
          end if
        end do
      end do deg_outer2
    end if

    ! =====================================================================
    ! SEQUENTIAL SEARCH over the chain t1 ... t6
    ! =====================================================================
    if ( .not. stop3opt ) then

      seq_outer: do t1 = 1, inst%ntot
        do dd1 = 1, 2
          dir1 = merge(1,-1,dd1 .eq. 1)
          t2 = nbpos(t1,dir1)
          dd = edge(x,t1,t2)

          do m1 = 1, inst%ncand
            if ( inst%canddist(m1,tour(t2)) .ge. dd ) exit
            t3 = pos(inst%cand(m1,tour(t2)))
            g1 = dd - edge(x,t2,t3)
            if ( g1 .le. tol ) cycle

            do dd2 = 1, 2
              dir2 = merge(1,-1,dd2 .eq. 1)
              t4 = nbpos(t3,dir2)

              do m2 = 1, inst%ncand
                if ( inst%canddist(m2,tour(t4)) .ge. edge(x,t3,t4) + g1 ) exit
                t5 = pos(inst%cand(m2,tour(t4)))
                g2 = g1 + edge(x,t3,t4) - edge(x,t4,t5)
                if ( g2 .le. tol ) cycle

                do dd3 = 1, 2
                  dir3 = merge(1,-1,dd3 .eq. 1)
                  t6 = nbpos(t5,dir3)

                  gain = g2 + edge(x,t5,t6) - edge(x,t6,t1)
                  if ( gain .le. 1.0d-9 ) cycle

                  b1 = brkpos(t1,dir1)
                  b2 = brkpos(t3,dir2)
                  b3 = brkpos(t5,dir3)
                  i = b1; j = b2; k = b3
                  if ( i .gt. j ) then
                    itmp = i; i = j; j = itmp
                  end if
                  if ( j .gt. k ) then
                    itmp = j; j = k; k = itmp
                  end if
                  if ( i .gt. j ) then
                    itmp = i; i = j; j = itmp
                  end if

                  if ( k .gt. inst%ntot - 1 ) cycle
                  if ( j .lt. i + 2 .or. k .lt. j + 2 ) cycle

                  ip = i + 1
                  jp = j + 1
                  kp = k + 1

                  addp(1,1) = t2; addp(2,1) = t3
                  addp(1,2) = t4; addp(2,2) = t5
                  addp(1,3) = t6; addp(2,3) = t1

                  mtype = 0
                  do mm = 1, 4
                    if ( mm .eq. 1 ) then
                      cnd(1,1)=i;  cnd(2,1)=j
                      cnd(1,2)=ip; cnd(2,2)=k
                      cnd(1,3)=jp; cnd(2,3)=kp
                    else if ( mm .eq. 2 ) then
                      cnd(1,1)=i;  cnd(2,1)=jp
                      cnd(1,2)=k;  cnd(2,2)=j
                      cnd(1,3)=ip; cnd(2,3)=kp
                    else if ( mm .eq. 3 ) then
                      cnd(1,1)=i;  cnd(2,1)=k
                      cnd(1,2)=jp; cnd(2,2)=ip
                      cnd(1,3)=j;  cnd(2,3)=kp
                    else
                      cnd(1,1)=i;  cnd(2,1)=jp
                      cnd(1,2)=k;  cnd(2,2)=ip
                      cnd(1,3)=j;  cnd(2,3)=kp
                    end if
                    if ( pairin(addp(1,1),addp(2,1),cnd) .and. &
                         pairin(addp(1,2),addp(2,2),cnd) .and. &
                         pairin(addp(1,3),addp(2,3),cnd) ) then
                      mtype = mm
                      exit
                    end if
                  end do

                  if ( mtype .eq. 0 ) cycle

                  if ( mtype .eq. 1 ) then
                    call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,i,j)
                    call TwoOptMove(inst%ntot,newtour,tour1,newx,sol1,newcost,cost1,j,k)
                    if ( alcallbcd ) call OptimizePoints(tour1,sol1,cost1)
                    if ( cost1 .lt. bestcost ) then
                      bestcost = cost1; bestour = tour1; bestsol = sol1
                      ! First improvement stops at the first improving move; best
                      ! improvement scans the whole (pruned) neighborhood.
                      if ( bestfirst .eq. 0 ) then
                        stop3opt = .true.
                        exit seq_outer
                      end if
                    end if
                  else if ( mtype .eq. 2 ) then
                    call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,j,k)
                    call TwoOptMove(inst%ntot,newtour,tour2,newx,sol2,newcost,cost2,i,k)
                    if ( alcallbcd ) call OptimizePoints(tour2,sol2,cost2)
                    if ( cost2 .lt. bestcost ) then
                      bestcost = cost2; bestour = tour2; bestsol = sol2
                      ! First improvement stops at the first improving move; best
                      ! improvement scans the whole (pruned) neighborhood.
                      if ( bestfirst .eq. 0 ) then
                        stop3opt = .true.
                        exit seq_outer
                      end if
                    end if
                  else if ( mtype .eq. 3 ) then
                    call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,i,j)
                    call TwoOptMove(inst%ntot,newtour,tour3,newx,sol3,newcost,cost3,i,k)
                    if ( alcallbcd ) call OptimizePoints(tour3,sol3,cost3)
                    if ( cost3 .lt. bestcost ) then
                      bestcost = cost3; bestour = tour3; bestsol = sol3
                      ! First improvement stops at the first improving move; best
                      ! improvement scans the whole (pruned) neighborhood.
                      if ( bestfirst .eq. 0 ) then
                        stop3opt = .true.
                        exit seq_outer
                      end if
                    end if
                  else
                    call TwoOptMove(inst%ntot,tour,newtour,x,newx,cost,newcost,j,k)
                    call TwoOptMove(inst%ntot,newtour,tour2,newx,sol2,newcost,cost2,i,k)
                    call TwoOptMove(inst%ntot,tour2,tour4,sol2,sol4,cost2,cost4,i + k - j,k)
                    if ( alcallbcd ) call OptimizePoints(tour4,sol4,cost4)
                    if ( cost4 .lt. bestcost ) then
                      bestcost = cost4; bestour = tour4; bestsol = sol4
                      ! First improvement stops at the first improving move; best
                      ! improvement scans the whole (pruned) neighborhood.
                      if ( bestfirst .eq. 0 ) then
                        stop3opt = .true.
                        exit seq_outer
                      end if
                    end if
                  end if

                end do
              end do
            end do
          end do
        end do
      end do seq_outer
    end if

    if ( .not.(alcallbcd) ) then
      call OptimizePoints(bestour,bestsol,bestcost)
    end if

    tour = bestour
    x = bestsol
    cost = bestcost

  contains

    ! Visit points of a fixed tour: BCD or CPLEX (usecplex).
    subroutine OptimizePoints(tvec,xvec,cvec)

      implicit none

      integer, intent(in) :: tvec(inst%ntot)
      real(kind=8), intent(inout) :: xvec(2,inst%ntot),cvec

      integer :: cplexstatus

      if ( usecplex ) then
        call SolveTourCPLEX(inst%ntot,tvec,xvec,cvec,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else
        call BCD(inst,count,tvec,xvec,cvec)
      end if

    end subroutine OptimizePoints

    pure function nbpos(p,d) result(q)
      integer, intent(in) :: p,d
      integer :: q
      if ( d .gt. 0 ) then
        q = merge(1, p + 1, p .eq. inst%ntot)
      else
        q = merge(inst%ntot, p - 1, p .eq. 1)
      end if
    end function nbpos

    pure function brkpos(p,d) result(b)
      integer, intent(in) :: p,d
      integer :: b
      if ( d .gt. 0 ) then
        b = p
      else
        b = merge(inst%ntot, p - 1, p .eq. 1)
      end if
    end function brkpos

    pure function pairin(p,q,s) result(ok)
      integer, intent(in) :: p,q,s(2,3)
      logical :: ok
      ok = any( ( s(1,:) .eq. p .and. s(2,:) .eq. q ) .or. &
                ( s(1,:) .eq. q .and. s(2,:) .eq. p ) )
    end function pairin

  end subroutine G3opt

  ! ------------------------------------------------------------

  subroutine Ginsert(inst,count,tour,x,cost,alcallbcd,bestfirst,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    integer, intent(in) :: bestfirst ! 0 = first improvement, 1 = best improvement
    logical, intent(in) :: alcallbcd,usecplex
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    logical :: stopinsert
    integer :: i,j,m,a,ap,im,ip,pass
    real(kind=8) :: insertcost,bestcost,rprev,rnext,rad

    ! LOCAL ARRAYS
    integer :: insertour(inst%ntot),bestour(inst%ntot),pos(inst%ntot)
    real(kind=8) :: insertsol(2,inst%ntot),bestsol(2,inst%ntot)

    count%insercall = count%insercall + 1

    bestour = tour
    bestsol = x
    bestcost = cost

    do i = 1, inst%ntot
      pos(tour(i)) = i
    end do

    stopinsert = .false.

    ! Position 1 holds the anchor; it is protected by never relocating it.
    do i = 2, inst%ntot
      if ( stopinsert ) exit

      im = i - 1
      ip = merge(1, i + 1, i .eq. inst%ntot)
      rprev = norm2( x(1:2,im) - x(1:2,i) )
      rnext = norm2( x(1:2,i) - x(1:2,ip) )
      rad = max( rprev, rnext )

      do m = 1, inst%ncand
        if ( stopinsert ) exit
        if ( inst%canddist(m,tour(i)) .ge. rad ) exit

        ap = pos(inst%cand(m,tour(i)))

        do pass = 1, 2
          if ( pass .eq. 1 ) then
            if ( inst%canddist(m,tour(i)) .ge. rprev ) cycle
            a = ap
          else
            if ( inst%canddist(m,tour(i)) .ge. rnext ) cycle
            a = merge(inst%ntot, ap - 1, ap .eq. 1)
          end if

          if ( a .eq. i .or. a .eq. im ) cycle

          if ( a .gt. i ) then
            j = a
          else
            j = a + 1
          end if
          
          ! Nothing is inserted at position 1, which would push the anchor to 2
          if ( j .lt. 2 .or. j .gt. inst%ntot .or. j .eq. i ) cycle

          call InsertMove(inst%ntot,tour,insertour,x,insertsol,cost,insertcost,i,j)
          if ( alcallbcd ) then
            call OptimizePoints(insertour,insertsol,insertcost)
          end if
          if ( insertcost .lt. bestcost  ) then
            bestcost = insertcost
            bestour = insertour
            bestsol = insertsol
            ! First improvement stops at the first improving move; best
            ! improvement scans the whole (pruned) neighborhood.
            if ( bestfirst .eq. 0 ) then
              stopinsert = .true.
              exit
            end if
          end if
        end do
      end do
    end do

    if ( .not.(alcallbcd) ) then
      call OptimizePoints(bestour,bestsol,bestcost)
    end if

    tour = bestour
    x = bestsol
    cost = bestcost

  contains

    ! Visit points of a fixed tour: BCD or CPLEX (usecplex).
    subroutine OptimizePoints(tvec,xvec,cvec)

      implicit none

      integer, intent(in) :: tvec(inst%ntot)
      real(kind=8), intent(inout) :: xvec(2,inst%ntot),cvec

      integer :: cplexstatus

      if ( usecplex ) then
        call SolveTourCPLEX(inst%ntot,tvec,xvec,cvec,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else
        call BCD(inst,count,tvec,xvec,cvec)
      end if

    end subroutine OptimizePoints

  end subroutine Ginsert

  ! ------------------------------------------------------------

  subroutine InsertMove(nballs,tour,newtour,x,newx,cost,newcost,i,j)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i,j,nballs
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: im,ip,jm,jp
    real(kind=8) :: brokencost,addedcost

    if ( j .lt. i ) then
      newtour(1:j - 1) = tour(1:j - 1)
      newtour(j) = tour(i)
      newtour(j + 1:i) = tour(j:i - 1)
      newtour(i + 1:nballs) = tour(i + 1:nballs)

      newx(1:2,1:j - 1) = x(1:2,1:j - 1)
      newx(1:2,j) = x(1:2,i)
      newx(1:2,j + 1:i) = x(1:2,j:i - 1)
      newx(1:2,i + 1:nballs) = x(1:2,i + 1:nballs)
    else
      newtour(1:i - 1) = tour(1:i - 1)
      newtour(i:j - 1) = tour(i + 1:j)
      newtour(j) = tour(i)
      newtour(j + 1:nballs) = tour(j + 1:nballs)

      newx(1:2,1:i - 1) = x(1:2,1:i - 1)
      newx(1:2,i:j - 1) = x(1:2,i + 1:j)
      newx(1:2,j) = x(1:2,i)
      newx(1:2,j + 1:nballs) = x(1:2,j + 1:nballs)
    end if

    if ( ( i .eq. 1 .and. j .eq. nballs ) .or. ( i .eq. nballs .and. j .eq. 1 ) ) then
      newcost = cost
    else if ( j .gt. i ) then
      im = merge(nballs, i - 1, i .eq. 1)
      ip = i + 1
      jp = merge(1, j + 1, j .eq. nballs)
      brokencost = edge(x, im,i) + edge(x, i,ip) + edge(x, j,jp)
      addedcost  = edge(x, im,ip) + edge(x, j,i) + edge(x, i,jp)
      newcost = cost + ( addedcost - brokencost )
    else
      jm = merge(nballs, j - 1, j .eq. 1)
      im = i - 1
      ip = merge(1, i + 1, i .eq. nballs)
      brokencost = edge(x, jm,j) + edge(x, im,i) + edge(x, i,ip)
      addedcost  = edge(x, jm,i) + edge(x, i,j) + edge(x, im,ip)
      newcost = cost + ( addedcost - brokencost )
    end if

  end subroutine InsertMove

  ! ------------------------------------------------------------

  subroutine Gutin(inst,count,tour,x,cost,seed,cardZ,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: usecplex
    integer, intent(in) :: cardZ
    real(kind=8), intent(inout) :: cost,seed

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    integer, parameter :: NMAX = 1000
    real(kind=8), parameter :: SCALE_FACTOR = 10000.0d0
    integer :: i,j,k,s,T,prev,next,num_objects,shift,r,temp,idx,nmov,cplexstatus
    real(kind=8) :: newcost,gmin

    ! LOCAL ARRAYS
    integer :: newtour(inst%ntot),A(NMAX,NMAX + 1),C(NMAX),Z(inst%ntot),&
    zbool(inst%ntot - 1),objects(inst%ntot - 1),config(inst%ntot - 1)
    real(kind=8) :: xmin(2),newx(2,inst%ntot),xvisit(2,cardZ,cardZ)

    ! EXTERNAL SUBROUTINES
    external :: ASSCT

    count%gutincall = count%gutincall + 1

    ! Position 1 (the anchor) is fixed, so only positions 2..ntot can be
    ! selected (nmov of them).
    nmov = inst%ntot - 1

    ! assct.f uses fixed-size arrays with NMAX = 1000, so s <= NMAX.
    s = min(cardZ, nmov / 2, NMAX)
    num_objects = nmov - s

    objects(1:num_objects) = 0
    if (s .gt. 0) objects(1:s) = 1

    do i = num_objects, 2, -1
      r = floor(drand(seed) * i) + 1
      temp = objects(i)
      objects(i) = objects(r)
      objects(r) = temp
    end do

    idx = 1
    do i = 1, num_objects
      if (objects(i) .eq. 1) then
        config(idx) = 1
        config(idx + 1) = 0
        idx = idx + 2
      else
        config(idx) = 0
        idx = idx + 1
      end if
    end do

    shift = floor(drand(seed) * nmov)
    zbool(1:nmov) = 0

    do i = 1, nmov
      idx = mod(i - 1 + shift, nmov) + 1
      zbool(idx) = config(i)
    end do

    s = 0
    do i = 1, nmov
      if (zbool(i) .eq. 1) then
        s = s + 1
        Z(s) = i + 1
      end if
    end do

    do j = 1, s
      prev = Z(j) - 1
      next = merge(1, Z(j) + 1, Z(j) == inst%ntot)
      do i = 1, s
        if ( i .eq. j ) then
          xvisit(1:2,i,j) = x(1:2,Z(i))
          gmin = norm2( xvisit(1:2,i,j) - x(1:2,prev) ) + &
          norm2( xvisit(1:2,i,j) - x(1:2,next) )
        else
          xmin = x(1:2,Z(i))
          call mindisk(inst,tour(Z(i)),x(1:2,prev),x(1:2,next),xmin,gmin)
          xvisit(1:2,i,j) = xmin(1:2)
        end if
        A(i,j) = nint(gmin * SCALE_FACTOR)
      end do
    end do

    ! Assignment costs are >= 0; a negative value (NaN/overflow) would make
    ! ASSCT return an invalid assignment, so the move is discarded.
    if ( any( A(1:s,1:s) .lt. 0 ) .or. any( A(1:s,1:s) .gt. 1000000000 ) ) then
      call GutinDiscard('invalid cost A(i,j)')
      return
    end if

    call ASSCT(s, A, C, T)

    ! Invalid assignment: discard the move instead of indexing Z and
    ! xvisit out of bounds.
    if ( any( C(1:s) .lt. 1 ) .or. any( C(1:s) .gt. s ) ) then
      call GutinDiscard('ASSCT returned C(k) outside 1..s')
      return
    end if
    
    newtour = tour
    newx = x
    do k = 1, s
      newtour(Z(k)) = tour(Z(C(k)))
      newx(1:2,Z(k)) = xvisit(1:2,C(k),k)
    end do

    call evalf(inst%ntot,newx,newcost)

    if ( usecplex ) then
      call SolveTourCPLEX(inst%ntot,newtour,newx,newcost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else
      call BCD(inst,count,newtour,newx,newcost)
    end if

    if ( newcost .lt. cost ) then
      tour = newtour
      x = newx
      cost = newcost
    end if

  contains

    subroutine GutinDiscard(reason)
      ! Counts a discarded Gutin move and reports it on the 1st, 10th, 100th,
      ! 1000th and every 10000th occurrence.
      character(len=*), intent(in) :: reason
      ngutindesc = ngutindesc + 1
      if ( ngutindesc .eq. 1 .or. ngutindesc .eq. 10 .or. &
           ngutindesc .eq. 100 .or. ngutindesc .eq. 1000 .or. &
           mod(ngutindesc,10000_8) .eq. 0 ) then
        write(*,'(A,A,A,I0)') ' [Gutin warning] move discarded (', &
        reason,')  total = ',ngutindesc
        flush( output_unit )
      end if
    end subroutine GutinDiscard

  end subroutine Gutin

  ! ------------------------------------------------------------

  subroutine StringRelocN(inst,count,tour,x,cost,strlen,alcallbcd,reverse,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: alcallbcd,usecplex
    integer, intent(in) :: strlen,reverse
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    integer :: i,j,k,kmax,inv_idx,num_inv,m,ss,se,jj,cnt_j
    real(kind=8) :: bestcost,srcost,rstart,rend
    logical :: is_inplace

    ! LOCAL ARRAYS
    integer :: bestour(inst%ntot),srtour(inst%ntot),pos(inst%ntot),jlist(4*inst%ncand)
    real(kind=8) :: string(2,strlen),bestsol(2,inst%ntot),srx(2,inst%ntot)
    logical :: do_inv(2)

    count%relocall = count%relocall + 1

    bestour = tour
    bestsol = x
    bestcost = cost

    do i = 1, inst%ntot
      pos(tour(i)) = i
    end do

    num_inv = 1
    do_inv(1) = .false.
    if ( reverse .eq. 1 ) then
      num_inv = 2
      do_inv(2) = .true.
    end if

    ! The anchor is never part of the string (i >= 2)
    do i = 2, inst%ntot - 1
      kmax = min(strlen - 1, inst%ntot - i)
      string(:,1:kmax + 1) = x(:,i:i + kmax)
      do k = 1, kmax
        ss = i
        se = i + k
        
        rstart = norm2( x(1:2,i - 1) - x(1:2,ss) )
        rend   = norm2( x(1:2,se) - x(1:2,merge(1,se + 1,se.eq.inst%ntot)) )

        cnt_j = 0
        do m = 1, inst%ncand
          if ( inst%canddist(m,tour(ss)) .ge. rstart ) exit
          jj = pos(inst%cand(m,tour(ss)))
          call pushj(jj)
          call pushj( merge(inst%ntot, jj - 1, jj.eq.1) )
        end do
        do m = 1, inst%ncand
          if ( inst%canddist(m,tour(se)) .ge. rend ) exit
          jj = pos(inst%cand(m,tour(se)))
          call pushj( merge(inst%ntot, jj - 1, jj.eq.1) )
          call pushj(jj)
        end do

        do m = 1, cnt_j
          j = jlist(m)
          if ( j .lt. i .or. j .gt. i + k ) then
            ! The insertion never happens at position 1
            if ( j .lt. 2 .or. j .gt. inst%ntot ) cycle 
            
            is_inplace = (j == i - 1) .or. (i == 1 .and. j == inst%ntot)
            do inv_idx = 1, num_inv
              if ( is_inplace .and. .not. do_inv(inv_idx) ) cycle
              call StringMove(inst%ntot,tour,srtour,x,srx,cost,srcost,i,k,j,&
              string,strlen,do_inv(inv_idx))
              if ( alcallbcd ) then
                call OptimizePoints(srtour,srx,srcost)
              end if
              if ( srcost .lt. bestcost  ) then
                bestcost = srcost
                bestour = srtour
                bestsol = srx
              end if
            end do
          end if
        end do
      end do
    end do

    if ( .not.(alcallbcd) ) then
      call OptimizePoints(bestour,bestsol,bestcost)
    end if

    tour = bestour
    x = bestsol
    cost = bestcost

  contains

    subroutine pushj(jv)
      integer, intent(in) :: jv
      integer :: t
      do t = 1, cnt_j
        if ( jlist(t) .eq. jv ) return
      end do
      cnt_j = cnt_j + 1
      jlist(cnt_j) = jv
    end subroutine pushj

    ! Visit points of a fixed tour: BCD or CPLEX (usecplex).
    subroutine OptimizePoints(tvec,xvec,cvec)

      implicit none

      integer, intent(in) :: tvec(inst%ntot)
      real(kind=8), intent(inout) :: xvec(2,inst%ntot),cvec

      integer :: cplexstatus

      if ( usecplex ) then
        call SolveTourCPLEX(inst%ntot,tvec,xvec,cvec,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else
        call BCD(inst,count,tvec,xvec,cvec)
      end if

    end subroutine OptimizePoints

  end subroutine StringRelocN

  ! ------------------------------------------------------------

  subroutine StringMove(nballs,tour,newtour,x,newx,cost,newcost,i,k,j,string,strlen,inverter)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in) :: inverter
    integer, intent(in) :: nballs,i,k,j,strlen
    real(kind=8), intent(in) :: cost
    real(kind=8), intent(out) :: newcost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(nballs)
    integer, intent(out) :: newtour(nballs)
    real(kind=8), intent(in) :: x(2,nballs),string(2,strlen)
    real(kind=8), intent(out) :: newx(2,nballs)

    ! LOCAL SCALARS
    integer :: leftout,rightout,leftin,rightin,strstart,strend
    real(kind=8) :: brokencost,addedcost,dx,dy,delta

    if ( j .lt. i ) then
      newtour(1:j) = tour(1:j)
      newx(1:2,1:j) = x(1:2,1:j)
      if ( inverter ) then
        newtour(j + 1:j + k + 1) = tour(i + k:i:-1)
        newx(1:2,j + 1:j + k + 1) = string(1:2,k + 1:1:-1)
      else
        newtour(j + 1:j + k + 1) = tour(i:i + k)
        newx(1:2,j + 1:j + k + 1) = string(1:2,1:k + 1)
      end if
      newtour(j + k + 2:i + k) = tour(j + 1:i - 1)
      newtour(i + k + 1:nballs) = tour(i + k + 1:nballs)
      newx(1:2,j + k + 2:i + k) = x(1:2,j + 1:i - 1)
      newx(1:2,i + k + 1:nballs) = x(1:2,i + k + 1:nballs)
    else
      newtour(1:i - 1) = tour(1:i - 1)
      newtour(i:j - k - 1) = tour(i + k + 1:j)
      newx(1:2,1:i - 1) = x(1:2,1:i - 1)
      newx(1:2,i:j - k - 1) = x(1:2,i + k + 1:j)
      if ( inverter ) then
        newtour(j - k:j) = tour(i + k:i:-1)
        newx(1:2,j - k:j) = string(1:2,k + 1:1:-1)
      else
        newtour(j - k:j) = tour(i:i + k)
        newx(1:2,j - k:j) = string(1:2,1:k + 1)
      end if
      newtour(j + 1:nballs) = tour(j + 1:nballs)
      newx(1:2,j + 1:nballs) = x(1:2,j + 1:nballs)
    end if

    leftout  = merge(nballs, i - 1, i == 1)
    rightout = merge(1, i + k + 1, i + k == nballs)
    strstart = i
    strend   = i + k
    leftin   = j
    rightin  = merge(1, j + 1, j == nballs)

    if (leftin .eq. leftout .or. leftin .eq. strend) then
      if (.not. inverter) then
        delta = 0.0d0
      else
        dx = x(1,leftout) - x(1,strstart)
        dy = x(2,leftout) - x(2,strstart)
        brokencost = hypot(dx,dy)
        dx = x(1,strend) - x(1,rightout)
        dy = x(2,strend) - x(2,rightout)
        brokencost = brokencost + hypot(dx,dy)
        dx = x(1,leftout) - x(1,strend)
        dy = x(2,leftout) - x(2,strend)
        addedcost = hypot(dx,dy)
        dx = x(1,strstart) - x(1,rightout)
        dy = x(2,strstart) - x(2,rightout)
        addedcost = addedcost + hypot(dx,dy)
        delta = addedcost - brokencost
      end if
      newcost = cost + delta
      return
    end if

    dx = x(1,leftout) - x(1,strstart)
    dy = x(2,leftout) - x(2,strstart)
    brokencost = hypot(dx,dy)
    dx = x(1,strend) - x(1,rightout)
    dy = x(2,strend) - x(2,rightout)
    brokencost = brokencost + hypot(dx,dy)
    dx = x(1,leftin) - x(1,rightin)
    dy = x(2,leftin) - x(2,rightin)
    brokencost = brokencost + hypot(dx,dy)

    dx = x(1,leftout) - x(1,rightout)
    dy = x(2,leftout) - x(2,rightout)
    addedcost = hypot(dx,dy)

    if (.not. inverter) then
      dx = x(1,leftin) - x(1,strstart)
      dy = x(2,leftin) - x(2,strstart)
      addedcost = addedcost + hypot(dx,dy)
      dx = x(1,strend) - x(1,rightin)
      dy = x(2,strend) - x(2,rightin)
      addedcost = addedcost + hypot(dx,dy)
    else
      dx = x(1,leftin) - x(1,strend)
      dy = x(2,leftin) - x(2,strend)
      addedcost = addedcost + hypot(dx,dy)
      dx = x(1,strstart) - x(1,rightin)
      dy = x(2,strstart) - x(2,rightin)
      addedcost = addedcost + hypot(dx,dy)
    end if

    delta = addedcost - brokencost
    newcost = cost + delta

  end subroutine StringMove

  ! ------------------------------------------------------------

  subroutine BalaSimonetti(inst,count,tour,x,cost,bsk,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: usecplex
    integer, intent(in) :: bsk
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)       
    real(kind=8), intent(inout) :: x(2,inst%ntot)   

    ! LOCAL SCALARS
    integer ::  k,i,nnodes,lsize,s,t,cplexstatus
    real(kind=8) :: djcost,newcost

    ! LOCAL ARRAYS 
    integer :: djtour(inst%ntot),citour(inst%ntot),newtour(inst%ntot),&
    V(inst%vsize),L(inst%lwsize),vcit(inst%vsize),prev_ptr(inst%vsize)
    real(kind=8) :: newsol(2,inst%ntot),W(inst%lwsize)

    count%bscall = count%bscall + 1

    nnodes = -1

    if ( bsk .eq. 2 ) then
      nnodes = 3 * inst%ntot - 1
      lsize = 5 * inst%ntot - 4
      call AdjLists(inst,x,nnodes,lsize,V,L,W,vcit,2)
    else if ( bsk .eq. 3 ) then
      nnodes = 8 * inst%ntot - 6
      lsize = 17 * inst%ntot - 22
      call AdjLists(inst,x,nnodes,lsize,V,L,W,vcit,3)
    end if

    s = 1 ! Source node
    t = nnodes ! Sink node
    call ShortestPath(inst,s,t,nnodes,V,L,W,prev_ptr,djcost)

    if ( djcost .ge. huge(1.0d0) ) return

    djtour(1) = s
    djtour(inst%ntot) = prev_ptr(t)
    do i = inst%ntot - 1, 2, -1
      djtour(i) = prev_ptr(djtour(i + 1))
    end do

    do k = 1, inst%ntot
      citour(k) = vcit(djtour(k))
    end do
    do k = 1, inst%ntot
      newtour(k) = tour(citour(k))
      newsol(1:2,k) = x(1:2,citour(k))
    end do

    call evalf(inst%ntot,newsol,newcost)

    if ( usecplex ) then
      call SolveTourCPLEX(inst%ntot,newtour,newsol,newcost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else
      call BCD(inst,count,newtour,newsol,newcost)
    end if

    if ( newcost .lt. cost ) then
      tour = newtour
      x = newsol
      cost = newcost
    end if
  end subroutine BalaSimonetti

  ! ------------------------------------------------------------

  subroutine AdjLists(inst,x,nnodes,lsize,V,L,W,vcit,algo)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nnodes,lsize,algo

    ! ARRAY ARGUMENTS
    integer, intent(out) :: V(inst%vsize),L(inst%lwsize),vcit(inst%vsize)
    real(kind=8), intent(in) :: x(2,inst%ntot) 
    real(kind=8), intent(out) :: W(inst%lwsize)

    ! LOCAL SCALARS
    integer ::  k,i,j
    real(kind=8), parameter :: infty = huge( 1.0d0 )

    if ( algo .eq. 2 ) then

      vcit(1) = 1 
      k = 2 
      do i = 2, inst%ntot - 1
        vcit(k) = i
        vcit(k + 1) = i + 1
        vcit(k + 2) = i - 1
        k = k + 3
      end do
      vcit(k) = inst%ntot
      vcit(k + 1) = 1
      vcit(k + 2) = inst%ntot - 1
      vcit(nnodes) = 1

      V(1) = 1
      L(1) = 2
      L(2) = 3
      W(1) = dist23(1, L(1))
      W(2) = dist23(1, L(2))

      k = 2 
      j = 3 

      do i = 2, inst%ntot - 1
        V(k) = j
        L(j) = k + 3
        L(j + 1) = k + 4
        W(j) = dist23(k, L(j))
        W(j + 1) = dist23(k, L(j + 1))

        j = j + 2

        V(k + 1) = j
        L(j) = k + 5
        W(j) = dist23(k + 1, L(j))

        j = j + 1

        V(k + 2) = j
        L(j) = k + 3
        L(j + 1) = k + 4
        W(j) = dist23(k + 2, L(j))
        W(j + 1) = dist23(k + 2, L(j + 1))

        j = j + 2
        k = k + 3
      end do

      V(k) = j
      L(j) = k + 3
      W(j) = dist23(k, L(j))

      j = j + 1

      V(k + 1) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 2) = j
      L(j) = k + 3
      W(j) = dist23(k + 2, L(j))

      V(nnodes) = lsize
      L(lsize) = -1
      W(lsize) = infty

    else if ( algo .eq. 3 ) then
      
      ! Layer 1
      vcit(1) = 1
      k = 2 

      ! Layer 2
      vcit(k) = 2
      vcit(k + 1) = 3
      vcit(k + 2) = 1
      vcit(k + 3) = 4
      vcit(k + 4) = 3
      vcit(k + 5) = 1
      vcit(k + 6) = 2
      vcit(k + 7) = 1 
      k = k + 8

      do i = 3, inst%ntot - 2
        vcit(k) = i
        vcit(k + 1) = i + 1
        vcit(k + 2) = i - 1
        vcit(k + 3) = i + 2
        vcit(k + 4) = i + 1
        vcit(k + 5) = i - 1
        vcit(k + 6) = i
        vcit(k + 7) = i - 2
        k = k + 8
      end do

      ! Layer inst%ntot-1
      vcit(k) = inst%ntot - 1
      vcit(k + 1) = inst%ntot
      vcit(k + 2) = inst%ntot - 2
      vcit(k + 3) = 1 
      vcit(k + 4) = inst%ntot
      vcit(k + 5) = inst%ntot - 2
      vcit(k + 6) = inst%ntot - 1
      vcit(k + 7) = inst%ntot - 3
      k = k + 8

      ! Layer inst%ntot
      vcit(k) = inst%ntot
      vcit(k + 1) = 1
      vcit(k + 2) = inst%ntot - 1
      vcit(k + 3) = 1
      vcit(k + 4) = 1
      vcit(k + 5) = inst%ntot - 1
      vcit(k + 6) = inst%ntot
      vcit(k + 7) = inst%ntot - 2

      ! Layer inst%ntot+1
      vcit(nnodes) = 1

      V(1) = 1
      L(1) = 2
      L(2) = 3
      L(3) = 5
      W(1) = dist23(1, L(1))
      W(2) = dist23(1, L(2))
      W(3) = dist23(1, L(3))

      k = 2 
      j = 4 
      do i = 2, inst%ntot - 1
        V(k) = j
        L(j) = k + 8
        L(j + 1) = k + 9
        L(j + 2) = k + 11
        W(j) = dist23(k, L(j))
        W(j + 1) = dist23(k, L(j + 1))
        W(j + 2) = dist23(k, L(j + 2))

        j = j + 3

        V(k + 1) = j
        L(j) = k + 10
        L(j + 1) = k + 12
        W(j) = dist23(k + 1, L(j))
        W(j + 1) = dist23(k + 1, L(j + 1))

        j = j + 2

        V(k + 2) = j
        L(j) = k + 8
        L(j + 1) = k + 9
        L(j + 2) = k + 11
        W(j) = dist23(k + 2, L(j))
        W(j + 1) = dist23(k + 2, L(j + 1))
        W(j + 2) = dist23(k + 2, L(j + 2))

        j = j + 3

        V(k + 3) = j
        L(j) = k + 13
        L(j + 1) = k + 14
        W(j) = dist23(k + 3, L(j))
        W(j + 1) = dist23(k + 3, L(j + 1))

        j = j + 2

        V(k + 4) = j
        L(j) = k + 15
        W(j) = dist23(k + 4, L(j))

        j = j + 1

        V(k + 5) = j
        L(j) = k + 10
        L(j + 1) = k + 12
        W(j) = dist23(k + 5, L(j))
        W(j + 1) = dist23(k + 5, L(j + 1))

        j = j + 2

        V(k + 6) = j
        L(j) = k + 15
        W(j) = dist23(k + 6, L(j))

        j = j + 1

        V(k + 7) = j
        L(j) = k + 8
        L(j + 1) = k + 9
        L(j + 2) = k + 11
        W(j) = dist23(k + 7, L(j))
        W(j + 1) = dist23(k + 7, L(j + 1))
        W(j + 2) = dist23(k + 7, L(j + 2))

        k = k + 8
        j = j + 3
      end do

      ! Layer Vn
      V(k) = j
      L(j) = k + 8
      W(j) = dist23(k, L(j))

      j = j + 1

      V(k + 1) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 2) = j
      L(j) = k + 8
      W(j) = dist23(k + 2, L(j))

      j = j + 1

      V(k + 3) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 4) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 5) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 6) = j
      L(j) = -1
      W(j) = infty

      j = j + 1

      V(k + 7) = j
      L(j) = k + 8
      W(j) = dist23(k + 7, L(j))

      ! Layer V_n+1
      V(nnodes) = lsize
      L(lsize) = -1
      W(lsize) = infty
    end if

    contains
      pure function dist23(a,b) result(d)
        integer, intent(in) :: a,b
        real(kind=8) :: d,dx,dy
        dx = x(1, vcit(a)) - x(1, vcit(b))
        dy = x(2, vcit(a)) - x(2, vcit(b))
        d = hypot(dx,dy)
      end function dist23

  end subroutine AdjLists

  ! ------------------------------------------------------------

  subroutine ShortestPath(inst,s,t,nnodes,V,L,W,prev_ptr,cost)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst

    ! SCALAR ARGUMENTS
    integer, intent(in) :: s,t,nnodes
    real(kind=8), intent(out) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: V(inst%vsize),L(inst%lwsize)
    integer, intent(out) :: prev_ptr(inst%vsize)
    real(kind=8), intent(in) :: W(inst%lwsize)

    ! LOCAL SCALARS
    integer ::  i,j,k
    real(kind=8), parameter :: infty = huge( 1.0d0 )
    real(kind=8) :: label

    ! LOCAL ARRAYS
    real(kind=8) :: adist(inst%vsize)

    adist(1:nnodes) = infty
    prev_ptr(1:nnodes) = -1
    adist(s) = 0.0d0

    do i = 1, nnodes - 1
      if ( adist(i) .eq. infty ) cycle
      do k = V(i), V(i + 1) - 1
        j = L(k)
        if ( j .gt. 0 ) then
          label = adist(i) + W(k)
          if ( label .lt. adist(j) ) then
            adist(j) = label
            prev_ptr(j) = i
          end if
        end if
      end do
    end do
    cost = adist(t)
  end subroutine ShortestPath

  ! ------------------------------------------------------------

  subroutine NeighborhoodSearch(inst,param,count,jviz,tour,x,cost,alcallbcd)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(parameter_type), intent(inout) :: param
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: alcallbcd
    integer, intent(in) :: jviz
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    count%nscal = count%nscal + 1

    if (jviz .eq. 1) then
      call BalaSimonetti(inst,count,tour,x,cost,param%bsk,param%usecplex)
    else if (jviz .eq. 2) then
      call DoubleBridgeN(inst,count,tour,x,cost,param%usecplex)
    else if (jviz .eq. 3) then
      call Gutin(inst,count,tour,x,cost,param%seed,param%cardZ,param%usecplex)
    else if (jviz .eq. 4) then
      call G2opt(inst,count,tour,x,cost,alcallbcd,param%bestfirst,param%usecplex)
    else if (jviz .eq. 5) then
      call SwapN(inst,count,tour,x,cost,alcallbcd,param%bestfirst,param%usecplex)
    else if (jviz .eq. 6) then
      call Ginsert(inst,count,tour,x,cost,alcallbcd,param%bestfirst,param%usecplex)
    else if (jviz .eq. 7) then
      call StringRelocN(inst,count,tour,x,cost,param%strlen,alcallbcd,param%reverse,param%usecplex)
    else if (jviz .eq. 8) then
      call G3opt(inst,count,tour,x,cost,alcallbcd,param%bestfirst,param%usecplex)
    end if

  end subroutine NeighborhoodSearch

  ! ------------------------------------------------------------

  subroutine PerturbDB(inst,count,tour,x,cost,seed,numdbmoves,usecplex)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    logical, intent(in) :: usecplex
    integer, intent(in) :: numdbmoves
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALAR
    integer :: i,j,iter,cplexstatus
    real(kind=8) :: newcost

    ! LOCAL ARRAY
    integer :: newtour(inst%ntot),orddb(4)
    real(kind=8) :: newx(2,inst%ntot)

    do iter = 1, numdbmoves

      orddb(1) = floor((inst%ntot - 2) * drand(seed)) + 2
      orddb(2) = floor((inst%ntot - 2) * drand(seed)) + 2
      orddb(3) = floor((inst%ntot - 2) * drand(seed)) + 2
      orddb(4) = floor((inst%ntot - 2) * drand(seed)) + 2

      do while ( (abs( orddb(1) - orddb(2) ) .lt. 1) .or. &
        (abs( orddb(1) - orddb(3) ) .lt. 1) .or. &
        (abs( orddb(1) - orddb(4) ) .lt. 1) .or. &
        (abs( orddb(2) - orddb(3) ) .lt. 1) .or. &
        (abs( orddb(2) - orddb(4) ) .lt. 1) .or. &
        (abs( orddb(3) - orddb(4) ) .lt. 1) )

        orddb(1) = floor((inst%ntot - 2) * drand(seed)) + 2
        orddb(2) = floor((inst%ntot - 2) * drand(seed)) + 2
        orddb(3) = floor((inst%ntot - 2) * drand(seed)) + 2
        orddb(4) = floor((inst%ntot - 2) * drand(seed)) + 2
      end do

      do i = 1,3
        do j = i + 1,4
          if ( orddb(j) .lt. orddb(i) ) then
            call dswapi(orddb(j),orddb(i))
          end if
        end do
      end do

      call DoubleBridgeMove(inst%ntot,tour,newtour,x,newx,cost,newcost,&
      orddb(1),orddb(2),orddb(3),orddb(4))

      tour = newtour
      x = newx
      cost = newcost
    end do

    if ( usecplex ) then
      call SolveTourCPLEX(inst%ntot,tour,x,cost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else
      call BCD(inst,count,tour,x,cost)
    end if

  end subroutine PerturbDB

  ! ------------------------------------------------------------

  subroutine VND(inst,param,count,tour,x,cost,alcallbcd,vnd_type)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(parameter_type), intent(inout) :: param
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    integer, intent(in) :: vnd_type
    logical, intent(in) :: alcallbcd
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! LOCAL SCALARS
    integer :: i,j,k
    real(kind=8) :: barcost,hatcost
    real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )

    ! LOCAL ARRAYS
    integer :: bartour(inst%ntot),hatour(inst%ntot),neigvec(param%jmax)
    real(kind=8) :: xbar(2,inst%ntot),xhat(2,inst%ntot),privec(param%jmax)

    count%vndcal = count%vndcal + 1

    ! Order of the neighborhoods: 0 = fixed order 1..jmax; 1 = random order
    ! at each VND call.
    if ( vnd_type .eq. 0 ) then
      do k = 1, param%jmax
        neigvec(k) = k
      end do
    else if ( vnd_type .eq. 1 ) then
      do k = 1, param%jmax
        privec(k) = drand(param%seed)
        neigvec(k) = k
      end do
      do i = 1, param%jmax - 1
        do k = i + 1, param%jmax
          if ( privec(k) .gt. privec(i) ) then
            call dswapr(privec(k),privec(i))
            call dswapi(neigvec(k),neigvec(i))
          end if
        end do
      end do
    end if

    bartour = tour
    xbar = x
    barcost = cost

    j = param%jmin

    do while ( j .le. param%jmax )
      hatour(1:inst%ntot) = bartour(1:inst%ntot)
      xhat = xbar
      hatcost = barcost

      call NeighborhoodSearch(inst,param,count,neigvec(j),hatour,xhat,hatcost,alcallbcd)

      if ( .not. all( hatour(1:inst%ntot) .eq. bartour(1:inst%ntot) ) .and. &
      hatcost .le. barcost - macheps12 * max( 1.0d0, abs(barcost) ) ) then
        bartour = hatour
        xbar = xhat
        barcost = hatcost
        ! Local improvement; ReportIncumbent checks whether it also improves
        ! the best solution.
        call ReportIncumbent(inst,count,barcost,bartour,xbar,param%usecplex)
        j = param%jmin
      else
        j = j + 1
      end if
    end do

    tour = bartour
    x = xbar
    cost = barcost

  end subroutine VND

  ! ------------------------------------------------------------

  subroutine ILS(inst,param,count,tour,x,cost)

    implicit none

    ! TYPE ARGUMENTS
    type(instance_type), intent(inout) :: inst
    type(parameter_type), intent(inout) :: param
    type(counter_type), intent(inout) :: count

    ! SCALAR ARGUMENTS
    real(kind=8), intent(inout) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: tour(inst%ntot)
    real(kind=8), intent(inout) :: x(2,inst%ntot)

    ! PARAMETERS
    real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )

    ! LOCAL SCALARS
    logical :: alcallbcd
    logical :: cooladapted
    integer :: kappa,ell,noimprov,inipol,ncool,status_env,cplexstatus
    integer :: ellmax
    integer :: ihalf
    real(kind=8) :: telap,thalf,iterremain,start_time,finish_time,tmax
    real(kind=8) :: bestcost,hatcost,acost,localcost
    character(len=256) :: dumcar

    ! LOCAL ARRAYS
    integer :: bestour(inst%ntot),hatour(inst%ntot),atour(inst%ntot),localtour(inst%ntot)
    real(kind=8) :: xbest(2,inst%ntot),xhat(2,inst%ntot),xatual(2,inst%ntot)

    call cpu_time(start_time)
       
    ! Time limit (environment variable ilstlim).
    tmax = 0.0d0
    call get_environment_variable('ilstlim', dumcar)
    if ( len_trim(dumcar) .gt. 0 ) then
        read(dumcar,*,iostat=status_env) tmax
        if ( status_env .ne. 0 ) tmax = 0.0d0
    end if

    alcallbcd = .true.
    count%maxbcdcic = param%maxbcdcic
    call VND(inst,param,count,tour,x,cost,alcallbcd,param%vnd_type)

    if ( param%usecplex ) then
      call SolveTourCPLEX(inst%ntot,tour,x,cost,cplexstatus)
      count%cplexcalls = count%cplexcalls + 1
      if ( cplexstatus .ne. 0 ) then
        write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
        stop
      end if
    else if ( param%refine .eq. 1 ) then
      count%maxbcdcic = 200
      call BCD(inst,count,tour,x,cost)
      count%maxbcdcic = param%maxbcdcic
    end if

    ! Recomputes the cost (removes accumulated rounding).
    call evalf(inst%ntot,x,cost)

    bestour = tour 
    xbest = x
    bestcost = cost

    atour = tour 
    xatual = x
    acost = cost

    localtour = tour
    localcost = cost

    noimprov = 0
    ell = 0
    kappa = 1
    cooladapted = .false.
    ihalf = 0
    thalf = 0.0d0

    ! Stagnation limit only for ntot <= 200 (ell_np_max consecutive
    ! iterations without improvement); larger instances stop by the time
    ! limit.
    if ( inst%ntot .le. 200 ) then
      ellmax = param%ell_np_max
    else
      ellmax = huge(ellmax)
    end if

    do while ( ell .lt. ellmax )
      count%ilsiter = count%ilsiter + 1

      ! Threshold reductions spread over the time limit
      if ( .not. cooladapted .and. tmax .gt. 0.0d0 ) then
          call cpu_time(finish_time)
          telap = finish_time - start_time

          if ( ihalf .eq. 0 ) then
              if ( count%ilsiter .ge. param%ncool_probe/2 .or. &
              telap .ge. 0.5d0 * param%cool_probe_frac * tmax ) then
                  ihalf = count%ilsiter
                  thalf = telap
              end if
          else if ( count%ilsiter .gt. ihalf .and. &
          ( count%ilsiter .ge. param%ncool_probe .or. &
          telap .ge. param%cool_probe_frac * tmax ) ) then

              if ( param%cool_marginal .and. telap .gt. thalf ) then
                  iterremain = dble( count%ilsiter - ihalf ) * &
                  ( tmax - telap ) / ( telap - thalf )
              else if ( telap .gt. 0.0d0 ) then
                  iterremain = dble( count%ilsiter ) * &
                  ( tmax - telap ) / telap
              else
                  iterremain = 0.0d0
              end if

              ncool = 1
              do while ( dble(2*ncool) .le. iterremain / dble(param%ncool_steps) )
                  ncool = 2 * ncool
              end do
              param%kappamax = ncool
              kappa = 0
              cooladapted = .true.
          end if
      end if

      hatour = atour
      xhat = xatual
      hatcost = acost

      if ( noimprov .eq. 3 ) then
        noimprov = 0
        inipol = floor( inst%ntot * drand(param%seed) ) + 1
        call nn(inst,count,inipol,hatour,xhat,hatcost,param%usecplex)
      else
        call PerturbDB(inst,count,hatour,xhat,hatcost,param%seed,param%numdbmoves,param%usecplex)
      end if

      call VND(inst,param,count,hatour,xhat,hatcost,alcallbcd,param%vnd_type)

      if ( param%usecplex ) then
        call SolveTourCPLEX(inst%ntot,hatour,xhat,hatcost,cplexstatus)
        count%cplexcalls = count%cplexcalls + 1
        if ( cplexstatus .ne. 0 ) then
          write(*,*) 'Error solving the fixed-tour CPLEX model. Status = ',cplexstatus
          stop
        end if
      else if ( param%refine .eq. 1 ) then
        count%maxbcdcic = 200
        call BCD(inst,count,hatour,xhat,hatcost)
        count%maxbcdcic = param%maxbcdcic
      end if

      call evalf(inst%ntot,xhat,hatcost)

      if ( hatcost .le. bestcost - macheps12 * max( 1.0d0, abs(bestcost) ) ) then
        bestour = hatour
        xbest = xhat
        bestcost = hatcost
        ! The complete BCD after VND may improve the best solution.
        call ReportIncumbent(inst,count,bestcost,bestour,xbest,param%usecplex)
        ell = 0
        noimprov = 0
      else
        if ( abs( hatcost - localcost ) .le. macheps12 * max( 1.0d0, abs(localcost) ) ) then
          if ( all( hatour(1:inst%ntot) .eq. localtour(1:inst%ntot) ) ) then
            noimprov = noimprov + 1 
          end if
        end if
        ell = ell + 1
      end if

      if ( (hatcost .le. acost - macheps12 * max( 1.0d0, abs(acost) )) .or. &
      (hatcost .le. (1.0d0 + param%eps) * bestcost) ) then
        atour = hatour
        xatual = xhat
        acost = hatcost

        localtour = atour
        localcost = acost
      end if

      kappa = kappa + 1

      if ( kappa .eq. param%kappamax ) then
        kappa = 0
        param%eps = param%h * param%eps
      end if
    end do

    tour = bestour
    x = xbest
    cost = bestcost

  end subroutine ILS

end module metaheuristic
