! =====================================================================
!  ILS with BCD for the Close Enough Traveling Salesman Problem (CETSP).
!  Reads the instance and the parameters, builds the initial solution
!  (multi-start nearest neighbor + BCD), runs the ILS and prints the result.
! =====================================================================

program ils_cetsp

  use instance, only: instance_type,parameter_type,counter_type,inputdata,drawsol
  use metaheuristic, only: nn,ILS,WriteBestRecord
  use cetsp_fixedtour_interface, only: InitializeCPLEXTour,FinalizeCPLEXTour,GetCPLEXTotalIter
  use iso_fortran_env, only: output_unit

  implicit none

  ! LOCAL SCALARS
  integer :: allocerr,inipol,cplexstatus
  real(kind=8) :: ccost,bestcost,startt,finisht
  type(instance_type) :: inst
  type(parameter_type) :: param
  type(counter_type) :: count

  ! LOCAL ARRAYS
  integer, allocatable :: ctour(:),besttour(:)
  real(kind=8), allocatable :: csol(:,:),bestsol(:,:)
  character(len=512) :: plotfn   ! figure file name
  character(len=5)  :: engine   ! 'bcd' or 'cplex'

  ! Instance and parameters. Coordinates are normalized; costs are
  ! reported in the original units (multiplied by inst%xymax).
  call inputdata(inst,count,param)

  ! The CPLEX fixed-tour solver is initialized only if it is used.
  if ( param%usecplex ) then
    engine = 'cplex'
    write(*,'(A)') ' Optimization engine: CPLEX'
    flush( output_unit )
    call InitializeCPLEXTour(inst%ntot,inst%Cx,inst%Cy,inst%Rad,cplexstatus)
    if ( cplexstatus .ne. 0 ) then
      write(*,*) 'Error initializing the CPLEX fixed-tour solver. Status = ', &
      cplexstatus
      stop
    end if
  else
    engine = 'bcd'
    write(*,'(A)') ' Optimization engine: BCD'
    flush( output_unit )
  end if

  allocate(ctour(inst%ntot),csol(2,inst%ntot),&
  besttour(inst%ntot),bestsol(2,inst%ntot),stat=allocerr)

  if ( allocerr .ne. 0 ) then
    write(*,*) 'Allocation error - Main.'
    stop
  end if

  ! Initial solution: nearest neighbor from every node, followed by BCD;
  ! the shortest tour is kept.
  call cpu_time(startt)
  bestcost = huge( 1.0d0 )
  do inipol = 1, inst%ntot
     call nn(inst,count,inipol,ctour,csol,ccost,param%usecplex)
     if ( ccost .lt. bestcost ) then
        bestcost = ccost
        besttour = ctour
        bestsol  = csol
     end if
  end do
  call cpu_time(finisht)

  ! The initial solution is the first incumbent. Each later improvement
  ! is reported when found, so a run stopped by the time limit keeps it.
  count%inccost = bestcost
  count%incnimp = 0
  write(*,'(A,F0.6,A,F0.6)') ' Initial solution (NN+BCD/CPLEX): ',bestcost * inst%xymax,' t(s) = ',finisht - startt    
  if ( param%usecplex ) then
     write(*,'(A,I0,A,I0)')   ' CPLEX calls: ',count%cplexcalls,'   CPLEX iters: ',GetCPLEXTotalIter()
  else
     write(*,'(A,I0,A,I0)')   ' BCD calls: ',count%bcdmcalls,'   BCD sweeps: ',count%bcdmcic
  end if
  write(*,'(A)')         ' Visiting order (tour):'
  write(*,'(1x,*(I0,1x))') besttour(1:inst%ntot)
  write(*,*)
  flush( output_unit )

  ! Figure of the initial solution.
  plotfn = 'figures/plot-nnsol-'//trim(engine)//'-'//trim(inst%basename)//'.mp'
  call drawsol(inst,bestsol,plotfn,0)

  ! Record of the best solution (updated at every improvement).
  call WriteBestRecord(inst,count,bestcost,besttour,bestsol,param%usecplex, &
  0,finisht - startt,0)

  ! ILS starting from the initial solution
  call cpu_time(startt)
  call ILS(inst,param,count,besttour,bestsol,bestcost)
  call cpu_time(finisht)

  ! Results (printed only if the run stops by stagnation; a run stopped
  ! by the time limit keeps its best solution in last_incumbent/)
  write(*,*)
  write(*,'(A)')        ' ======================= RESULT ======================'
  write(*,'(A,A)')      ' Instance             : ', trim(adjustl(inst%filename))
  write(*,'(A,I0)')     ' Number of discs      : ', inst%nballs
  if ( inst%hasdepot ) then
    write(*,'(A,I0,A)') ' Depot                : yes (node ', inst%idepot, &
    ', radius 0, anchored at position 1)'
  else
    write(*,'(A,I0,A)') ' Depot                : no (disc ', inst%ianchor, &
    ' anchored at position 1)'
  end if
  write(*,'(A,F0.6)')   ' Tour cost (ILS)      : ', bestcost * inst%xymax
  write(*,'(A,F0.6,A)') ' CPU time (ILS)       : ', finisht - startt, ' s'
  write(*,'(A,I0)')     ' ILS iterations       : ', count%ilsiter
  write(*,'(A,I0)')     ' VND calls            : ', count%vndcal
  write(*,'(A,I0)')     ' N. Searchs           : ', count%nscal
  write(*,'(A,I0)')     ' N. 2-Opt             : ', count%twoptcall
  write(*,'(A,I0)')     ' N. 3-Opt             : ', count%treeoptcall
  write(*,'(A,I0)')     ' N. Swap              : ', count%swapcall
  write(*,'(A,I0)')     ' N. Insertion         : ', count%insercall
  write(*,'(A,I0)')     ' N. String Reloc      : ', count%relocall
  write(*,'(A,I0)')     ' N. Double-Bridge     : ', count%dbcall
  write(*,'(A,I0)')     ' N. Gutin             : ', count%gutincall
  write(*,'(A,I0)')     ' N. Balas-Simonetti   : ', count%bscall
  write(*,'(A,I0)')     ' HC calls             : ', count%hcalls
  write(*,'(A,I0)')     ' BCD calls (total)    : ', count%bcdmcalls
  write(*,'(A,I0)')     ' BCD sweeps (total)   : ', count%bcdmcic
  write(*,'(A,I0)')     ' CPLEX calls (total)  : ', count%cplexcalls
  write(*,'(A,I0)')     ' CPLEX iters (total)  : ', GetCPLEXTotalIter()
  write(*,'(A)')        ' Visiting order (tour):'
  write(*,'(1x,*(I0,1x))') besttour(1:inst%ntot)
  write(*,'(A)')        ' ======================================================'

  deallocate(inst%Cx,inst%Cy,inst%Rad,ctour,csol,besttour,bestsol,stat=allocerr)
  if ( allocerr .ne. 0 ) then
    write(*,*) 'Deallocation error - Main.'
    stop
  end if

  if ( param%usecplex ) call FinalizeCPLEXTour()

end program ils_cetsp
