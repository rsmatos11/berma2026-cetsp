module cetsp_fixedtour_interface

  ! ============================================================================
  ! Interface between Fortran and the C++/CPLEX solver of the fixed-tour
  ! problem (cetsp_fixedtour_interface.cpp). InitializeCPLEXTour must be
  ! called once (with the normalized Cx, Cy, Rad) before SolveTourCPLEX, and
  ! FinalizeCPLEXTour once at the end.
  ! ============================================================================

  use iso_c_binding, only: c_int,c_double

  implicit none

  public :: InitializeCPLEXTour,FinalizeCPLEXTour,SolveTourCPLEX,SolveTourCPLEXStats, &
            GetCPLEXTotalIter

  ! Total number of CPLEX barrier iterations (see GetCPLEXTotalIter).
  integer, save :: cplex_total_iter = 0

  interface

    ! Initializes the CPLEX solver with the disk data (called once).

    subroutine initialize_fixed_tour_solver(n,Cx,Cy,Rad,status) &
    bind(C,name='initialize_fixed_tour_solver')

      use iso_c_binding, only: c_int,c_double

      implicit none

      integer(c_int), value, intent(in) :: n
      real(c_double), intent(in) :: Cx(*)
      real(c_double), intent(in) :: Cy(*)
      real(c_double), intent(in) :: Rad(*)
      integer(c_int), intent(out) :: status

    end subroutine initialize_fixed_tour_solver

    ! Releases the CPLEX solver (called once).

    subroutine finalize_fixed_tour_solver() bind(C,name='finalize_fixed_tour_solver')

      implicit none

    end subroutine finalize_fixed_tour_solver

    ! Solves the fixed-tour problem for the visiting order tour(1:ntour).

    subroutine solve_fixed_tour(ntour,tour,pathcoo,result,bestbound,gap, &
    bbnodes,niterations,totaltime,status) bind(C,name='solve_fixed_tour')

      use iso_c_binding, only: c_int,c_double

      implicit none

      integer(c_int), value, intent(in) :: ntour
      integer(c_int), intent(in) :: tour(*)
      real(c_double), intent(out) :: pathcoo(2,*)
      real(c_double), intent(out) :: result
      real(c_double), intent(out) :: bestbound
      real(c_double), intent(out) :: gap
      integer(c_int), intent(out) :: bbnodes
      integer(c_int), intent(out) :: niterations
      real(c_double), intent(out) :: totaltime
      integer(c_int), intent(out) :: status

    end subroutine solve_fixed_tour

  end interface

contains

  subroutine InitializeCPLEXTour(n,Cx,Cy,Rad,status)

    ! Copies the disk data (normalized units) into the CPLEX solver.

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(out) :: status

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: Cx(n),Cy(n),Rad(n)

    ! LOCAL SCALARS
    integer(c_int) :: n_c,status_c

    ! LOCAL ARRAYS
    real(c_double) :: Cx_c(n),Cy_c(n),Rad_c(n)

    n_c  = int(n,kind=c_int)
    Cx_c = real(Cx,kind=c_double)
    Cy_c = real(Cy,kind=c_double)
    Rad_c = real(Rad,kind=c_double)

    call initialize_fixed_tour_solver(n_c,Cx_c,Cy_c,Rad_c,status_c)

    status = int(status_c)

  end subroutine InitializeCPLEXTour

  subroutine FinalizeCPLEXTour()

    ! Releases the CPLEX solver.

    implicit none

    call finalize_fixed_tour_solver()

  end subroutine FinalizeCPLEXTour

  subroutine SolveTourCPLEX(ntour,tour,x,cost,status)

    ! Optimal visit points of a fixed tour computed by CPLEX (alternative
    ! to BCD, with the same arguments x and cost). status /= 0 if CPLEX fails.

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ntour
    integer, intent(out) :: status
    real(kind=8), intent(out) :: cost

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(ntour)
    real(kind=8), intent(out) :: x(2,ntour)

    ! LOCAL SCALARS
    integer(c_int) :: ntour_c,bbnodes_c,niterations_c,status_c
    real(c_double) :: result_c,bestbound_c,gap_c,totaltime_c

    ! LOCAL ARRAYS
    integer(c_int) :: tour_c(ntour)
    real(c_double) :: pathcoo_c(2,ntour)

    ntour_c = int(ntour,kind=c_int)
    tour_c  = int(tour,kind=c_int)

    call solve_fixed_tour(ntour_c,tour_c,pathcoo_c,result_c,bestbound_c,gap_c, &
    bbnodes_c,niterations_c,totaltime_c,status_c)

    x(1:2,1:ntour) = real(pathcoo_c(1:2,1:ntour),kind=8)
    cost            = real(result_c,kind=8)
    status          = int(status_c)

    cplex_total_iter = cplex_total_iter + int(niterations_c)

  end subroutine SolveTourCPLEX

  integer function GetCPLEXTotalIter()

    ! Total number of CPLEX barrier iterations in this run.

    implicit none

    GetCPLEXTotalIter = cplex_total_iter

  end function GetCPLEXTotalIter

  subroutine SolveTourCPLEXStats(ntour,tour,x,cost,bestbound,gap,bbnodes, &
  niterations,totaltime,status)

    ! Same as SolveTourCPLEX, also returning CPLEX statistics (bound, gap,
    ! nodes, iterations, time).

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ntour
    integer, intent(out) :: bbnodes,niterations,status
    real(kind=8), intent(out) :: cost,bestbound,gap,totaltime

    ! ARRAY ARGUMENTS
    integer, intent(in) :: tour(ntour)
    real(kind=8), intent(out) :: x(2,ntour)

    ! LOCAL SCALARS
    integer(c_int) :: ntour_c,bbnodes_c,niterations_c,status_c
    real(c_double) :: result_c,bestbound_c,gap_c,totaltime_c

    ! LOCAL ARRAYS
    integer(c_int) :: tour_c(ntour)
    real(c_double) :: pathcoo_c(2,ntour)

    ntour_c = int(ntour,kind=c_int)
    tour_c  = int(tour,kind=c_int)

    call solve_fixed_tour(ntour_c,tour_c,pathcoo_c,result_c,bestbound_c,gap_c, &
    bbnodes_c,niterations_c,totaltime_c,status_c)

    x(1:2,1:ntour) = real(pathcoo_c(1:2,1:ntour),kind=8)
    cost            = real(result_c,kind=8)
    bestbound       = real(bestbound_c,kind=8)
    gap             = real(gap_c,kind=8)
    bbnodes         = int(bbnodes_c)
    niterations     = int(niterations_c)
    totaltime       = real(totaltime_c,kind=8)
    status          = int(status_c)

  end subroutine SolveTourCPLEXStats

end module cetsp_fixedtour_interface
