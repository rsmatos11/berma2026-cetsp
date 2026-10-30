/*
 * Replacement for cetsp_fixedtour_interface.cpp when CPLEX is not available.
 *
 * It provides the three functions called by cetsp_fixedtour_interface.f90,
 * so that the program can be compiled and run with the BCD engine
 * (usecplex = 0) without CPLEX. If the CPLEX engine is requested
 * (usecplex = 1), initialization fails with a clear message.
 */
#include <stdio.h>

void initialize_fixed_tour_solver(int n, const double* Cx, const double* Cy,
                                  const double* Rad, int* status)
{
  (void)n; (void)Cx; (void)Cy; (void)Rad;
  fprintf(stderr, "This executable was compiled without CPLEX: "
                  "the CPLEX engine (usecplex = 1) is not available.\n");
  *status = 1;
}

void finalize_fixed_tour_solver(void)
{
}

void solve_fixed_tour(int ntour, const int* tour, double* pathcoo,
                      double* result, double* bestbound, double* gap,
                      int* bbnodes, int* niterations, double* totaltime,
                      int* status)
{
  (void)ntour; (void)tour; (void)pathcoo;
  *result = 0.0; *bestbound = 0.0; *gap = 0.0;
  *bbnodes = 0; *niterations = 0; *totaltime = 0.0;
  *status = 1;
}
