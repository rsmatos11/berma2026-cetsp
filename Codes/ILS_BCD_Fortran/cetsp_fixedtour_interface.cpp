#include <ilcplex/ilocplex.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

//============================================================================//
//   CPLEX solver for the visit points of a fixed tour, called from Fortran:
//
//       min  sum_k d_k
//       s.t. (yx_k - Cx_k)^2 + (yy_k - Cy_k)^2 <= Rad_k^2         for all k
//            (yx_k - yx_{k+1})^2 + (yy_k - yy_{k+1})^2 <= d_k^2   cyclic
//
//   (the same problem solved by BCD). The disk data is copied once by
//   initialize_fixed_tour_solver; each call to solve_fixed_tour only
//   receives the visiting order.
//============================================================================//

namespace {

  bool persistentReady = false;

  std::vector<double> persistentCx;
  std::vector<double> persistentCy;
  std::vector<double> persistentRad;

  double persistentXmin = 0.0;
  double persistentXmax = 0.0;
  double persistentYmin = 0.0;
  double persistentYmax = 0.0;

  // ---- optional call-timing accounting ----
  const bool debugTiming = false;

  double totalBuildTime = 0.0;
  double totalSolveTime = 0.0;
  int    totalCalls     = 0;

} // unnamed namespace

extern "C" {

  // Copies (n,Cx,Cy,Rad) into the solver (called once). status: 0 = ok,
  // 1 = invalid n.
  void initialize_fixed_tour_solver(int n, const double* Cx, const double* Cy,
                                     const double* Rad, int* status)
  {
    *status = 1;

    if ( n <= 0 ) {
      persistentReady = false;
      return;
    }

    persistentCx.assign(Cx, Cx + n);
    persistentCy.assign(Cy, Cy + n);
    persistentRad.assign(Rad, Rad + n);

    persistentXmin = persistentCx[0] - persistentRad[0];
    persistentXmax = persistentCx[0] + persistentRad[0];
    persistentYmin = persistentCy[0] - persistentRad[0];
    persistentYmax = persistentCy[0] + persistentRad[0];

    for ( int i = 1; i < n; ++i ) {
      persistentXmin = std::min(persistentXmin,persistentCx[i] - persistentRad[i]);
      persistentXmax = std::max(persistentXmax,persistentCx[i] + persistentRad[i]);
      persistentYmin = std::min(persistentYmin,persistentCy[i] - persistentRad[i]);
      persistentYmax = std::max(persistentYmax,persistentCy[i] + persistentRad[i]);
    }

    persistentReady = true;

    *status = 0;
  }

  // Releases the solver (called once).
  void finalize_fixed_tour_solver()
  {
    persistentCx.clear();
    persistentCy.clear();
    persistentRad.clear();
    persistentReady = false;

    if ( debugTiming ) {
      std::cerr << "==================================================" << std::endl;
      std::cerr << "[CPLEX TIMING] Total calls: "          << totalCalls << std::endl;
      std::cerr << "[CPLEX TIMING] Total model-build time: " << totalBuildTime << " s" << std::endl;
      std::cerr << "[CPLEX TIMING] Total solve time:        " << totalSolveTime << " s" << std::endl;
      if ( totalCalls > 0 ) {
        std::cerr << "[CPLEX TIMING] Avg build per call: " << (totalBuildTime / totalCalls) * 1000.0 << " ms" << std::endl;
        std::cerr << "[CPLEX TIMING] Avg solve per call: " << (totalSolveTime / totalCalls) * 1000.0 << " ms" << std::endl;
      }
      std::cerr << "==================================================" << std::endl;
    }
  }

  // Solves the problem for the visiting order tour(1:ntour) (1-based, ntour
  // = n). pathcoo(2*k, 2*k+1) receives the visit point of tour[k].
  // status: 0 = ok, 2 = CPLEX could not solve, 3 = IloException,
  //         4 = std::exception, 5 = unknown exception,
  //         6 = solver not initialized or ntour != n.
  void solve_fixed_tour(int ntour, const int* tour, double* pathcoo,
                         double* result, double* bestbound, double* gap,
                         int* bbnodes, int* niterations, double* totaltime,
                         int* status)
  {
    *result      = 0.0;
    *bestbound   = 0.0;
    *gap         = 0.0;
    *bbnodes     = 0;
    *niterations = 0;
    *totaltime   = 0.0;
    *status      = 1;

    const auto totalStart = std::chrono::steady_clock::now();

    if ( !persistentReady || ntour != static_cast<int>(persistentCx.size()) ) {
      std::cerr << "The persistent CETSP CPLEX instance has not been initialized "
                << "or ntour does not match n." << std::endl;
      *status = 6;
      return;
    }

    IloEnv env;

    try {

      const auto buildStart = std::chrono::steady_clock::now();

      IloModel model(env);

      IloNumVarArray d(env,ntour,0.0,IloInfinity);
      IloNumVarArray yx(env,ntour,persistentXmin,persistentXmax);
      IloNumVarArray yy(env,ntour,persistentYmin,persistentYmax);

      IloExpr obj(env);
      for ( int k = 0; k < ntour; ++k ) {
        obj += d[k];
      }
      model.add(IloMinimize(env,obj));
      obj.end();

      // ---- disk constraints ----
      // A zero-radius disk (the depot) is fixed at its center by variable
      // bounds (a quadratic constraint would have no interior point and
      // stall the barrier method).
      const double radius_eps = 1.0e-12;

      for ( int k = 0; k < ntour; ++k ) {

        const int disk = tour[k] - 1;

        if ( persistentRad[disk] <= radius_eps ) {

          yx[k].setBounds(persistentCx[disk],persistentCx[disk]);
          yy[k].setBounds(persistentCy[disk],persistentCy[disk]);

          continue;
        }

        IloExpr expr = (yx[k] - persistentCx[disk]) * (yx[k] - persistentCx[disk]) +
                       (yy[k] - persistentCy[disk]) * (yy[k] - persistentCy[disk]);

        model.add(expr <= persistentRad[disk] * persistentRad[disk]);

        expr.end();
      }

      // ---- distance (arc-length) constraints, cyclic ----
      for ( int k = 0; k < ntour; ++k ) {

        const int knext = (k + 1 == ntour) ? 0 : k + 1;

        IloExpr lhs = (yx[k] - yx[knext]) * (yx[k] - yx[knext]) +
                      (yy[k] - yy[knext]) * (yy[k] - yy[knext]);

        model.add(lhs <= d[k] * d[k]);

        lhs.end();
      }

      IloCplex cplex(model);

      cplex.setOut(env.getNullStream());
      cplex.setWarning(env.getNullStream());
      cplex.setParam(IloCplex::Param::Threads,1);

      const auto buildEnd = std::chrono::steady_clock::now();

      const bool solved = cplex.solve();

      const auto solveEnd = std::chrono::steady_clock::now();

      totalBuildTime += std::chrono::duration<double>(buildEnd - buildStart).count();
      totalSolveTime += std::chrono::duration<double>(solveEnd - buildEnd).count();
      totalCalls++;

      if ( !solved ) {

        std::cerr << "CPLEX could not solve the fixed-tour CETSP model." << std::endl;
        *status = 2;

      } else {

        *result      = cplex.getObjValue();
        *bestbound   = *result;
        *gap         = 0.0;
        *bbnodes     = 0;
        *niterations = cplex.getNiterations();

        for ( int k = 0; k < ntour; ++k ) {
          pathcoo[2 * k]     = cplex.getValue(yx[k]);
          pathcoo[2 * k + 1] = cplex.getValue(yy[k]);
        }

        *status = 0;
      }

    } catch ( const IloException& error ) {

      std::cerr << "CPLEX exception: " << error.getMessage() << std::endl;
      *status = 3;

    } catch ( const std::exception& error ) {

      std::cerr << "C++ exception: " << error.what() << std::endl;
      *status = 4;

    } catch ( ... ) {

      std::cerr << "Unknown exception in solve_fixed_tour." << std::endl;
      *status = 5;
    }

    env.end();

    const auto totalEnd = std::chrono::steady_clock::now();

    *totaltime = std::chrono::duration<double>(totalEnd - totalStart).count();
  }

} // extern "C"
