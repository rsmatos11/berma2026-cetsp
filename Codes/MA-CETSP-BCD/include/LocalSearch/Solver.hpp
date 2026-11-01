/**
 * Solver.hpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#ifndef CETSP_SOLVER_HPP
#define CETSP_SOLVER_HPP

#include "Defs.hpp"
#include "Genetic/List.hpp"
#include <chrono>

class Solver {
private:
    std::vector<std::vector<double>> centers;    // (x, y, r)
public:
    Solver();
    ~Solver();
    void setContext(Centers &centers);
    void solve(List* solution);
    void polish(List* solution);  // final polishing of the best solution
};

#endif // CETSP_SOLVER_HPP