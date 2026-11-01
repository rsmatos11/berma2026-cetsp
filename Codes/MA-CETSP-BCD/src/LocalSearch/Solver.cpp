/**
 * Solver.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "LocalSearch/Solver.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "bcd_cetsp.hpp"

Solver::Solver() {}
Solver::~Solver() {}

void Solver::setContext(Centers &centers) {
    this->centers = centers;
}

// The SOCP of the visit points of a fixed sequence is solved by BCD (bcd_cetsp.cpp).

namespace {

// Number of calls and total time of BCD, printed at exit.
struct SolverStats {
    long calls = 0;
    double time = 0.0;
    long sweeps = 0;      // total sweeps
    long capped = 0;      // calls stopped by the sweep limit
    ~SolverStats() {
        std::cout << "[SOLVER] engine: BCD calls: " << calls
                  << " total_time: " << time << " s"
                  << " max_sweeps: " << (std::getenv("BCD_MAX_SWEEPS") ? std::getenv("BCD_MAX_SWEEPS") : "200")
                  << " rel_tol: " << (std::getenv("BCD_REL_TOL") ? std::getenv("BCD_REL_TOL") : "1e-6")
                  << " mean_sweeps: " << (calls > 0 ? double(sweeps) / calls : 0.0)
                  << " capped: " << capped << std::endl;
    }
} stats;

}  // namespace

static bcd_cetsp::Result runBCD(Centers &centers, List* solution,
                                const bcd_cetsp::Options &opt) {
    const int n = centers.size();
    std::vector<double> cx(n), cy(n), r(n), x(n), y(n);
    Node *p = solution->head();
    for (int i = 0; i < n; ++i) {
        cx[i] = centers[p->id][0];
        cy[i] = centers[p->id][1];
        r[i] = centers[p->id][2];
        x[i] = p->x;
        y[i] = p->y;
        // start from a feasible point: pull it onto the disk if needed
        const double d = std::hypot(x[i] - cx[i], y[i] - cy[i]);
        if (d > r[i]) {
            if (r[i] > 0.0) {
                x[i] = cx[i] + r[i] * (x[i] - cx[i]) / d;
                y[i] = cy[i] + r[i] * (y[i] - cy[i]) / d;
            } else {
                x[i] = cx[i];
                y[i] = cy[i];
            }
        }
        p = p->next;
    }

    const bcd_cetsp::Result res =
        bcd_cetsp::bcd(n, cx.data(), cy.data(), r.data(), x.data(), y.data(), opt);

    p = solution->head();
    for (int i = 0; i < n; ++i) {
        p->x = x[i];
        p->y = y[i];
        p = p->next;
    }
    solution->setValue(bcd_cetsp::tour_length(n, x.data(), y.data()));
    return res;
}

static void solveEngine(Centers &centers, List* solution) {
    // Environment: BCD_MAX_SWEEPS (default 200), BCD_REL_TOL (default 1e-6).
    static const int max_sweeps = std::getenv("BCD_MAX_SWEEPS")
        ? std::atoi(std::getenv("BCD_MAX_SWEEPS")) : 200;
    static const double rel_tol = std::getenv("BCD_REL_TOL")
        ? std::atof(std::getenv("BCD_REL_TOL")) : 1e-6;
    bcd_cetsp::Options opt;
    opt.max_sweeps = max_sweeps;
    opt.rel_tol = rel_tol;
    const bcd_cetsp::Result res = runBCD(centers, solution, opt);
    stats.sweeps += res.sweeps;
    if (res.sweeps >= opt.max_sweeps) stats.capped += 1;
}

void Solver::solve(List* solution) {
    auto t0 = std::chrono::high_resolution_clock::now();
    solveEngine(centers, solution);
    auto t1 = std::chrono::high_resolution_clock::now();
    stats.calls += 1;
    stats.time += std::chrono::duration<double>(t1 - t0).count();
}

// Final polishing of the best solution, after the search: BCD up to
// POLISH_MAX_SWEEPS sweeps (default 1000000), stopping when POLISH_WINDOW
// (default 1000) sweeps together reduce the length by a relative amount below
// POLISH_WINDOW_TOL (default 1e-6). POLISH_WINDOW=0: stop when one sweep reduces
// it by less than sqrt(machine epsilon). Not counted in the [SOLVER] statistics.
void Solver::polish(List* solution) {
    auto t0 = std::chrono::high_resolution_clock::now();
    const double before = solution->getValue();
    long sweeps = 0;
    bcd_cetsp::Options opt;
    opt.max_sweeps = std::getenv("POLISH_MAX_SWEEPS")
        ? std::atoi(std::getenv("POLISH_MAX_SWEEPS")) : 1000000;
    opt.window = std::getenv("POLISH_WINDOW")
        ? std::atoi(std::getenv("POLISH_WINDOW")) : 1000;
    opt.window_tol = std::getenv("POLISH_WINDOW_TOL")
        ? std::atof(std::getenv("POLISH_WINDOW_TOL")) : 1e-6;
    sweeps = runBCD(centers, solution, opt).sweeps;
    auto t1 = std::chrono::high_resolution_clock::now();
    const auto prec = std::cout.precision(10);
    std::cout << "[POLISH] before: " << before << " after: " << solution->getValue()
              << " sweeps: " << sweeps
              << " time: " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
    std::cout.precision(prec);
}
