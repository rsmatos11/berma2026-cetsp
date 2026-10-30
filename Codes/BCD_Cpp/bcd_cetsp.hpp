// =====================================================================
//  bcd_cetsp: Block Coordinate Descent (BCD) for the visit points of a
//  fixed visiting sequence in the Close Enough Traveling Salesman Problem.
//
//  Given disks B(c_k, r_k), k = 0..n-1, in visiting order, BCD minimizes
//  the length of the closed tour x_0 -> x_1 -> ... -> x_{n-1} -> x_0 with
//  x_k in B(c_k, r_k). It optimizes one point at a time, in cyclic sweeps;
//  each single-disk subproblem is solved in closed form (when the segment
//  joining the two neighbors meets the disk) or by a multi-start 1D Newton
//  method on the boundary angle.
//
//  C++ translation of the Fortran code used in:
//    E. G. Birgin, P. C. R. Ertel, R. S. Matos. Using block coordinate
//    descent for solving the close enough traveling salesman problem.
// =====================================================================

#ifndef BCD_CETSP_HPP
#define BCD_CETSP_HPP

namespace bcd_cetsp {

// Choice of the point when the segment [p,q] joining the two neighbors
// meets the disk (every point of I = [p,q] cap B is then optimal):
//   Maximin    point of I farthest from both p and q (default; keeps
//              consecutive visit points apart)
//   Midpoint   midpoint of I
//   Projection point of [p,q] closest to the center
enum class TieRule { Maximin, Midpoint, Projection };

struct Options {
    int max_sweeps = 200;              // maximum number of sweeps
    double rel_tol = -1.0;             // stop when a sweep reduces the length
                                       // by a relative amount below rel_tol;
                                       // <= 0 means sqrt(machine epsilon)
    TieRule tie_rule = TieRule::Maximin;
    int window = 0;                    // > 0: replaces the rule above; stop
                                       // when the last `window` sweeps
                                       // together reduce the length by a
                                       // relative amount below window_tol
    double window_tol = 1e-6;
};

struct Result {
    double length;                     // tour length after BCD
    int sweeps;                        // sweeps performed
};

// Runs BCD on the disks (cx[k], cy[k], r[k]), k = 0..n-1, in visiting order.
// On entry, (x[k], y[k]) must lie in disk k (the centers are a valid
// choice); on exit they are the optimized visit points. Stops when a sweep
// reduces the length by a relative amount below opt.rel_tol (default
// sqrt(machine epsilon)), or after opt.max_sweeps sweeps. Coordinates are normalized internally, so
// any unit may be used. A depot is a disk with radius 0.
Result bcd(int n, const double* cx, const double* cy, const double* r,
           double* x, double* y, const Options& opt = Options());

// Single-disk subproblem: point z of B((cx,cy), R) minimizing
// |z - p| + |z - q|. Writes z to (zx, zy) and returns the minimum value.
double mindisk(double cx, double cy, double R,
               double px, double py, double qx, double qy,
               double& zx, double& zy, TieRule rule = TieRule::Maximin);

// Length of the closed tour through (x[k], y[k]), k = 0..n-1.
double tour_length(int n, const double* x, const double* y);

}  // namespace bcd_cetsp

// Plain C interface (for C, Fortran, Python, ...). tie_rule: 0 = maximin,
// 1 = midpoint, 2 = projection. Returns the tour length; *sweeps receives
// the number of sweeps (may be NULL).
extern "C" double bcd_cetsp_c(int n, const double* cx, const double* cy,
                              const double* r, double* x, double* y,
                              int max_sweeps, int tie_rule, int* sweeps);

#endif  // BCD_CETSP_HPP
