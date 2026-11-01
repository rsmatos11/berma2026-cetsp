// =====================================================================
//  bcd_cetsp.cpp: BCD for the visit points of a fixed visiting sequence
//  (CETSP) and the single-disk subproblem. See bcd_cetsp.hpp.
// =====================================================================

#include "bcd_cetsp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace bcd_cetsp {

namespace {

const double PI = 3.14159265358979323846;

// Distance from the origin to the segment [(a1,b1),(a2,b2)].
double seg_dist_origin(double a1, double b1, double a2, double b2) {
    const double v1 = a2 - a1, v2 = b2 - b1;
    const double vv = v1 * v1 + v2 * v2;
    if (vv < 1e-14) {
        // Degenerate segment: distance to the nearer endpoint.
        return std::min(std::hypot(a1, b1), std::hypot(a2, b2));
    }
    double t = -(a1 * v1 + b1 * v2) / vv;
    t = std::max(0.0, std::min(1.0, t));
    return std::hypot(a1 + t * v1, b1 + t * v2);
}

// Point of segment [P1,P2] closest to C.
void closest_point_on_segment(double p1x, double p1y, double p2x, double p2y,
                              double cx, double cy, double& ox, double& oy) {
    const double vx = p2x - p1x, vy = p2y - p1y;
    const double wx = cx - p1x, wy = cy - p1y;
    const double vv = vx * vx + vy * vy;
    double t = 0.0;
    if (vv >= 1e-14) {
        t = (wx * vx + wy * vy) / vv;
        t = std::max(0.0, std::min(1.0, t));
    }
    ox = p1x + t * vx;
    oy = p1y + t * vy;
}

// f(theta) = |x(theta)-P1| + |x(theta)-P2|, x(theta) = R (cos, sin),
// disk-centered frame. sqrt is used instead of hypot (the
// coordinates are moderate here, no overflow is possible), and f' and f''
// are computed together in one Newton step (newton_step).
double f_theta(double th, double a1, double b1, double a2, double b2, double R) {
    const double c = std::cos(th), s = std::sin(th);
    const double x1 = R * c - a1, y1 = R * s - b1;
    const double x2 = R * c - a2, y2 = R * s - b2;
    return std::sqrt(x1 * x1 + y1 * y1) + std::sqrt(x2 * x2 + y2 * y2);
}

// Newton step f'(theta) / f''(theta).
double newton_step(double th, double a1, double b1, double a2, double b2, double R) {
    const double c = std::cos(th), s = std::sin(th);
    const double x1 = R * c - a1, y1 = R * s - b1;
    const double x2 = R * c - a2, y2 = R * s - b2;
    const double g1 = std::sqrt(x1 * x1 + y1 * y1);
    const double g2 = std::sqrt(x2 * x2 + y2 * y2);
    const double n1 = a1 * s - b1 * c;
    const double n2 = a2 * s - b2 * c;
    const double g1p = R * n1 / g1;
    const double g2p = R * n2 / g2;
    const double n1p = a1 * c + b1 * s;
    const double n2p = a2 * c + b2 * s;
    const double fp = g1p + g2p;
    const double fpp = (R * n1p - g1p * g1p) / g1 + (R * n2p - g2p * g2p) / g2;
    return fp / fpp;
}

// Newton's method in one variable (tol = 1e-8 on the step, at most 50
// iterations; the convergence is quadratic, so the error after the last
// step is about 1e-16).
// Returns the angle reduced to [0, 2 pi).
double newton_one(double th0, double a1, double b1, double a2, double b2, double R) {
    double th = th0;
    for (int k = 0; k < 50; ++k) {
        const double step = newton_step(th, a1, b1, a2, b2, R);
        th -= step;
        if (std::fabs(step) < 1e-8) break;
    }
    const double p = 2.0 * PI;
    return th - std::floor(th / p) * p;   // Fortran modulo(th, 2 pi)
}

// Adaptive multi-start Newton. Anchors: theta1 = direction of P1, theta2 =
// direction of P2, theta3 = direction of uhat1 + uhat2 (bisector). Newton
// starts from whichever of theta1, theta2 has the smaller f; also from the
// other one if the two values of f are close (relative gap < 1e-4), and
// from theta3 if f is nearly flat (rho = R |uhat1+uhat2| / (|P1|+|P2|) <
// 0.05). The best result is kept. Here |P1|, |P2| > R > 0.
double solve_newton(double a1, double b1, double a2, double b2, double R) {
    const double RHO_THRESHOLD = 5e-2;
    const double F_AMBIG_THRESHOLD = 1e-4;

    const double theta1 = std::atan2(b1, a1);
    const double theta2 = std::atan2(b2, a2);
    const double f1 = f_theta(theta1, a1, b1, a2, b2, R);
    const double f2 = f_theta(theta2, a1, b1, a2, b2, R);

    // primary anchor: smaller f
    double th = newton_one(f1 <= f2 ? theta1 : theta2, a1, b1, a2, b2, R);
    double fbest = f_theta(th, a1, b1, a2, b2, R);
    double theta_best = th;

    // second anchor if the two values of f are close
    if (std::fabs(f1 - f2) < F_AMBIG_THRESHOLD * std::min(f1, f2)) {
        th = newton_one(f1 <= f2 ? theta2 : theta1, a1, b1, a2, b2, R);
        const double fval = f_theta(th, a1, b1, a2, b2, R);
        if (fval < fbest) { fbest = fval; theta_best = th; }
    }

    // bisector anchor if f is nearly flat
    const double nrm1 = std::sqrt(a1 * a1 + b1 * b1);
    const double nrm2 = std::sqrt(a2 * a2 + b2 * b2);
    const double sx = a1 / nrm1 + a2 / nrm2;
    const double sy = b1 / nrm1 + b2 / nrm2;
    const double rho = R * std::sqrt(sx * sx + sy * sy) / (nrm1 + nrm2);
    if (rho < RHO_THRESHOLD) {
        th = newton_one(std::atan2(sy, sx), a1, b1, a2, b2, R);
        const double fval = f_theta(th, a1, b1, a2, b2, R);
        if (fval < fbest) { fbest = fval; theta_best = th; }
    }
    return theta_best;
}

// Interval [tlo,thi] of t such that P1 + t (P2-P1), t in [0,1], lies in
// the disk. Returns false if the segment is degenerate or the interval is
// empty.
bool segment_disk_roots(double p1x, double p1y, double p2x, double p2y,
                        double cx, double cy, double R, double& tlo, double& thi) {
    tlo = 0.0;
    thi = 0.0;
    const double vx = p2x - p1x, vy = p2y - p1y;
    const double fx = p1x - cx, fy = p1y - cy;
    const double a = vx * vx + vy * vy;
    if (a < 1e-14) return false;

    const double b = 2.0 * (fx * vx + fy * vy);
    const double c = fx * fx + fy * fy - R * R;

    // Discriminant computed from the distance of C to the line (stable).
    const double t0 = -(fx * vx + fy * vy) / a;
    const double footx = p1x + t0 * vx;
    const double footy = p1y + t0 * vy;
    const double dperp = std::hypot(footx - cx, footy - cy);
    const double disc = 4.0 * a * (R - dperp) * (R + dperp);
    if (disc < 0.0) return false;

    // Stable quadratic formula.
    const double sq = std::sqrt(disc);
    const double q = (b >= 0.0) ? -0.5 * (b + sq) : -0.5 * (b - sq);
    double t1, t2;
    if (std::fabs(q) > 1e-14) {
        t1 = q / a;
        t2 = c / q;
    } else {
        t1 = (-b - sq) / (2.0 * a);
        t2 = (-b + sq) / (2.0 * a);
    }
    tlo = std::max(0.0, std::min(t1, t2));
    thi = std::min(1.0, std::max(t1, t2));
    return tlo <= thi;
}

// Point of I = [P1,P2] cap B(C,R) chosen by the tie rule (closest point of
// the segment to C if I is empty). Maximin: the segment midpoint clamped
// to [tlo,thi].
void segment_disk_point(double p1x, double p1y, double p2x, double p2y,
                        double cx, double cy, double R, TieRule rule,
                        double& ox, double& oy) {
    if (rule == TieRule::Projection) {
        closest_point_on_segment(p1x, p1y, p2x, p2y, cx, cy, ox, oy);
        return;
    }
    double tlo, thi;
    if (!segment_disk_roots(p1x, p1y, p2x, p2y, cx, cy, R, tlo, thi)) {
        closest_point_on_segment(p1x, p1y, p2x, p2y, cx, cy, ox, oy);
        return;
    }
    const double t = (rule == TieRule::Midpoint)
                         ? 0.5 * (tlo + thi)
                         : std::max(tlo, std::min(0.5, thi));
    ox = p1x + t * (p2x - p1x);
    oy = p1y + t * (p2y - p1y);
}

}  // namespace

double mindisk(double cx, double cy, double R,
               double px, double py, double qx, double qy,
               double& zx, double& zy, TieRule rule) {
    if (R <= 0.0) {
        zx = cx;
        zy = cy;
        return std::hypot(zx - px, zy - py) + std::hypot(zx - qx, zy - qy);
    }

    // Shift to the disk-centered frame.
    const double a1 = px - cx, b1 = py - cy;
    const double a2 = qx - cx, b2 = qy - cy;
    const double dseg = seg_dist_origin(a1, b1, a2, b2);

    // Tolerance on the Case 1 / Case 2 boundary (avoids sending a neighbor
    // lying on the boundary to Newton, which would divide by zero).
    const double tol = 1e-12 * std::max({1.0, std::hypot(a1, b1), std::hypot(a2, b2)});

    if (dseg <= R + tol) {
        // ---- Case 1: segment meets the disk, closed form ----
        segment_disk_point(px, py, qx, qy, cx, cy, R, rule, zx, zy);
        double dmin = std::hypot(px - qx, py - qy);
        // If the segment passes slightly outside the disk (within the
        // tolerance), the point is pulled onto the boundary.
        const double dcen = std::hypot(zx - cx, zy - cy);
        if (dcen > R) {
            zx = cx + R * (zx - cx) / dcen;
            zy = cy + R * (zy - cy) / dcen;
            dmin = std::hypot(zx - px, zy - py) + std::hypot(zx - qx, zy - qy);
        }
        return dmin;
    }

    // ---- Case 2: boundary optimum by 1D Newton ----
    const double th = solve_newton(a1, b1, a2, b2, R);
    zx = cx + R * std::cos(th);
    zy = cy + R * std::sin(th);
    return std::hypot(zx - px, zy - py) + std::hypot(zx - qx, zy - qy);
}

double tour_length(int n, const double* x, const double* y) {
    double f = 0.0;
    for (int k = 0; k < n - 1; ++k) f += std::hypot(x[k] - x[k + 1], y[k] - y[k + 1]);
    if (n > 0) f += std::hypot(x[n - 1] - x[0], y[n - 1] - y[0]);
    return f;
}

Result bcd(int n, const double* cx, const double* cy, const double* r,
           double* x, double* y, const Options& opt) {
    Result res{0.0, 0};
    if (n <= 0) return res;

    // Normalization (as in the Fortran code): divide by the largest
    // |center coordinate| + radius.
    double s = -std::numeric_limits<double>::max();
    for (int k = 0; k < n; ++k)
        s = std::max({s, std::fabs(cx[k]) + r[k], std::fabs(cy[k]) + r[k]});
    if (s <= 0.0) s = 1.0;

    std::vector<double> ncx(n), ncy(n), nr(n), nx(n), ny(n);
    for (int k = 0; k < n; ++k) {
        ncx[k] = cx[k] / s;  ncy[k] = cy[k] / s;  nr[k] = r[k] / s;
        nx[k] = x[k] / s;    ny[k] = y[k] / s;
    }

    const double macheps12 = (opt.rel_tol > 0.0)
        ? opt.rel_tol
        : std::sqrt(std::numeric_limits<double>::epsilon());
    double cur = tour_length(n, nx.data(), ny.data());
    double prev = cur;
    int sweeps = 0;
    // window rule: hist[j % (window + 1)] = length after sweep j
    std::vector<double> hist(opt.window > 0 ? opt.window + 1 : 0);
    if (opt.window > 0) hist[0] = cur;

    while (true) {
        ++sweeps;
        for (int k = 0; k < n; ++k) {
            const int kp = (k == 0) ? n - 1 : k - 1;
            const int kn = (k == n - 1) ? 0 : k + 1;
            const double px = nx[kp], py = ny[kp], qx = nx[kn], qy = ny[kn];
            // Local cost of the two edges adjacent to point k
            const double oldloc = std::hypot(nx[k] - px, ny[k] - py) +
                                  std::hypot(nx[k] - qx, ny[k] - qy);
            const double gmin = mindisk(ncx[k], ncy[k], nr[k], px, py, qx, qy,
                                        nx[k], ny[k], opt.tie_rule);
            cur += gmin - oldloc;
        }
        if (opt.window > 0) {
            const int w = opt.window + 1;
            hist[sweeps % w] = cur;
            if (sweeps >= opt.window &&
                hist[(sweeps - opt.window) % w] - cur < opt.window_tol * cur) break;
        } else if (cur >= prev * (1.0 - macheps12)) break;
        if (sweeps >= opt.max_sweeps) break;
        prev = cur;
    }

    for (int k = 0; k < n; ++k) {
        x[k] = nx[k] * s;
        y[k] = ny[k] * s;
    }
    res.length = cur * s;
    res.sweeps = sweeps;
    return res;
}

}  // namespace bcd_cetsp

extern "C" double bcd_cetsp_c(int n, const double* cx, const double* cy,
                              const double* r, double* x, double* y,
                              int max_sweeps, int tie_rule, int* sweeps) {
    bcd_cetsp::Options opt;
    opt.max_sweeps = max_sweeps;
    opt.tie_rule = (tie_rule == 1)   ? bcd_cetsp::TieRule::Midpoint
                   : (tie_rule == 2) ? bcd_cetsp::TieRule::Projection
                                     : bcd_cetsp::TieRule::Maximin;
    const bcd_cetsp::Result res = bcd_cetsp::bcd(n, cx, cy, r, x, y, opt);
    if (sweeps) *sweeps = res.sweeps;
    return res.length;
}
