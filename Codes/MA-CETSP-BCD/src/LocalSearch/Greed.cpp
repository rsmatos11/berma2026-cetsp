/**
 * Greed.cpp
 * created on : March 03 2022
 * author : Z.LEI
 **/

#include "LocalSearch/Greed.hpp"

// Greedy pass of MA-CETSP-BCD. After LKH: up to GREEDY_BCD_SWEEPS BCD sweeps
// (default 2; 0 = original "squeeze" pass, used in the ablation study).
// GREEDY_BCD_TOL: relative tolerance (default sqrt(machine epsilon)). Number of
// calls and time are printed at exit.
#include "bcd_cetsp.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

const int GREEDY_BCD_SWEEPS = std::getenv("GREEDY_BCD_SWEEPS")
    ? std::atoi(std::getenv("GREEDY_BCD_SWEEPS")) : 2;
const double GREEDY_BCD_TOL = std::getenv("GREEDY_BCD_TOL")
    ? std::atof(std::getenv("GREEDY_BCD_TOL")) : -1.0;
// GREEDY_BCD_WHERE = after (default: only after LKH), before or both; the other
// call keeps the greedy pass below.
const std::string GREEDY_BCD_WHERE = std::getenv("GREEDY_BCD_WHERE")
    ? std::getenv("GREEDY_BCD_WHERE") : "after";
// GREEDY_BEFORE = mindisk (default: greedyMindisk below) or squeeze (original
// pass) for the greedy pass BEFORE LKH (stage 1).
// GREEDY_BEFORE_TIE = minmin (default) or maximin: choice of the point when the
// segment (a, b) meets the disk (case 3).
const std::string GREEDY_BEFORE = std::getenv("GREEDY_BEFORE")
    ? std::getenv("GREEDY_BEFORE") : "mindisk";
const std::string GREEDY_BEFORE_TIE = std::getenv("GREEDY_BEFORE_TIE")
    ? std::getenv("GREEDY_BEFORE_TIE") : "minmin";

struct GreedyStats {
    long calls = 0;
    double time = 0.0;
    long before_calls = 0;     // passes of greedyMindisk
    double before_time = 0.0;
    long case3 = 0;            // points chosen in case 3 (segment meets the disk)
    ~GreedyStats() {
        std::cout << "[GREEDY] mode: " << (GREEDY_BCD_SWEEPS > 0 ? "BCD" : "ORIGINAL")
                  << " sweeps: " << GREEDY_BCD_SWEEPS << " where: " << GREEDY_BCD_WHERE
                  << " calls: " << calls
                  << " total_time: " << time << " s"
                  << " before: " << GREEDY_BEFORE
                  << " tie: " << GREEDY_BEFORE_TIE
                  << " before_calls: " << before_calls
                  << " case3: " << case3
                  << " before_time: " << before_time << " s" << std::endl;
    }
} greedy_stats;

// Points of the segment a + t (b - a), t in [0, 1], inside the disk
// (cx, cy, r): the chord [tlo, thi]. Returns false if the segment misses it.
bool chordOfSegment(double ax, double ay, double bx, double by, double cx, double cy, double r,
                    double &tlo, double &thi) {
    const double vx = bx - ax, vy = by - ay, fx = ax - cx, fy = ay - cy;
    const double A = vx * vx + vy * vy;
    if (A < 1e-14) return false;
    const double B = 2.0 * (fx * vx + fy * vy), C = fx * fx + fy * fy - r * r;
    const double disc = B * B - 4.0 * A * C;
    if (disc < 0.0) return false;
    const double sq = std::sqrt(disc);
    tlo = std::max(0.0, (-B - sq) / (2.0 * A));
    thi = std::min(1.0, (-B + sq) / (2.0 * A));
    return tlo <= thi;
}

// Variant of the original squeeze pass, used before LKH. One
// pass in route order, depot fixed, each point updated with the previous
// point already updated (as in the original). For each node: if the previous
// point is in the disk, take it; else, if the next point is in the disk, take
// it (same test as the squeeze, Geometry::inCircle); else, if the segment
// meets the disk (case 3): with GREEDY_BEFORE_TIE=minmin, the end of the chord
// nearest to the previous or to the next point (the opposite of maximin; the
// squeeze took the end nearest to the previous point), otherwise the maximin
// point of mindisk; if the segment misses the disk (case 4), the optimal point
// on the boundary (mindisk, Newton). The new route replaces s only if shorter.
void greedyMindisk(Centers &centers, List* &s) {
    if (s->size() < 3) return;
    List* ns = new List(*s);
    double value = 0.0;
    Node* p = ns->head()->next;
    while (p != ns->head()) {
        const double cx = centers[p->id][0], cy = centers[p->id][1], r = centers[p->id][2];
        Node* pre = p->pre;
        Node* next = p->next;
        double tlo, thi;
        if (Geometry::inCircle(pre->x, pre->y, cx, cy, r)) {
            p->x = pre->x, p->y = pre->y;
        } else if (Geometry::inCircle(next->x, next->y, cx, cy, r)) {
            p->x = next->x, p->y = next->y;
        } else if (GREEDY_BEFORE_TIE == "minmin" &&
                   chordOfSegment(pre->x, pre->y, next->x, next->y, cx, cy, r, tlo, thi)) {
            // case 3, minmin: |z - pre| = tlo L at one end, |z - next| = (1 - thi) L at the other
            const double t = (tlo <= 1.0 - thi) ? tlo : thi;
            double zx = pre->x + t * (next->x - pre->x), zy = pre->y + t * (next->y - pre->y);
            const double d = std::sqrt((zx - cx) * (zx - cx) + (zy - cy) * (zy - cy));
            if (d > r) {   // rounding: pull the point onto the disk
                zx = cx + r * (zx - cx) / d;
                zy = cy + r * (zy - cy) / d;
            }
            p->x = zx, p->y = zy;
            greedy_stats.case3 += 1;
        } else {
            double zx, zy;
            bcd_cetsp::mindisk(cx, cy, r, pre->x, pre->y, next->x, next->y, zx, zy);
            p->x = zx, p->y = zy;
        }
        value += Node::distance(pre, p);
        p = p->next;
    }
    value += Node::distance(p->pre, p);
    ns->setValue(value);
    if (ns->getValue() < s->getValue()) {
        delete s;
        s = ns;
    } else {
        delete ns;
    }
}

// Up to k BCD sweeps on a copy of s; replaces s if the length decreases
// (same acceptance rule as the original greedy pass).
void greedyBCD(Centers &centers, List* &s, int k) {
    const int n = s->size();
    if (n < 3) return;
    List* ns = new List(*s);
    std::vector<double> cx(n), cy(n), r(n), x(n), y(n);
    Node *p = ns->head();
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
    bcd_cetsp::Options opt;
    opt.max_sweeps = k;
    opt.rel_tol = GREEDY_BCD_TOL;
    bcd_cetsp::bcd(n, cx.data(), cy.data(), r.data(), x.data(), y.data(), opt);
    p = ns->head();
    for (int i = 0; i < n; ++i) {
        p->x = x[i];
        p->y = y[i];
        p = p->next;
    }
    ns->setValue(bcd_cetsp::tour_length(n, x.data(), y.data()));
    if (ns->getValue() < s->getValue()) {
        delete s;
        s = ns;
    } else {
        delete ns;
    }
}

}  // namespace
// --------------------------------------------------------------------------

Greed::Greed(std::string greed_type) {
    this->greed_type = greed_type;
}
Greed::~Greed() {}

void Greed::setContext(Centers &centers) {
    this->centers = centers;
}

void Greed::run(List* &s) {
    if (GREEDY_BEFORE == "mindisk" && stage == 1) {
        auto t0 = std::chrono::high_resolution_clock::now();
        greedyMindisk(centers, s);
        auto t1 = std::chrono::high_resolution_clock::now();
        greedy_stats.before_calls += 1;
        greedy_stats.before_time += std::chrono::duration<double>(t1 - t0).count();
        return;
    }
    if (GREEDY_BCD_SWEEPS > 0 &&
        (GREEDY_BCD_WHERE == "both" ||
         (GREEDY_BCD_WHERE == "after" && stage == 2) ||
         (GREEDY_BCD_WHERE == "before" && stage == 1))) {
        auto t0 = std::chrono::high_resolution_clock::now();
        greedyBCD(centers, s, GREEDY_BCD_SWEEPS);
        auto t1 = std::chrono::high_resolution_clock::now();
        greedy_stats.calls += 1;
        greedy_stats.time += std::chrono::duration<double>(t1 - t0).count();
        return;
    }
    auto start = std::chrono::high_resolution_clock::now();

    List* ns = new List(*s);
    // size must be equal to or greater than 3
    if (ns->size() < 3) return;
    double value = 0;
    Node* p = ns->head()->next;
    while (p != ns->head()) {
        int id = p->id;
        value += updatePosition(p, centers[id][0], centers[id][1], centers[id][2]);
        p = p->next;
    }
    value += Node::distance(p->pre, p);
    ns->setValue(value);
    if (ns->getValue() < s->getValue()) {
        delete s;
        s = ns;
    } else if (LOG){
        std::cout << "[GREED] not accept" << std::endl;
    };

    auto end = std::chrono::high_resolution_clock::now();
    if (LOG) std::cout << "greed solution : " << s->getValue() << " time : " << std::chrono::duration<double> (end - start).count() << " s" << std::endl;
}

double Greed::updatePosition(Node* node, double x0, double y0, double r) {
    Node* pre = node->pre;
    Node* next = node->next;
    double x1 = pre->x, y1 = pre->y, x2 = next->x, y2 = next->y;
    // judge points in or out circle
    bool in_circle1 = Geometry::inCircle(x1, y1, x0, y0, r);
    bool in_circle2 = Geometry::inCircle(x2, y2, x0, y0, r);
    if (in_circle1 && in_circle2) {
        if (greed_type == "SQUEEZE") {
            node->x = x1, node->y = y1;
            return 0;
        } else if (greed_type == "SPARSE") {
            node->x = (x1 + x2) / 2, node->y = (y1 + y2) / 2;
            return Node::distance(pre, node);
        }
    } else if (in_circle1) {
        node->x = x1, node->y = y1;
        return 0;
    } else if (in_circle2) {
        node->x = x2, node->y = y2;
        return Node::distance(pre, node);
    } else {
        auto intersections = Geometry::solveLineIntersectSphere(x1, y1, x2, y2, x0, y0, r);
        if (intersections.size() == 2) {
            node->x = intersections[0][0];
            node->y = intersections[0][1];
        } else if (intersections.size() == 1){
            node->x = intersections[0][0];
            node->y = intersections[0][1];
        } else {
            AlhazenProblem ap(x1, y1, x2, y2, x0, y0, r);
            auto position = ap.solve();
            node->x = position[0];
            node->y = position[1];
        }
        return Node::distance(pre, node);
    }
    return Node::distance(pre, node);
}

std::vector<double> Greed::approxPosition(double x0, double y0, double r, Node* pre, Node* next) {
    double x, y;
    double x1 = pre->x, y1 = pre->y, x2 = next->x, y2 = next->y;
    // judge points in or out circle
    bool in_circle1 = Geometry::inCircle(x1, y1, x0, y0, r);
    bool in_circle2 = Geometry::inCircle(x2, y2, x0, y0, r);
    if (in_circle1 && in_circle2) {
        x = x1, y = y1;
        if (greed_type == "SQUEEZE") {
            x = x1, y = y1;
        } else if (greed_type == "SPARSE") {
            x = (x1 + x2) / 2, y = (y1 + y2) / 2;
        }
    } else if (in_circle1) {
        x = x1, y = y1;
    } else if (in_circle2) {
        x = x2, y = y2;
    } else {
        auto intersections = Geometry::solveLineIntersectSphere(x1, y1, x2, y2, x0, y0, r);
        if (intersections.size() == 2) {
            x = (intersections[0][0] + intersections[1][0]) / 2;
            y = (intersections[0][1] + intersections[1][1]) / 2;
        } else if (intersections.size() == 1){
            x = intersections[0][0];
            y = intersections[0][1];
        } else {
            double mid_point_x = (pre->x + next->x) / 2;
            double mid_point_y = (pre->y + next->y) / 2;
            auto intersections = Geometry::solveLineIntersectSphere(mid_point_x, mid_point_y, x0, y0, x0, y0, r);
            if (intersections.size() != 1) std::cout << "ERROR : only one intersection !!!" << std::endl;
            x = intersections[0][0];
            y = intersections[0][1];
        }
    }
    return std::vector<double> {x, y};
}

bool Greed::inLine(double x0, double y0, double r, Node* pre, Node* next) {
    double x1 = pre->x, y1 = pre->y, x2 = next->x, y2 = next->y;
    // in circle
    bool in_circle1 = Geometry::inCircle(x1, y1, x0, y0, r);
    bool in_circle2 = Geometry::inCircle(x2, y2, x0, y0, r);
    if (in_circle1 || in_circle2) return true;
    // have intersections
    auto intersections = Geometry::solveLineIntersectSphere(x1, y1, x2, y2, x0, y0, r);
    if (intersections.size() > 0) return true;
    return false;
}
