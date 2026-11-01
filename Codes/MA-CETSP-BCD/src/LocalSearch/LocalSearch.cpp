/**
 * LocalSearch.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "LocalSearch/LocalSearch.hpp"

// --- Exact single-disk subproblem in the local search ---------------------
// In bestImproveApproxRelocate and bestImproveApproxSwap, the position of a
// moved node and the cost of the move use bcd_cetsp::mindisk (optimal point of
// the disk for the two new neighbors, maximin tie rule) instead of the
// original approximations (projection of the midpoint, middle of the chord,
// squeeze). LS_MINDISK=0 restores the original code. Calls are printed at exit.
// Before calling mindisk, a cheap lower bound of the insertion
// cost discards candidates that cannot beat the best move found so far
// (same moves, fewer mindisk calls). LS_PRUNE=0 turns it off.
#include "bcd_cetsp.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace {

const bool LS_MINDISK = !std::getenv("LS_MINDISK") || std::atoi(std::getenv("LS_MINDISK")) != 0;
const bool LS_PRUNE = !std::getenv("LS_PRUNE") || std::atoi(std::getenv("LS_PRUNE")) != 0;

struct MindiskStats {
    long calls = 0;
    long pruned = 0;      // candidates discarded by the lower bound
    ~MindiskStats() {
        std::cout << "[LS] mindisk: " << (LS_MINDISK ? "on" : "off")
                  << " calls: " << calls
                  << " prune: " << (LS_PRUNE ? "on" : "off")
                  << " pruned: " << pruned << std::endl;
    }
} mindisk_stats;

// Lower bound of |z - pre| + |z - next| over the disk c, the
// largest of (i) |pre - next|; (ii) the distances of pre and next to the
// disk, summed; (iii) 2 sqrt((d/2)^2 + h^2), the smallest sum at distance h
// from the line pre-next, h = distance of the disk to that line.
double lowerBound(const std::vector<double> &c, Node* pre, Node* next) {
    const double ax = pre->x, ay = pre->y, bx = next->x, by = next->y;
    const double ux = bx - ax, uy = by - ay;
    const double d = std::sqrt(ux * ux + uy * uy);
    double lb = std::max(0.0, std::sqrt((ax - c[0]) * (ax - c[0]) + (ay - c[1]) * (ay - c[1])) - c[2])
              + std::max(0.0, std::sqrt((bx - c[0]) * (bx - c[0]) + (by - c[1]) * (by - c[1])) - c[2]);
    if (d > 0.0) {
        const double h = std::fabs(ux * (c[1] - ay) - uy * (c[0] - ax)) / d - c[2];
        if (h > 0.0) lb = std::max(lb, 2.0 * std::sqrt(0.25 * d * d + h * h));
    }
    return std::max(lb, d);
}

// Discard only when the bound is clearly not better (rounding margin), so
// that the chosen moves are those without the bound.
inline bool cannotImprove(double lb_delta, double best_delta, double scale) {
    return lb_delta >= best_delta + 1e-10 * (1.0 + scale);
}

// Optimal point (zx, zy) of disk c for neighbors pre and next; returns
// |z - pre| + |z - next|.
double exactPosition(const std::vector<double> &c, Node* pre, Node* next, double &zx, double &zy) {
    mindisk_stats.calls += 1;
    return bcd_cetsp::mindisk(c[0], c[1], c[2], pre->x, pre->y, next->x, next->y, zx, zy);
}

}  // namespace
// --------------------------------------------------------------------------

// LocalSearch::LocalSearch(std::string improvement, std::string greed_type) : greed(greed_type) {
//     this->improvement = improvement;
// }

LocalSearch::LocalSearch(Parameters *params) : greed(params->greed) {
    this->lkh_random_num = params->random_num;
    this->improvement = params->improvement;
}

LocalSearch::~LocalSearch() {}

void LocalSearch::setContext(Random *random, Centers &centers, Neighbor* neighbor, std::string timestamp) {
    this->random = random;
    this->centers = centers;
    this->neighbor = neighbor;
    solver.setContext(centers);
    if (lkh_random_num == 0) {
        lkh_random_num = random->randomInt(INT_MAX);
    }
    lkh.setContext(timestamp, lkh_random_num);
    greed.setContext(centers);
}


// Stage times (see LocalSearch.hpp).
StageTimes stage_times;

static double secondsSince(std::chrono::high_resolution_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
}

List *LocalSearch::initSolOpt(List *s) {
    auto t0 = std::chrono::high_resolution_clock::now();
    s = lkh.run(s, false);
    stage_times.lkh += secondsSince(t0);
    t0 = std::chrono::high_resolution_clock::now();
    greed.stage = 2;   // after LKH
    greed.run(s);
    stage_times.greedy += secondsSince(t0);
    return s;
}


List *LocalSearch::VND(List *s) {
    s->evaluate();
    if (LOG) std::cout << "offspring solution : " << s->getValue() << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    greed.stage = 1;   // before LKH
    greed.run(s);
    stage_times.greedy += secondsSince(t0);
    t0 = std::chrono::high_resolution_clock::now();
    s = lkh.run(s, true);
    stage_times.lkh += secondsSince(t0);
    t0 = std::chrono::high_resolution_clock::now();
    greed.stage = 2;   // after LKH
    greed.run(s);
    stage_times.greedy += secondsSince(t0);
    t0 = std::chrono::high_resolution_clock::now();
    jointOpt(s);
    stage_times.ls += secondsSince(t0);
    t0 = std::chrono::high_resolution_clock::now();
    solver.solve(s);
    stage_times.socp += secondsSince(t0);
    stage_times.vnd_calls += 1;
    return s;
}

void LocalSearch::jointOpt(List* s) {
    this->neighbors = this->neighbor->getNeighbors();
    auto start = std::chrono::high_resolution_clock::now();
    if (improvement == "FIRST") {
        firstImproveApproxRelocate(s);
        firstImproveApproxSwap(s);
    } else if (improvement == "BEST") {
        bestImproveApproxRelocate(s);
        bestImproveApproxSwap(s);
    } else {
        std::cerr << "[ERROR] unknown improvement type" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    if (LOG) std::cout << "joint solution : " << s->getValue() << " time : " << std::chrono::duration<double> (end - start).count() << " s" << std::endl;
}

void LocalSearch::randomTwoOpt(List* s) {
    if (s->size() < 4) return;
    auto start = std::chrono::high_resolution_clock::now();
    // first-improvement 2-opt with random order
    int size = s->size();
    std::vector<Node*> nodes(size);
    Node* p = s->head();
    for (int i = 0; i < size; ++i) {
        nodes[i] = p;
        p = p->next;
    }
    bool improved = true;
    while (improved) {
        improved = false;
        random->permutation(nodes);
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (nodes[i]->next != nodes[j] && nodes[j]->next != nodes[i] && nodes[i]->next->next != nodes[j] && nodes[j]->next->next != nodes[i]) {
                    Node* pi = nodes[i];
                    Node* pj = nodes[j];
                    Node* pi_pre = pi->pre;
                    Node* pj_pre = pj->pre;
                    double delta = Node::distance(pi_pre, pj_pre) + Node::distance(pi, pj) - Node::distance(pi_pre, pi) - Node::distance(pj_pre, pj);
                    if (delta < -EPSILON) {
                        // reverse i->j-1
                        s->reverse(pi, pj);
                        pi_pre->next = pj_pre;
                        pj_pre->pre = pi_pre;
                        pi->next = pj;
                        pj->pre = pi;
                        s->setValue(s->getValue() + delta);
                        improved = true;
                    }
                }
                if (improved) break;
            }
            if (improved) break;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    if (LOG) std::cout << "2-opt solution: " << s->getValue() << " time : " << std::chrono::duration<double> (end - start).count() << " s" << std::endl;
}

void LocalSearch::bestImproveApproxRelocate(List *s) {
    int size = s->size();
    std::vector<Node*> nodes(size);
    Node* p = s->head();
    for (int i = 0; i < size; ++i) {
        nodes[i] = p;
        p = p->next;
    }
    bool improved = true;
    while (improved) {
        improved = false;
        random->permutation(nodes);

        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i]->id == 0) continue;
            int best_j;
            double node_x = 0, node_y = 0;
            double best_delta = 0;
            bool best_is_in_line = true;
            for (int j = 0; j < nodes.size(); ++j) {
                if (i == j || nodes[i] == nodes[j]->next) continue;
                if (!neighbors[nodes[i]->id][nodes[j]->id]) continue;        // neighbors or global
                // move nodes[i] after nodes[j]
                double delta = 0;
                double approx_point_x = 0, approx_point_y = 0;
                bool is_in_line = true;
                if (LS_MINDISK) {
                    // skip candidates whose lower bound cannot beat best_delta
                    if (LS_PRUNE) {
                        const double dpq = Node::distance(nodes[j], nodes[j]->next);
                        const double lb = lowerBound(centers[nodes[i]->id], nodes[j], nodes[j]->next);
                        const double lb_delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                                + lb - dpq;
                        if (cannotImprove(lb_delta, best_delta, lb)) {
                            mindisk_stats.pruned += 1;
                            continue;
                        }
                    }
                    // exact insertion cost and point (both cases)
                    const double ins = exactPosition(centers[nodes[i]->id], nodes[j], nodes[j]->next, approx_point_x, approx_point_y);
                    delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                            + ins - Node::distance(nodes[j], nodes[j]->next);
                    is_in_line = false;
                } else if (greed.inLine(centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2], nodes[j], nodes[j]->next)) {
                    delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next);
                    is_in_line = true;
                } else {
                    double mid_point_x = (nodes[j]->x + nodes[j]->next->x) / 2;
                    double mid_point_y = (nodes[j]->y + nodes[j]->next->y) / 2;
                    auto intersections = Geometry::solveLineIntersectSphere(mid_point_x, mid_point_y, centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2]);
                    if (intersections.size() != 1) {
                        std::cout << "ERROR : " << intersections.size() << " intersections" << std::endl;
                    }
                    approx_point_x = intersections[0][0];
                    approx_point_y = intersections[0][1];
                    delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                            + Geometry::EucDistance(approx_point_x, approx_point_y, nodes[j]->x, nodes[j]->y) + Geometry::EucDistance(approx_point_x, approx_point_y, nodes[j]->next->x, nodes[j]->next->y)
                            - Node::distance(nodes[j], nodes[j]->next);
                    is_in_line = false;
                }
                if (delta < best_delta) {
                    best_delta = delta;
                    best_j = j;
                    node_x = approx_point_x;
                    node_y = approx_point_y;
                    best_is_in_line = is_in_line;
                }
            }

            if (best_delta < 0) {
                if (best_delta < -EPSILON) improved = true;
                else if (random->randomInt(100) < 50) continue;
                Node* pi = nodes[i];
                Node* pj = nodes[best_j];
                Node* pi_pre = pi->pre;
                Node* pi_next = pi->next;
                Node* pj_next = pj->next;
                pi_pre->next = pi_next;
                pi_next->pre = pi_pre;
                pi->pre = pj;
                pi->next = pj_next;
                pj->next = pi;
                pj_next->pre = pi;
                if (best_is_in_line) {
                    greed.updatePosition(pi, centers[pi->id][0], centers[pi->id][1], centers[pi->id][2]);
                } else {
                    pi->x = node_x;
                    pi->y = node_y;
                }
                s->setValue(s->getValue() + best_delta);
            }
        }
    }
}

void LocalSearch::bestImproveApproxSwap(List* s) {
    int size = s->size();
    std::vector<Node*> nodes(size);
    Node* p = s->head();
    for (int i = 0; i < size; ++i) {
        nodes[i] = p;
        p = p->next;
    }
    bool improved = true;
    while (improved) {
        improved = false;
        random->permutation(nodes);
        for (int i = 0; i < size; ++i) {
            if (nodes[i]->id == 0) continue;

            int best_j;
            double nix, niy, njx, njy;
            double best_delta = 0;

            for (int j = i + 1; j < size; ++j) {
                if (nodes[j]->id == 0 || nodes[i]->next == nodes[j] || nodes[j]->next == nodes[i]) continue;
                if (!neighbors[nodes[i]->id][nodes[j]->id]) continue;
                // swap nodes[i] and nodes[j]

                std::vector<double> posi(2), posj(2);
                if (LS_MINDISK) {
                    // skip pairs whose lower bound cannot beat best_delta
                    if (LS_PRUNE) {
                        const double lb = lowerBound(centers[nodes[i]->id], nodes[j]->pre, nodes[j]->next)
                                        + lowerBound(centers[nodes[j]->id], nodes[i]->pre, nodes[i]->next);
                        const double lb_delta = lb
                                - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                                - Node::distance(nodes[j], nodes[j]->pre) - Node::distance(nodes[j], nodes[j]->next);
                        if (cannotImprove(lb_delta, best_delta, lb)) {
                            mindisk_stats.pruned += 1;
                            continue;
                        }
                    }
                    // exact points (i and j are not adjacent, so the two are independent)
                    exactPosition(centers[nodes[i]->id], nodes[j]->pre, nodes[j]->next, posi[0], posi[1]);
                    exactPosition(centers[nodes[j]->id], nodes[i]->pre, nodes[i]->next, posj[0], posj[1]);
                } else {
                    posi = greed.approxPosition(centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2], nodes[j]->pre, nodes[j]->next);
                    posj = greed.approxPosition(centers[nodes[j]->id][0], centers[nodes[j]->id][1], centers[nodes[j]->id][2], nodes[i]->pre, nodes[i]->next);
                }

                double delta = Geometry::EucDistance(posi[0], posi[1], nodes[j]->pre->x, nodes[j]->pre->y) + Geometry::EucDistance(posi[0], posi[1], nodes[j]->next->x, nodes[j]->next->y)
                               + Geometry::EucDistance(posj[0], posj[1], nodes[i]->pre->x, nodes[i]->pre->y) + Geometry::EucDistance(posj[0], posj[1], nodes[i]->next->x, nodes[i]->next->y)
                               - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                               - Node::distance(nodes[j], nodes[j]->pre) - Node::distance(nodes[j], nodes[j]->next);

                if (delta < best_delta) {
                    best_delta = delta;
                    best_j = j;
                    nix = posi[0];
                    niy = posi[1];
                    njx = posj[0];
                    njy = posj[1];
                }
            }
            if (best_delta < 0) {
                if (best_delta < -EPSILON) improved = true;
                else if (random->randomInt(100) < 50) continue;
                Node *pi = nodes[i];
                Node *pj = nodes[best_j];
                Node *pi_pre = pi->pre;
                Node *pj_pre = pj->pre;
                Node *pi_next = pi->next;
                Node *pj_next = pj->next;
                pi_pre->next = pj;
                pj->pre = pi_pre;
                pi_next->pre = pj;
                pj->next = pi_next;
                pj_pre->next = pi;
                pi->pre = pj_pre;
                pj_next->pre = pi;
                pi->next = pj_next;
                pi->x = nix, pi->y = niy;
                pj->x = njx, pj->y = njy;
                s->setValue(s->getValue() + best_delta);
            }
        }
    }

}

void LocalSearch::firstImproveApproxRelocate(List* s) {
    int size = s->size();
    std::vector<Node*> nodes(size);
    Node* p = s->head();
    for (int i = 0; i < size; ++i) {
        nodes[i] = p;
        p = p->next;
    }
    bool improved = true;
    while (improved) {
        improved = false;
        random->permutation(nodes);
        bool isInLine = true;
        for (int i = 0; i < size; ++i) {
            if (nodes[i]->id == 0) continue;
            for (int j = 0; j < size; ++j) {
                if (i == j || nodes[i] == nodes[j]->next) continue;
                if (!neighbors[nodes[i]->id][nodes[j]->id]) continue;        // neighbors or global
                // move nodes[i] after nodes[j]
                double delta;
                double approx_point_x ,approx_point_y;
                if (greed.inLine(centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2], nodes[j], nodes[j]->next)) {
                    delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next);
                    isInLine = true;
                } else {
                    double mid_point_x = (nodes[j]->x + nodes[j]->next->x) / 2;
                    double mid_point_y = (nodes[j]->y + nodes[j]->next->y) / 2;
                    auto intersections = Geometry::solveLineIntersectSphere(mid_point_x, mid_point_y, centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2]);
                    if (intersections.size() != 1) {
                        std::cout << "ERROR : " << intersections.size() << " intersections" << std::endl;
                    }
                    approx_point_x = intersections[0][0];
                    approx_point_y = intersections[0][1];
                    delta = Node::distance(nodes[i]->pre, nodes[i]->next) - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                            + Geometry::EucDistance(approx_point_x, approx_point_y, nodes[j]->x, nodes[j]->y) + Geometry::EucDistance(approx_point_x, approx_point_y, nodes[j]->next->x, nodes[j]->next->y)
                            - Node::distance(nodes[j], nodes[j]->next);
                    isInLine = false;
                }
                if (delta < -EPSILON) {
                    Node* pi = nodes[i];
                    Node* pj = nodes[j];
                    Node* pi_pre = pi->pre;
                    Node* pi_next = pi->next;
                    Node* pj_next = pj->next;
                    pi_pre->next = pi_next;
                    pi_next->pre = pi_pre;
                    pi->pre = pj;
                    pi->next = pj_next;
                    pj->next = pi;
                    pj_next->pre = pi;
                    if (isInLine) {
                        greed.updatePosition(pi, centers[pi->id][0], centers[pi->id][1], centers[pi->id][2]);
                    } else {
                        pi->x = approx_point_x;
                        pi->y = approx_point_y;
                    }

                    s->setValue(s->getValue() + delta);
                    improved = true;
                }
                if (improved) break;
            }
            if (improved) break;
        }
    }
}

void LocalSearch::firstImproveApproxSwap(List* s) {
    int size = s->size();
    std::vector<Node*> nodes(size);
    Node* p = s->head();
    for (int i = 0; i < size; ++i) {
        nodes[i] = p;
        p = p->next;
    }
    bool improved = true;
    while (improved) {
        improved = false;
        random->permutation(nodes);
        for (int i = 0; i < size; ++i) {
            if (nodes[i]->id == 0) continue;
            for (int j = i + 1; j < size; ++j) {
                if (nodes[j]->id == 0) continue;
                // swap nodes[i] and nodes[j]
                if (nodes[i]->next == nodes[j] || nodes[j]->next == nodes[i]) continue;
                if (!neighbors[nodes[i]->id][nodes[j]->id]) continue;
                auto posi = greed.approxPosition(centers[nodes[i]->id][0], centers[nodes[i]->id][1], centers[nodes[i]->id][2], nodes[j]->pre, nodes[j]->next);
                auto posj = greed.approxPosition(centers[nodes[j]->id][0], centers[nodes[j]->id][1], centers[nodes[j]->id][2], nodes[i]->pre, nodes[i]->next);
                double delta = Geometry::EucDistance(posi[0], posi[1], nodes[j]->pre->x, nodes[j]->pre->y) + Geometry::EucDistance(posi[0], posi[1], nodes[j]->next->x, nodes[j]->next->y)
                        + Geometry::EucDistance(posj[0], posj[1], nodes[i]->pre->x, nodes[i]->pre->y) + Geometry::EucDistance(posj[0], posj[1], nodes[i]->next->x, nodes[i]->next->y)
                        - Node::distance(nodes[i], nodes[i]->pre) - Node::distance(nodes[i], nodes[i]->next)
                        - Node::distance(nodes[j], nodes[j]->pre) - Node::distance(nodes[j], nodes[j]->next);

                if (delta < -EPSILON) {
                    Node *pi = nodes[i];
                    Node *pj = nodes[j];
                    Node *pi_pre = pi->pre;
                    Node *pj_pre = pj->pre;
                    Node *pi_next = pi->next;
                    Node *pj_next = pj->next;
                    pi_pre->next = pj;
                    pj->pre = pi_pre;
                    pi_next->pre = pj;
                    pj->next = pi_next;
                    pj_pre->next = pi;
                    pi->pre = pj_pre;
                    pj_next->pre = pi;
                    pi->next = pj_next;
                    pi->x = posi[0];
                    pi->y = posi[1];
                    pj->x = posj[0];
                    pj->y = posj[1];
                    s->setValue(s->getValue() + delta);
                    improved = true;
                }
                if (improved) break;
            }
            if (improved) break;
        }
    }
}
