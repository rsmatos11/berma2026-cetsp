// =====================================================================
//  bcd_rota: example program. Reads an instance and a visiting sequence,
//  runs BCD (up to 200 sweeps) and writes the visit points.
//
//  Usage:
//    bcd_rota instance.dat withdepot sequence.txt points_in.txt points_out.txt
//
//  instance.dat    n; depot "x y" (only if withdepot = 1); one line
//                  "x y r" per disk
//  sequence.txt    one node per line, numbered 1..ntot (disks 1..n and, if
//                  withdepot = 1, depot = n+1); any rotation
//  points_in.txt   initial point of each node ("x y", in node order
//                  1..ntot); must be feasible (e.g. the centers)
//  points_out.txt  visit points after BCD, same format
//
//  Prints the tour length and the number of sweeps. Same interface as the
//  Fortran program of the same name, for comparison.
// =====================================================================

#include "bcd_cetsp.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: bcd_rota instance.dat withdepot sequence.txt "
                     "points_in.txt points_out.txt\n";
        return 1;
    }
    const int withdepot = std::atoi(argv[2]);

    // ---- instance (nodes 0..ntot-1; depot, if any, is the last node) ----
    std::ifstream fi(argv[1]);
    if (!fi) { std::cerr << "Cannot open " << argv[1] << "\n"; return 1; }
    int nballs;
    fi >> nballs;
    const int ntot = withdepot ? nballs + 1 : nballs;
    std::vector<double> cx(ntot), cy(ntot), r(ntot);
    if (withdepot) {
        fi >> cx[ntot - 1] >> cy[ntot - 1];
        r[ntot - 1] = 0.0;
    }
    for (int i = 0; i < nballs; ++i) fi >> cx[i] >> cy[i] >> r[i];
    if (!fi) { std::cerr << "Error reading " << argv[1] << "\n"; return 1; }
    const int anchor = withdepot ? ntot - 1 : 0;

    // ---- sequence, rotated so that the anchor (depot or disk 1) is first ----
    std::ifstream fs(argv[3]);
    std::vector<int> seq(ntot);
    int pos = -1;
    for (int k = 0; k < ntot; ++k) {
        if (!(fs >> seq[k]) || seq[k] < 1 || seq[k] > ntot) {
            std::cerr << "Invalid or short sequence: " << argv[3] << "\n";
            return 1;
        }
        --seq[k];
        if (seq[k] == anchor) pos = k;
    }
    if (pos < 0) { std::cerr << "Anchor node missing from the sequence\n"; return 1; }
    std::vector<int> tour(ntot);
    for (int k = 0; k < ntot; ++k) tour[k] = seq[(pos + k) % ntot];

    // ---- initial points ----
    std::ifstream fp(argv[4]);
    std::vector<double> px(ntot), py(ntot);
    for (int i = 0; i < ntot; ++i) fp >> px[i] >> py[i];
    if (!fp) { std::cerr << "Error reading " << argv[4] << "\n"; return 1; }

    // ---- data in visiting order ----
    std::vector<double> tcx(ntot), tcy(ntot), tr(ntot), x(ntot), y(ntot);
    for (int k = 0; k < ntot; ++k) {
        const int i = tour[k];
        tcx[k] = cx[i]; tcy[k] = cy[i]; tr[k] = r[i];
        x[k] = px[i];   y[k] = py[i];
    }

    bcd_cetsp::Options opt;
    opt.max_sweeps = 200;
    const bcd_cetsp::Result res =
        bcd_cetsp::bcd(ntot, tcx.data(), tcy.data(), tr.data(), x.data(), y.data(), opt);
    const double len = bcd_cetsp::tour_length(ntot, x.data(), y.data());

    // ---- output, in node order ----
    for (int k = 0; k < ntot; ++k) { px[tour[k]] = x[k]; py[tour[k]] = y[k]; }
    FILE* fo = std::fopen(argv[5], "w");
    if (!fo) { std::cerr << "Cannot write " << argv[5] << "\n"; return 1; }
    for (int i = 0; i < ntot; ++i) std::fprintf(fo, " %24.16E %24.16E\n", px[i], py[i]);
    std::fclose(fo);

    std::printf("%20.6f %8d\n", len, res.sweeps);
    return 0;
}
