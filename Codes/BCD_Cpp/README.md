# BCD for the CETSP (C++)

Block Coordinate Descent (BCD) for computing the optimal visit points of a
fixed visiting sequence in the Close Enough Traveling Salesman Problem
(CETSP). Given disks B(c_k, r_k), k = 0..n-1, in visiting order, BCD
minimizes the length of the closed tour x_0 -> x_1 -> ... -> x_{n-1} -> x_0
with x_k in B(c_k, r_k), optimizing one visit point at a time in cyclic
sweeps. Each single-disk subproblem is solved in closed form, when the
segment joining the two neighbors meets the disk, or by the adaptive
multi-start 1D Newton method on the boundary angle.

This is the C++ implementation described in:

> E. G. Birgin, P. C. R. Ertel, R. S. Matos.
> *Using block coordinate descent for solving the close enough traveling
> salesman problem.* (journal and DOI to be added)

It is the implementation used in the MA-CETSP-BCD experiments of the paper.
It implements the same method as the Fortran code used in the ILS
experiments (see "Relation to the Fortran code" below).

## Files

| File | Content |
|------|---------|
| `bcd_cetsp.hpp`, `bcd_cetsp.cpp` | the library |
| `bcd_rota.cpp` | example program: reads an instance, a visiting sequence and initial points, runs BCD and writes the visit points |
| `example/` | an example: instance kroD100 with overlap ratio 0.02, a visiting sequence and initial points |
| `Makefile` | builds the example program |

## Build

```bash
make            # builds ./bcd_rota (g++ -O3 -std=c++17)
```

To use the library in your own code, compile `bcd_cetsp.cpp` together with
it and include `bcd_cetsp.hpp`.

## Usage

```cpp
#include "bcd_cetsp.hpp"

// disks in visiting order: centers (cx[k], cy[k]), radii r[k];
// (x[k], y[k]) must start inside disk k (the centers are a valid choice)
bcd_cetsp::Options opt;                 // defaults: see below
bcd_cetsp::Result res = bcd_cetsp::bcd(n, cx, cy, r, x, y, opt);
// res.length: tour length; res.sweeps: number of sweeps;
// (x[k], y[k]) now hold the optimized visit points
```

Other functions: `mindisk` (single-disk subproblem: the point z of a disk
minimizing |z - p| + |z - q|), `tour_length`, and the plain C interface
`bcd_cetsp_c` (for use from C, Fortran, Python, ...). Coordinates are
normalized internally, so any unit may be used. A depot is a disk with
radius 0.

## Options

| Field | Default | Meaning |
|-------|---------|---------|
| `max_sweeps` | 200 | maximum number of sweeps |
| `rel_tol` | sqrt(machine epsilon) | stop when one sweep reduces the length by a relative amount below `rel_tol` |
| `window` | 0 (off) | if > 0, replaces the rule above: stop when the last `window` sweeps together reduce the length by a relative amount below `window_tol` |
| `window_tol` | 1e-6 | see `window` |
| `tie_rule` | `Maximin` | point chosen when the segment joining the two neighbors meets the disk (`Maximin`, `Midpoint` or `Projection`) |

The defaults reproduce the stopping rule of the Fortran code. In the
MA-CETSP-BCD experiments of the paper, each visiting sequence was evaluated
with `rel_tol = 1e-6` and `max_sweeps = 200`, and the best solution was
polished at the end with `window = 1000` and `window_tol = 1e-6`. In
instances with highly overlapping disks, BCD may reduce the tour length
very slowly but steadily, and the window rule avoids stopping too early.

## Example

```bash
./bcd_rota example/kroD100.dat 1 example/sequence.txt example/centers.txt points.txt
```

Expected output (tour length and number of sweeps):

```
          159.037364       32
```

The visit points are written to `points.txt`. See the header of
`bcd_rota.cpp` for the file formats. The visiting sequence is the best tour
found by the ILS of the paper for this instance, and the initial points are
the disk centers.

## Relation to the Fortran code

The method is the same as in the Fortran code of the ILS. The differences
are only in the implementation of the 1D Newton method: `sqrt` instead of
`hypot`, first and second derivatives computed together, and a step
tolerance of 1e-8 instead of 1e-13 (since convergence is quadratic, the
error after the last step is about 1e-16). On six benchmark instances
(kroD100, rat195, lin318, pcb442 and dsj1000 with overlap ratio 0.02, and
bubbles9), starting from the disk centers, both codes give the same tour
lengths (to the printed precision) and the same number of sweeps; the
visit points differ by at most 4e-7, which occurs where the tour length is
insensitive to their position.
