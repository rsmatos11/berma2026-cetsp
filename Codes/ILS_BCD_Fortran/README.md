# ILS with BCD for the Close Enough Traveling Salesman Problem

Iterated Local Search (ILS) for the Close Enough Traveling Salesman Problem
(CETSP), in which the visit points of each candidate visiting sequence are
computed by a Block Coordinate Descent (BCD) method with a multi-start 1D
Newton solver for the single-disk subproblem.

This is the code used in:

> E. G. Birgin, P. C. R. Ertel, R. S. Matos.
> *Using block coordinate descent for solving the close enough traveling
> salesman problem.* (journal and DOI to be added)

## Requirements

- `gfortran` and `gcc`
- Optional: IBM ILOG CPLEX Studio, only to compute the visit points with
  CPLEX instead of BCD (`usecplex=1`). Set `CPLEX_STUDIO_DIR` to the CPLEX
  installation folder; this also requires `g++`.

## Quick start

```bash
./run-cetsp.sh
```

The script compiles the program and runs it on every instance in
`inst-cetsp/` (two examples are included). Edit the variables at the top of
`run-cetsp.sh` to change the instance folder, the time limit, the seed or
the parameters.

## Instance format

Plain text (`.dat`):

```
n            number of disks
x0 y0        depot (only if withdepot=1)
x  y  r      one line per disk: center and radius
...
```

With `withdepot=1` the depot is a radius-zero disk fixed at the first
position of the tour; with `withdepot=0`, disk 1 is used instead.

## Parameters

Set in `run-cetsp.sh` (values used in the paper):

| Variable      | Symbol      | Value  | Meaning |
|---------------|-------------|--------|---------|
| `tlim`        |             | 36000  | time limit per instance (s) |
| `seed`        |             | 123456 | random seed |
| `ncool_steps` | κ           | 20     | number of threshold reductions |
| `eps`         | η₀          | 0.02   | initial acceptance threshold |
| `h`           | h           | 0.85   | threshold reduction factor (η ← hη) |
| `strlen`      | L           | 3      | maximum string length (string relocation) |
| `reverse`     |             | 0      | 1 = relocated strings may be reversed |
| `cardZ`       | \|Z\|       | 10     | positions reassigned by the Gutin neighborhood |
| `numdbmoves`  | ℓ           | 5      | double-bridge moves per perturbation |
| `bsk`         | k           | 3      | Balas–Simonetti parameter (2 or 3) |
| `vnd_type`    |             | 0      | 0 = fixed neighborhood order, 1 = random |
| `bestfirst`   |             | 0      | 0 = first improvement, 1 = best improvement |
| `ncand_max`   | K           | 16     | size of the candidate lists |
| `maxbcdcic`   | s_max       | 1      | maximum BCD sweeps per evaluation in VND |
| `refine`      |             | 1      | 1 = complete BCD (up to 200 sweeps) after each VND |
| `withdepot`   |             | 1      | 1 = instance files have a depot line |
| `usecplex`    |             | 0      | 0 = BCD, 1 = CPLEX for the visit points |

The run stops at the time limit or, for instances with at most 200 nodes,
after 5,000 consecutive iterations without improvement.

The κ threshold reductions are scheduled from the time limit: the first
iterations are used to estimate the time per iteration, and the reductions
are spread evenly over the iterations expected to fit in the remaining
time. Therefore, the results depend on the speed of the machine, and runs
with the same seed and time limit may differ between machines (and,
occasionally, between runs on the same machine).

## Output

For each instance:

- `solutions_cetsp/<instance>-bcd.txt` — log of the run (every improvement
  of the best solution and, if the run stops by stagnation, a final summary)
- `last_incumbent/<instance>-bcd.csv` — best solution found: tour length,
  CPU time until it was found, counters and visiting order
- `solution_report/<instance>.txt` — best tour: length, visiting order and
  visit points (original units, full precision)
- `figures/<instance>-bcd.mp` — MetaPost figure of the best solution
- `results.csv` — one line per instance (same columns as `last_incumbent/`)

## Source files

| File | Content |
|------|---------|
| `main.f90` | main program |
| `instance.f90` | instance, parameters, figure output |
| `blockcoordinatedescent.f90` | BCD and the single-disk solver (1D Newton) |
| `metaheuristic.f90` | nearest neighbor, VND neighborhoods, VND, ILS |
| `assct.f` | assignment problem (Gutin neighborhood) |
| `cetsp_fixedtour_interface.f90`, `.cpp` | interface to CPLEX (optional) |
| `cetsp_fixedtour_interface_nocplex.c` | replaces the CPLEX interface when CPLEX is not available |

## Third-party code

`assct.f` is Algorithm 548 of the ACM Collected Algorithms
(G. Carpaneto and P. Toth, *Solution of the assignment problem*,
ACM Trans. Math. Software 6(1):104–111, 1980), with the maximum
dimension changed to 1000 and the "infinity" constant reduced to fit a
32-bit integer. It is distributed under the ACM Software License, which
allows non-commercial use only.
