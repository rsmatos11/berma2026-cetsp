# MA-CETSP-BCD

Memetic algorithm for the Close Enough Traveling Salesman Problem (CETSP) in
which the visit points are computed by Block Coordinate Descent (BCD). Code used
in:

> E. G. Birgin, P. C. R. Ertel, R. S. Matos.
> *Using block coordinate descent for solving the close enough traveling
> salesman problem.* (journal and DOI to be added)

It is a modified version of the code of
[MA-CETSP](https://github.com/leizy1008/MA-CETSP) (Z. Lei, J.-K. Hao, Applied Soft
Computing 153:111266, 2024). The original copyright notice and the MIT license
are in `LICENSE`.

Main changes: the SOCP is solved by BCD instead of Gurobi; the local search uses
the exact solution of the single-disk subproblem; two BCD sweeps replace the
greedy pass after LKH; and the best solution is polished by BCD at the end of
the search. See Section 5.2 of the paper.

## Requirements

- `g++` with C++17 support and `make`
- [LKH](http://akira.ruc.dk/~keld/research/LKH/) version 2.0.10, by Keld Helsgaun,
  not included here: place the folder `LKH-2.0.10` in this folder. Another version
  may give results different from those of the paper.

## Build and run

```bash
./build.sh
./MA-CETSP-BCD -f ../../Instances/Behdani2014/radius_0.50/CETSP-12-0.txt -s 1 -r 5000 -t 7200
```

Run the program from this folder. Options: `-f` instance file, `-s` random seed,
`-r` maximum number of generations (the run also stops after `r`/10 generations
without improvement), `-t` time limit in seconds. The best solution is written
to `results/`, with the nodes numbered after the removal of the disks that
contain the depot or another disk (0 is the depot).

## Reproducing the paper

The paper used seeds 1 to 20 with the options above. For the instances of
Mennell, add `DEPOT_FLOAT=1`, which reads the depot in single precision, as the
original MA-CETSP does for these instances:

```bash
DEPOT_FLOAT=1 ./MA-CETSP-BCD -f ../../Instances/Mennell2009/equal_radius/Overlap_ratio_002/kroD100.txt -s 1 -r 5000 -t 7200
```

This run gives tour length 159.03735801628068.

The variants of the ablation study use the same program, with these additional
environment variables:

| Variant | Environment variables |
|---------|-----------------------|
| V1 | `LS_MINDISK=0 LS_PRUNE=0 GREEDY_BCD_SWEEPS=0 GREEDY_BEFORE=squeeze` |
| V2 | `GREEDY_BCD_SWEEPS=0 GREEDY_BEFORE=squeeze` |
| V3 | `GREEDY_BEFORE=squeeze` |
