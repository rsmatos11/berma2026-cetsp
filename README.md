# CETSP benchmark instances and solutions

Benchmark instances of the Close Enough Traveling Salesman Problem (CETSP),
with the best solutions obtained in:

> E. G. Birgin, P. C. R. Ertel, R. S. Matos.
> *Using block coordinate descent for solving the close enough traveling
> salesman problem.* (journal and DOI to be added)

Web page with the tables of results: https://rsmatos11.github.io/berma2026-cetsp/

## Contents

| Folder | Content |
|--------|---------|
| `Instances/Behdani2014/` | 720 instances of Behdani and Smith (2014), with radii 0.25, 0.50 and 1.00 |
| `Instances/Mennell2009/` | 62 instances of Mennell (2009): 48 with a common radius and 14 with variable radii |
| `Solutions/Files/` | best solution found by MA-CETSP-BCD for each instance; `Mennell2009_ILS-BCD/`: solutions of ILS-BCD for the instances of Mennell (2009) |
| `Solutions/Figures/` | figure of each solution (PDF; MetaPost sources in `mpfiles/`), organized in the same way |
| `Codes/` | source code (see `Codes/README.md`) |

## Instance format

```
n            number of disks
x0 y0        depot
x  y  r      one line per disk: center and radius
...
```

Disk k is the k-th disk line, after the depot line.

The instance files were reformatted (aligned columns); the values are those
of the original instances.

## Solution format

```
instance = kroD100
overlap ratio = 0.02
method = MA-CETSP-BCD
tour length = 159.037358003615
number of disks = 99

 disk                x                y
depot   2.780000000000   1.650000000000
    2   2.794857497633   2.242342151742
  ...
```

The lines after the header give the visit points in visiting order; the
tour starts and ends at the depot. Disks are numbered as in the instance
file.

## References

- B. Behdani, J. C. Smith. An integer-programming-based approach to the
  close-enough traveling salesman problem. INFORMS Journal on Computing
  26(3):415–432, 2014.
- W. P. Coutinho, R. Q. do Nascimento, A. A. Pessoa, A. Subramanian.
  A branch-and-bound algorithm for the close-enough traveling salesman
  problem. INFORMS Journal on Computing 28(4):752–765, 2016.
- W. K. Mennell. Heuristics for solving three routing problems:
  close-enough traveling salesman problem, close-enough vehicle routing
  problem, and sequence-dependent team orienteering problem. PhD thesis,
  University of Maryland, 2009.
- Z. Lei, J.-K. Hao. An effective memetic algorithm for the close-enough
  traveling salesman problem. Applied Soft Computing 153:111266, 2024.
