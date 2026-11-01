#!/bin/bash
# Builds LKH (LKH-2.0.10 must be in this folder; see README.md) and MA-CETSP-BCD.
set -euo pipefail
cd "$(dirname "$0")"
[ -d LKH-2.0.10 ] || { echo "LKH-2.0.10 not found in this folder: see README.md"; exit 1; }
[ -x LKH-2.0.10/LKH ] || make -C LKH-2.0.10
INC="-Iinclude -Iinclude/Genetic -Iinclude/Genetic/Crossover -Iinclude/LocalSearch -Iinclude/Utils -Ibcd"
SRC="$(find src -name '*.cpp' | sort)"
g++ -O3 -flto -std=c++17 $INC $SRC bcd/bcd_cetsp.cpp -o MA-CETSP-BCD
echo "OK: ./MA-CETSP-BCD"
