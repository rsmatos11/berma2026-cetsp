#!/bin/bash
set -euo pipefail

# =====================================================================
#  ILS with BCD for the Close Enough Traveling Salesman Problem (CETSP)
#
#  Compiles the program and runs it, with one seed, on every instance
#  in INSTDIR. For each instance it writes:
#    solutions_cetsp/<instance>-<engine>.txt   full log of the run
#    last_incumbent/<instance>-<engine>.csv    best solution found
#    results.csv                               one line per instance
#
#  Source files (compiled in this order, due to module dependencies):
#    cetsp_fixedtour_interface.f90/.cpp   interface to the CPLEX fixed-tour
#                                         solver (used only if usecplex=1)
#    cetsp_fixedtour_interface_nocplex.c  used instead of the .cpp file when
#                                         CPLEX is not available
#    assct.f                              assignment problem (Gutin neighborhood)
#    instance.f90                         instance and parameters
#    blockcoordinatedescent.f90           BCD with multi-start 1D Newton
#    metaheuristic.f90                    ILS and VND
#    main.f90                             main program
#
#  Instance format (.dat): first line, the number of disks; if withdepot=1,
#  second line, the depot "x y"; then one line "x y r" per disk.
#
#  Requirements: gfortran and gcc. CPLEX is optional: it is used only if
#  the environment variable CPLEX_STUDIO_DIR points to an IBM ILOG CPLEX
#  Studio installation (needed only for usecplex=1; also requires g++).
#
#  Usage:  ./run-cetsp.sh
# =====================================================================

EXECUTABLE=ils_cetsp
INSTDIR=inst-cetsp

# ---- instance and engine ----
withdepot=1      # 1 = instance files have a depot line; 0 = they do not
usecplex=0       # 0 = visit points computed by BCD; 1 = by CPLEX (SOCP)

# ---- run ----
tlim=36000       # time limit per instance (seconds)
seed=123456      # random seed
seedrun=1        # run index (used together with the seed)

# ---- ILS parameters (values used in the paper) ----
ncool_steps=20   # kappa: number of threshold reductions over the run
eps=0.02         # eta_0: initial acceptance threshold
h=0.85           # reduction factor of the threshold (eta <- h * eta)
strlen=3         # L: maximum string length (string relocation neighborhood)
reverse=0        # 1 = relocated strings may be reversed; 0 = no
cardZ=10         # |Z|: positions reassigned by the Gutin neighborhood
numdbmoves=5     # ell: double-bridge moves per perturbation
bsk=3            # k: Balas-Simonetti parameter (2 or 3)
vnd_type=0       # 0 = neighborhoods in fixed order; 1 = random order at each call
bestfirst=0      # 0 = first improvement; 1 = best improvement
ncand_max=16     # K: size of the candidate lists
maxbcdcic=1      # s_max: maximum BCD sweeps per evaluation in VND
refine=1         # 1 = complete BCD (up to 200 sweeps) after each VND; 0 = no

# ---- compilation ----
FLAGS="-O3 -Wall"
FSRC="cetsp_fixedtour_interface.f90 assct.f instance.f90 blockcoordinatedescent.f90 metaheuristic.f90 main.f90"

if [ -n "${CPLEX_STUDIO_DIR:-}" ]; then
    CPLEXINC="-I${CPLEX_STUDIO_DIR}/cplex/include -I${CPLEX_STUDIO_DIR}/concert/include"
    CPLEXLIB="-L${CPLEX_STUDIO_DIR}/cplex/lib/x86-64_linux/static_pic \
-L${CPLEX_STUDIO_DIR}/concert/lib/x86-64_linux/static_pic \
-lilocplex -lconcert -lcplex -lpthread -ldl -lm -lstdc++"
    g++ -O2 -c cetsp_fixedtour_interface.cpp -o cetsp_fixedtour_interface_cpp.o $CPLEXINC -DIL_STD
else
    if [ "$usecplex" = "1" ]; then
        echo "ERROR: usecplex=1 requires CPLEX (set CPLEX_STUDIO_DIR)." >&2
        exit 1
    fi
    echo "CPLEX_STUDIO_DIR not set: compiling without CPLEX (BCD engine only)."
    CPLEXLIB=""
    gcc -O2 -c cetsp_fixedtour_interface_nocplex.c -o cetsp_fixedtour_interface_nocplex.o
fi
gfortran $FLAGS -c $FSRC
gfortran $FLAGS *.o -o "$EXECUTABLE" $CPLEXLIB
rm -f ./*.o ./*.mod

# ---- runs ----
export ilstlim=$tlim
export ILSUSECPLEX=$usecplex

enginedir="bcd"
if [ "$usecplex" = "1" ]; then
    enginedir="cplex"
fi

mkdir -p solutions_cetsp last_incumbent figures solution_report

RESULTS_CSV="results.csv"
if [ ! -f "$RESULTS_CSV" ]; then
    echo "Instance,NDiscs,IncumbentN,TourCost,CPUTime,ILSiter,VNDcalls,NSearchs,N2opt,N3opt,NSwap,NInsertion,NStringReloc,NDoubleBridge,NGutin,NBalasSimonetti,HCcalls,BCDcalls,BCDsweeps,CPLEXcalls,CPLEXiter,Tour" > "$RESULTS_CSV"
fi

mapfile -t instfiles < <(printf '%s\n' "${INSTDIR}"/* | sort -V)

for inst in "${instfiles[@]}"
do
    instname=$(basename "$inst" .dat)
    outfile="solutions_cetsp/${instname}-${enginedir}.txt"

    echo "== Running: \"$inst\" (seed $seed, run $seedrun, tlim ${tlim}s) =="
    set +e
    {
        echo \"$inst\"
        echo "${withdepot}"
        echo "${seed} ${seedrun}"
        echo "${ncool_steps} ${eps} ${h} ${strlen} ${reverse} ${cardZ} ${numdbmoves} ${bsk} ${vnd_type} ${bestfirst} ${ncand_max} ${maxbcdcic} ${refine}"
    } | timeout ${tlim}s ./"$EXECUTABLE" > "$outfile" 2>&1
    run_status=$?
    set -e

    # 0 = stopped by stagnation; 124 = stopped by the time limit.
    if [ "$run_status" -ne 0 ] && [ "$run_status" -ne 124 ]; then
        echo "ERROR: $EXECUTABLE exited with status $run_status for $inst" >&2
        exit 1
    fi

    lastfile="last_incumbent/${instname}-${enginedir}.csv"
    if [ -f "$lastfile" ]; then
        tail -n +2 "$lastfile" >> "$RESULTS_CSV"
    else
        echo "WARNING: no last_incumbent record found for $inst" >&2
    fi
done
