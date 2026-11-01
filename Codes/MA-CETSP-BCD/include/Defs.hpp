/**
 * Defs.hpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#ifndef CETSP_DEFS_HPP
#define CETSP_DEFS_HPP

#include <vector>
#include <string>

typedef std::vector<std::vector<double>> Centers;
typedef std::vector<std::vector<int>> Neighbors;

const int ITERATION = 5000;             // max iteration
const double MAX_TIME = 36000;          // max running time
const int POPULATION_SIZE = 20;         // population size
const int DISTANCE_THRESHOLD = 5;       // min distance in population
const double FIT_BETA = 0.96;           // fitness function distance coef
const int NEIGHBOR_SIZE = 50;           // neighbors size of a target
const double EPSILON = 1e-4;            // approximation in geometry and local search
const double DELTA = 1e-3;              // approximation in accepting best solution
const double PI = 3.14159265358979323846;

const std::string INITIALIZATION = "KMEANS";    // RANDOM, KMEANS
const std::string SELECTION = "RANDOM";         // RANDOM, ROULETTE
const std::string CROSSOVER = "KSX";            // KSX, GAX, EAX
const std::string IMPROVEMENT = "BEST";         // FIRST, BEST
const std::string GREEDY_ALGO = "SQUEEZE";      // SQUEEZE, SPARSE
const std::string DISTANCE = "EDIT";

const bool LOG = false;                          // log more details
const int SEED = 0;                              // random seed

// paths, relative to the folder where the program is run
const std::string RES_DIR = "results/";            // directory of results
const std::string LKH_EXE = "LKH-2.0.10/LKH";      // path for LKH solver
const std::string LKH_TMP_ROOT = "tmp_lkh/lkh";    // LKH tmp output directory

#endif //CETSP_DEFS_HPP
