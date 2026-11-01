/**
 * Algo.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "Algo.hpp"
#include <algorithm>
#include <iomanip>
#include "LocalSearch/Solver.hpp"

Algo::Algo(Parameters *params): population(params), data(params) {
    this->params = params;
    this->random = new Random(params->seed);
    this->iteration = params->iteration;
    this->patience_threshold = params->patience;
    this->timestamp = params->timestamp;
}

Algo::~Algo(){
    delete random;
}

void Algo::run() {
    // read data;
    Centers centers = data.getData();
    // set context
    population.setContext(centers, random, timestamp);

    auto start_run = std::chrono::high_resolution_clock::now();
    std::cout << std::setprecision(10);   // values in the log with 10 digits

    int iter = 0;
    int best_iter = 0;
    int patience = 0;
    bool improved = false;
    std::chrono::duration<double> best_running_time{0.0};
    std::string stop = "iterations";      // stopping rule that ended the run
    // init population
    List *best_solution = population.initPopulation();
    double best_solution_value = best_solution->getValue();

    // iteration
    while (iter++ < iteration) {
        std::cout << std::endl << "Iteration " << iter << " : " << std::endl;

        auto start_iter = std::chrono::high_resolution_clock::now();
        // next population
        best_solution = population.nextPopulation(patience);
        improved = best_solution->getValue() - best_solution_value < -DELTA ? true : false;
        best_solution_value = best_solution->getValue();

        auto end_iter = std::chrono::high_resolution_clock::now();

        if (improved) {
            best_iter = iter;
            patience = 0;
            best_running_time = end_iter - start_run;
            if (LOG) data.write(best_solution, iter, std::to_string(best_running_time.count()));
        } else {
            ++patience;
        }

        std::cout << "[LOG] iter: " << iter
                  << " best_iter: " << best_iter
                  << " best_value: " << best_solution->getValue()
                  << " iter_time: " << std::chrono::duration<double> (end_iter - start_iter).count()
                  << " total_time: " << std::chrono::duration<double> (end_iter - start_run).count()
                  << std::endl;

        if (patience >= patience_threshold) {
            std::cout << std::endl << "[STOP] best solution hasn't been improved since " << patience_threshold << " iterations" << std::endl;
            stop = "patience";
            break;
        }

        if (std::chrono::duration<double> (end_iter - start_run).count() > params->max_time) {
            std::cout << std::endl << "[STOP] running time exceeds " << params->max_time << " seconds" << std::endl;
            stop = "time";
            break;
        }
    }

    // final polishing of the best solution, after the search
    const double value_before_polish = best_solution->getValue();
    auto start_polish = std::chrono::high_resolution_clock::now();
    Solver polisher;
    polisher.setContext(centers);
    polisher.polish(best_solution);
    const double polish_time = std::chrono::duration<double> (std::chrono::high_resolution_clock::now() - start_polish).count();

    auto end_run = std::chrono::high_resolution_clock::now();

    if (!LOG)  data.write(best_solution, best_iter, std::to_string(best_running_time.count()));

    std::cout << std::endl
            << "[SUMMARY] instance: " << data.getFilename()
            << " best_value: " << best_solution->getValue()
            << " best_time: " << std::to_string(best_running_time.count())
            << " total_time: " << std::chrono::duration<double> (end_run - start_run).count()
            << " result_file: " << data.getResultFilename()
            << std::endl;

    // one line with the results of the run
    // (key: value, values with full precision)
    std::cout << std::setprecision(17)
              << "[RESULT] instance: " << data.getFilename()
              << " seed: " << params->seed
              << " nodes: " << centers.size()
              << " best_value: " << best_solution->getValue()
              << " value_before_polish: " << value_before_polish
              << " best_iter: " << best_iter
              << " iterations: " << std::min(iter, iteration)
              << " stop: " << stop
              << " best_time: " << best_running_time.count()
              << " total_time: " << std::chrono::duration<double> (end_run - start_run).count()
              << " polish_time: " << polish_time
              << " result_file: " << data.getResultFilename()
              << std::endl;
    std::cout << "[TIMES] greedy: " << stage_times.greedy
              << " lkh: " << stage_times.lkh
              << " ls: " << stage_times.ls
              << " socp: " << stage_times.socp
              << " vnd_calls: " << stage_times.vnd_calls
              << std::endl;
}
