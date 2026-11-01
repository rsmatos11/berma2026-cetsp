/**
 * Data.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "Utils/Data.hpp"
#include <iomanip>
#include <cstdlib>

Data::Data(Parameters *params) {
    this->instance_file = params->instance_file;
    this->filename = std::filesystem::path(instance_file).stem().string();
    // the seed in the name keeps runs started in the same second apart
    this->timestamp = params->timestamp + "s" + std::to_string(params->seed);
    this->result_filename = "";
    res_dir = RES_DIR;
}

void Data::read() {
    std::cout << "instance : " << filename << std::endl;
    // instance format: n; depot x y; n lines x y r
    std::ifstream in(instance_file);
    int n = 0;
    std::string sx, sy;
    double x, y, r;
    in >> n >> sx >> sy;
    if (!in) {
        std::cout << "ERROR : read data " << instance_file << std::endl;
        std::exit(1);
    }
    // DEPOT_FLOAT=1: depot in single precision, as in the reader of the original
    // MA-CETSP for the instances of Mennell (needed to reproduce the paper exactly)
    if (std::getenv("DEPOT_FLOAT") && std::atoi(std::getenv("DEPOT_FLOAT")) != 0) {
        x = std::stof(sx);
        y = std::stof(sy);
    } else {
        x = std::stod(sx);
        y = std::stod(sy);
    }
    centers.push_back(std::vector<double> {x, y, 0});
    for (int i = 0; i < n && in >> x >> y >> r; ++i)
        centers.push_back(std::vector<double> {x, y, r});
    if (!in || (int)centers.size() != n + 1) {
        std::cout << "ERROR : read data " << instance_file << std::endl;
        std::exit(1);
    }
    in.close();
    reduceSize();
}

void Data::write(List* solution, int iter, std::string running_time) {
    try {
        if (! std::filesystem::exists(res_dir)) {
            std::filesystem::create_directory(res_dir);
        }
        std::string running_dir = res_dir + filename + "-" + this->timestamp + "/";
        if (! std::filesystem::exists(running_dir)) {
            std::filesystem::create_directory(running_dir);
        }
        std::ofstream out;
        result_filename = filename + "-" + this->timestamp + "-"+ std::to_string(iter);
        if (LOG){
            std::cout << "file name : " << result_filename << std::endl;
        }
        out.open(running_dir + result_filename +  ".txt");
        out << std::setprecision(17);  // full precision
        Node* p = solution->head();
        for (int i = 0; i < solution->size(); ++i) {
            out << p->id << ",";
            p = p->next;
        }
        out << "\n";
        out << "value : " << solution->getValue() << " " << "running time : " << running_time <<  "\n";
        p = solution->head();
        for (int i = 0; i < solution->size(); ++i) {
            out << p->id << " " << p->x << " " << p->y << "\n";
            p = p->next;
        }
        out.close();
    } catch (std::exception e) {
        std::cout << "ERROR : write result" << std::endl;
    }
}

void Data::reduceSize() {
    // reduce the vertices which are covered by depot
    double x0 = centers[0][0], y0 = centers[0][1];
    auto it = centers.begin() + 1;
    while (it != centers.end()) {
        double x = it->at(0), y = it->at(1), r = it->at(2);
        if (pow(x-x0, 2)  + pow(y-y0, 2) < pow(r, 2)) {
            it = centers.erase(it);
        } else {
            ++it;
        }
    }
    // reduce the centers which are covered by others
    auto iti = centers.begin() + 1;
    while (iti != centers.end()) {
        auto itj = iti + 1;
        while (itj != centers.end()) {
            double x1 = iti->at(0), y1 = iti->at(1), r1 = iti->at(2);
            double x2 = itj->at(0), y2 = itj->at(1), r2 = itj->at(2);
            double dist = sqrt(pow(x1-x2, 2)  + pow(y1-y2, 2));
            if (dist <= r2 - r1) {
                itj = centers.erase(itj);
            } else if (dist <= r1 - r2) {
                iti = centers.erase(iti);
                --iti;
                break;
            } else {
                ++itj;
            }
        }
        ++iti;
    }
}

Centers Data::getData() {
    read();
    return centers;
}

std::string Data::getResultFilename() {
    return result_filename;
}

std::string Data::getFilename() {
    return filename;
}