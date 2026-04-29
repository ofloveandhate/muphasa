//
//  main.cpp
//  mph
//
//  Created by Oliver on 2020-04-09.
//  Copyright © 2020 Oliver. All rights reserved.
//

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <iostream>
#include <time.h>


#include "main.h"
#include "presentation.h"



double res_memory=0;
double virt_memory=0;



/*
 Testing
 */

std::vector<double> run_gb_singatures_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(image, input_matrix, debug);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << minimal_presentation.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = minimal_presentation.first.size();
    return std::vector<double>{time, res_memory, size};
}

std::vector<double> run_gb_shreyer_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    Matrix kernel = computeGroebnerBases_gradeopt(input_matrix).second;
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation = compute_presentation_schreyer(image, kernel);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << presentation.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = presentation.first.size();
    return std::vector<double>{time, res_memory, size};
}


/* Input */

void print_usage_and_exit(int exit_code) {
    std::cerr
    << "Usage: "
    << "mph "
    << "[options] [filename]" << std::endl
    << std::endl
    << "Options:" << std::endl
    << std::endl
    << "  --help           print this screen" << std::endl
    << "  --dim <k>        compute presentation matrix of the k-th persistent homology module" << std::endl
    << "  --firep <k>      input file with rivet <firep> file format" << std::endl
    << std::endl;
    exit(exit_code);
}



int main(int argc, char** argv) {
    const char* filename = nullptr;

    int dim_max = 1;
    bool firep = false;

    if(argc==1){
        print_usage_and_exit(0);
    }

    for (index_t i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage_and_exit(0);
        } else if (arg == "--dim") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            dim_max = (int)std::stol(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--firep"){
            firep = true;
        } else {
            if (filename) { print_usage_and_exit(-1); }
            filename = argv[i];
        }
    }

    std::pair<Matrix, std::vector<grade_t>> presentation;
    if(firep){
        Matrix high_matrix, low_matrix;
        std::ifstream file_stream(filename);
        if (file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }
        read_input_file<SyzColumn, SyzColumn>(file_stream, high_matrix, low_matrix);
        presentation = computeMinimalPresentation_3p(high_matrix, low_matrix);
    } else{
        std::ifstream file_stream(filename);
        if (filename && file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }

        // TODO: how should choice of metrics and filters be specified in input?

        std::vector<Metric*> metrics;
        metrics.push_back(new SquaredEuclideanMetric());
        std::vector<Filter*> filters;
        filters.push_back(new Filter());
        std::vector<input_t> max_metric_values;
        max_metric_values.push_back(10);
        std::vector<std::vector<input_t>> points = read_point_cloud(file_stream, 0);
        std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(points, metrics, filters, max_metric_values, dim_max);
        verify_kernel(boundary_matrices.first, boundary_matrices.second);
        for(auto& metric : metrics){
            delete metric;
        }
        for(auto& filter : filters){
            delete filter;
        }
        for(size_t i=0; i<boundary_matrices.second.size(); i++){
            boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
        }
        presentation = computeMinimalPresentation_3p(boundary_matrices.first, boundary_matrices.second);
    }


    exit(0);
}
