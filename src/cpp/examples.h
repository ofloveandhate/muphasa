//
//  examples.h
//  mph
//
//  Created by Oliver on 2020-10-26.
//  Copyright © 2020 Oliver. All rights reserved.
//

#ifndef examples_h
#define examples_h
#include <iostream>
#include <chrono>
#include <random>
#include <stdexcept>
#include "utils.h"
#include "grade.h"
#include "signatureColumn.h"
#include "IO.h"

grade_t get_unique_grade(int p, std::vector<hash_set<index_t>>& visited_grades){
    grade_t grade;
    for(size_t i=0; i<p; i++){
        while(true){
            int val = rand()%10000;
            if(visited_grades[i].find(val) == visited_grades[i].end()){
                grade.push_back(val);
                visited_grades[i].insert(val);
                break;
            }
        }
    }
    return grade;
}

grade_t get_random_grade(int p, std::uniform_int_distribution<int>& dist, std::mt19937& seq){
    grade_t grade;
    for(size_t i=0; i<p; i++){
        int val = dist(seq);
        grade.push_back(val);
    }
    return grade;
}

void get_random_column(int n, SignatureColumn& column, double density){
    grade_t grade;
    double lower_bound = 0;
    double upper_bound = 1;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;
    for(size_t i=0; i<n; i++){
        if(unif(re)< density){
            column.push(column_entry_t(1, i));
        }
    }
    return;
}

std::pair<Matrix, std::vector<grade_t>> random_boundary_matrix(int n, int m, int dim, int n_parameters, uint32_t seed=-1){
    /*
     Returns a random subset of the n-simplex boundary matrix at dimension 'dim' with random grades on the (dim-1)-simplices.
     */
    Matrix boundary;
    
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> uni(0,m-1);
    
    std::vector<grade_t> row_grades;
    for(size_t i=0; i<m; i++){
        row_grades.push_back(get_random_grade(n_parameters, uni, gen));
    }
    std::vector<grade_t> point_grades;
    for(size_t i=0; i<m; i++){
        grade_t column_grade = get_random_grade(n_parameters, uni, gen);
        point_grades.push_back(row_grades[i].join(column_grade));
    }
    for(size_t i=0; i<n; i++){
        hash_set<int> visited;
        int index = uni(gen);
        visited.insert(index);
        SignatureColumn column(point_grades[index], i);
        column.push(column_entry_t(1, index));
        grade_t row_grade = row_grades[index];
        for(size_t j=0; j<dim-1; j++){
            index = uni(gen);
            while(visited.find(index) != visited.end()){
                index = uni(gen);
            }
            visited.insert(index);
            column.push(column_entry_t(1, index));
            row_grade = row_grade.join(row_grades[index]);
            //column.grade = column.grade.join(point_grades[index]);
        }
        column.grade = row_grade;
        boundary.push_back(column);
    }
    return std::pair<Matrix, std::vector<grade_t>>(boundary, row_grades);
}


std::vector<input_t> get_random_next_step(std::vector<input_t> point, std::normal_distribution<double>& dist, std::mt19937& seq){
    /* Assume last coordinate is time */
    std::vector<input_t> point_;
    for(size_t i=0; i<point.size()-1; i++){
        point_.push_back(point[i] + dist(seq));
    }
    point_.push_back(point[point.size()-1] + 1);
    return point_;
}

std::vector<std::vector<input_t>> trajectory(size_t n_samples, size_t dim, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                                             (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                                         std::chrono::system_clock::now().time_since_epoch()
                                                                                              ).count() +
                                             (std::mt19937::result_type)
                                             std::chrono::duration_cast<std::chrono::microseconds>(
                                                                                                   std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                                                   ).count() );
    }
    
    std::mt19937 gen(seed);
    std::normal_distribution<double> d{0,1};
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    
    std::vector<std::vector<input_t>> points;
    std::vector<input_t> point;
    for(size_t i=0; i<dim+1; i++){
        point.push_back(uniform(gen));
    }
    for(size_t i=0; i<n_samples; i++){
        point = get_random_next_step(point, d, gen);
        points.push_back(point);
    }
    return points;
}

std::vector<std::vector<input_t>> trajectory_box(size_t n_samples, size_t dim, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::normal_distribution<double> d{0,0.2};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    std::vector<std::vector<input_t>> points;
    std::vector<input_t> point;
    for(size_t i=0; i<dim; i++){
        point.push_back(uniform(gen));
    }
    for(size_t i=0; i<n_samples; i++){
        point = get_random_next_step(point, d, gen);
        for(size_t j=0; j<point.size(); j++){
            if(point[j] < 0){
                point[j] = -point[j];
            }else if(point[j] > 1){
                point[j] = 2-point[j];
            }
        }
        points.push_back(point);
    }
    return points;
}


std::vector<std::vector<input_t>> time_varying_point_cloud(size_t n_points, size_t n_timesteps, size_t dim, uint32_t seed=-1){
    std::vector<std::vector<input_t>> points;
    for(size_t i=0; i<n_points; i++){
        std::vector<std::vector<input_t>> t;
        if(seed==-1){
            t = trajectory_box(n_timesteps, dim);
        }else{
            t = trajectory_box(n_timesteps, dim, seed+(uint32_t)i);
        }
        for(auto& p:t){
            points.push_back(p);
        }
    }
    return points;
}

std::vector<std::vector<std::vector<input_t>>> time_varying_point_cloud_trajectories(size_t n_points, size_t n_timesteps, size_t dim, uint32_t seed=-1){
    std::vector<std::vector<std::vector<input_t>>> trajectories;
    for(size_t i=0; i<n_points; i++){
        std::vector<std::vector<input_t>> t;
        if(seed==-1){
            t = trajectory_box(n_timesteps, dim);
        }else{
            t = trajectory_box(n_timesteps, dim, seed+(uint32_t)i);
        }
        trajectories.push_back(t);
    }
    return trajectories;
}

std::pair<Matrix, Matrix> get_random_boundary_matrix(size_t n_points, size_t n_dim, size_t n_parameters, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> uni(0,100);
    
    grade_t max_grade;
    for(size_t i=0; i<n_parameters; i++){
        max_grade.push_back(100);
    }
    
    
    std::vector<std::vector<input_t>> points;
    for(size_t i=0; i<n_points; i++){
        points.push_back(std::vector<input_t>{1});
    }
    std::vector<grade_t> grades;
    for(size_t i=0; i<n_points; i++){
        grades.push_back(get_random_grade(n_parameters, uni, gen));
    }
    return compute_boundary_matrices_grades(points, grades, max_grade, n_dim);
}




void critical_points_geometric_pres(){
    get_mem_usage(virt_memory, res_memory);
    double start_mem = res_memory;
    double metric_thresh = 3.5;
    int n = 10;
    std::vector<std::vector<double>> times;
    std::vector<std::vector<double>> resident_memory;
    std::vector<std::vector<double>> output_size;
    std::vector<double> boundary_matrix_size;
    while(n < 201){
        std::cout << "\n Iteration: " << n << std::endl;
        uint32_t seed = 1;
        std::vector<std::vector<input_t>> points = time_varying_point_cloud(10, 10*n, 3, seed);//time_varying_point_cloud(60, 100, 3, 1);
        std::vector<std::vector<input_t>> xpoints = points;
        for(auto& p:xpoints){
            p.pop_back();
        }
        
        std::vector<Metric*> metrics;
        metrics.push_back(new SquaredEuclideanMetric());
        std::vector<Filter*> filters;
        filters.push_back(new XFilter(3, 1));
        filters.push_back(new XFilter(3, -1));
        //filters.push_back(new XFilter(0, 1));
        //filters.push_back(new XFilter(0, -1));
        //filters.push_back(new XFilter(1, 1));
        //filters.push_back(new XFilter(1, -1));
        //filters.push_back(new XFilter(2, 1));
        //filters.push_back(new XFilter(2, -1));
        //filters.push_back(new DensityFilter(xpoints));
        
        std::vector<input_t> max_metric_values;
        max_metric_values.push_back(metric_thresh);
        
        std::pair<Matrix, Matrix> boundaries = compute_boundary_matrices(points, metrics, filters, max_metric_values, 1);
        
        for(size_t i=0; i<boundaries.second.size(); i++){
            boundaries.second[i].signature_index = i;
        }
        
        sort(boundaries.second.begin(), boundaries.second.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
             {
                 return lhs.grade.lt_colex(rhs.grade);
             });
        
        //Reindex image columns
        hash_map<size_t, size_t> map;
        for(size_t i=0; i<boundaries.second.size();i++){
            map[boundaries.second[i].signature_index] = i;
        }
        
        for(size_t i=0; i<boundaries.first.size(); i++){
            SyzColumn c;
            while(boundaries.first[i].get_pivot().get_index() != -1){
                c.push(column_entry_t(1, map[boundaries.first[i].get_pivot().get_index()]));
                boundaries.first[i].pop_pivot();
            }
            while(c.get_pivot().get_index() != -1){
                boundaries.first[i].push(c.get_pivot());
                c.pop_pivot();
            }
        }
        
        sort(boundaries.first.begin(), boundaries.first.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
             {
                 return lhs.grade < rhs.grade;
             });
        
        Matrix boundary = boundaries.second;
        Matrix image = compute_minimal_generating_set(boundaries.first);
        
        std::cout << "Size of boundary matrix: " << boundary.size() << " Size of image matrix: " << image.size() << std::endl;
        
        for(size_t i=0; i<boundary.size(); i++){
            boundary[i].signature_index = i;
            boundary[i].syzygy = SyzColumn();
            boundary[i].syzygy.push(column_entry_t(1, i));
        }
        
        //Test if image maps into kernel
        for(auto& column : image){
            SignatureColumn c(column);
            SignatureColumn z(column.grade, -1);
            while(c.get_pivot_index() != -1){
                z.plus(boundary[c.get_pivot_index()]);
                c.pop_pivot();
            }
            if(z.get_pivot_index() != -1){
                throw std::runtime_error("Image does not map to kernel");
            }
        }
        
        Matrix input_copy;
        for(auto& c : boundary){
            input_copy.push_back(SignatureColumn(c));
        }
        Matrix image_copy;
        for(auto& c : image){
            image_copy.push_back(SignatureColumn(c));
        }
        std::vector<double> gb_res = run_gb_singatures_pres(input_copy, image_copy);
        input_copy.clear();
        for(auto& c : boundary){
            input_copy.push_back(SignatureColumn(c));
        }
        image_copy.clear();
        for(auto& c : image){
            image_copy.push_back(SignatureColumn(c));
        }
        std::vector<double> gb_schreyer = run_gb_shreyer_pres(input_copy, image_copy);
        
        times.push_back(std::vector<double>{gb_schreyer[0], gb_res[0]});
        resident_memory.push_back(std::vector<double>{gb_schreyer[1]-start_mem, gb_res[1]-start_mem});
        output_size.push_back(std::vector<double>{gb_schreyer[2], gb_res[2]});
        boundary_matrix_size.push_back(boundary.size());
        
        std::cout << "$"<< metric_thresh << "$ & " << times[times.size()-1][0]/1000 << " & " << resident_memory[times.size()-1][0] << " & " << times[times.size()-1][1]/1000 << " & " << resident_memory[times.size()-1][1]  << "\n";
        
        n+= 10;
    }
    for(size_t i=0; i<times.size(); i++){
        std::cout << "$"<< metric_thresh+i*0.05 << "$ & " << times[i][0]/1000 << " & " << resident_memory[i][0] << " & " << times[i][1]/1000 << " & " << resident_memory[i][1] << "\n";
    }
    std::cout << "\n";
    std::cout << "\n Times GBS+Schreyer \n";
    for(auto& v : times){
        std::cout << v[0] << ", ";
    }
    std::cout << "\n Times PresentationPair \n";
    for(auto& v : times){
        std::cout << v[1] << ", ";
    }
    std::cout << "\n Resident memory GBS+Schreyer\n";
    for(auto& v : resident_memory){
        std::cout << v[0] << ", ";
    }
    std::cout << "\n Resident memory PresentationPair\n";
    for(auto& v : resident_memory){
        std::cout << v[1] << ", ";
    }
    
    std::cout << "\n Input size \n";
    for(auto& v : boundary_matrix_size){
        std::cout << v << ", ";
    }
}



#endif /* examples_h */
