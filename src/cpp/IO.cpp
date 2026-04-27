#include "IO.h"



void build_VR_subcomplex(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     points {std::vector<std::vector<input_t>>} -- a list of points in space.
     metrics {std::vector<Metric>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter>} -- a list of filter functions on the points that should be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    for (size_t j = vertices.back() + 1; j < points.size(); j++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        for (size_t k = 0; k < vertices.size(); k++) {
            size_t p = vertices[k];
            for(size_t metric_index = 0; metric_index<metrics.size(); metric_index++){
                input_t d = (*metrics[metric_index]).eval(points[p], points[j]);
                if(d < max_metric_values[metric_index]){
                    curr_values[metric_index] = d < curr_values[metric_index] ? curr_values[metric_index] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(!should_add){
                break;
            }
        }
        
        if(should_add){
            for(size_t filter_index=0; filter_index<filters.size(); filter_index++){
                curr_values[metrics.size()+filter_index] = (*filters[filter_index]).eval(points[j]) < curr_values[metrics.size()+filter_index] ? curr_values[metrics.size()+filter_index] : (*filters[filter_index]).eval(points[j]);
            }
            vertices.push_back(j);
            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}


void build_VR_subcomplex_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, Metric* metric, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, input_t& max_metric_value, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex for a dynamic metric space in the sense of Memoli/Kim (2020).
     
     trajectory_dms {std::vector<std::vector<std::vector<input_t>>>} -- a list of trajectories for each point
     metric {Metric} -- a metric on the points that should be used to contruct the complex
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed. KATE: not sure what this is used for
     max_metric_value {std::vector<input_t>} -- the maximum allowed distance for the metric
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    
    // Add the current simplex to the relevant list and, if it's big enough, return
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    for (size_t new_vertex = vertices.back() + 1; new_vertex < trajectory_dms[0].size(); new_vertex++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        
        for(size_t time_index = 0; time_index<trajectory_dms.size(); time_index++){
            should_add = true;
            for (size_t old_vertex_index = 0; old_vertex_index < vertices.size(); old_vertex_index++) {
                size_t old_vertex = vertices[old_vertex_index];
                input_t d = trajectory_dms[time_index][new_vertex][old_vertex];
                if(d < max_metric_value){
                    curr_values[0] = d < curr_values[0] ? curr_values[0] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(should_add){
                break;
            }
        }
        
        if(should_add){
            vertices.push_back(new_vertex);
            build_VR_subcomplex_spatiotemporal(trajectory_dms, metric, vertices, curr_values, max_metric_value, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}


void build_VR_subcomplex_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& point_grades, std::vector<size_t>& vertices, grade_t& prev_grade, grade_t& max_grade, std::vector<Simplex<index_t>>& low_simplices, std::vector<Simplex<index_t>>& mid_simplices, std::vector<Simplex<index_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     points {std::vector<std::vector<input_t>>} -- a list of points in space.
     metrics {std::vector<Metric>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter>} -- a list of filter functions on the points that should be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
        return;
    }
    
    for (size_t j = vertices.back() + 1; j < points.size(); j++) {
        grade_t grade =  prev_grade.join(point_grades[j]);
        bool should_add = grade.leq_poset(max_grade);
        
        if(should_add){
            vertices.push_back(j);
            build_VR_subcomplex_grades(points, point_grades, vertices, grade, max_grade, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}



void build_VR_subcomplex_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<std::vector<input_t>>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     distance_matrices {std::vector<std::vector<std::vector<input_t>>>} -- a list of distance matrices.
     filters {std::vector<std::vector<input_t>>} -- a list of evaluations of filter functions on the points that will be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    size_t n_points = distance_matrices[0].size();
    
    for (size_t j = vertices.back() + 1; j < n_points; j++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        for (size_t k = 0; k < vertices.size(); k++) {
            size_t p = vertices[k];
            for(size_t metric_index = 0; metric_index<distance_matrices.size(); metric_index++){
                input_t d = distance_matrices[metric_index][p][j];
                if(d < max_metric_values[metric_index]){
                    curr_values[metric_index] = d < curr_values[metric_index] ? curr_values[metric_index] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(!should_add){
                break;
            }
        }
        
        if(should_add){
            for(size_t filter_index=0; filter_index<filters.size(); filter_index++){
                curr_values[distance_matrices.size()+filter_index] = filters[filter_index][j] < curr_values[distance_matrices.size()+filter_index] ? curr_values[distance_matrices.size()+filter_index] : filters[filter_index][j];
            }
            vertices.push_back(j);
            build_VR_subcomplex_dm(distance_matrices, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}


std::pair<Matrix, Matrix> compute_boundary_matrices(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<input_t>& max_metric_values, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    if(points.size() > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        for (size_t i = 0; i < points.size(); i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            std::vector<input_t> curr_values;
            for(size_t metric_index=0; metric_index < metrics.size(); metric_index++){
                curr_values.push_back(0);
            }
            for(size_t filter_index=0; filter_index < filters.size(); filter_index++){
                curr_values.push_back((*filters[filter_index]).eval(points[i]));
            }
            
            std::vector<std::vector<input_t>> points_bak = points;
            
            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices, mid_simplices, low_simplices).first;
        add_boundary_columns<input_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}


std::vector<std::vector<input_t>> compute_distance_matrix(std::vector<std::vector<input_t>>& point_cloud, Metric* metric) {
    /* Compressed distance matrix */
    std::vector<std::vector<input_t>> values;
    for (size_t i = 0; i < point_cloud.size(); i++) {
        std::vector<input_t> row;
        for (size_t j = 0; j <= i ; j++) {
            // input_t val = ((int)(1*(*metric).eval(point_cloud[i], point_cloud[j])));
            input_t val = (((*metric).eval(point_cloud[i], point_cloud[j])));
            row.push_back(val);
        }
        values.push_back(row);
    }
    return values;
}


std::vector<std::vector<std::vector<input_t>>> get_trajectory_dms(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric){
    /* Computing the distance matrix at each time-step of the trajectories */
    std::vector<std::vector<std::vector<input_t>>> distance_matrices;

    for (size_t timestep = 0; timestep < trajectories[0].size(); timestep++) {
        std::vector<std::vector<input_t>> point_cloud;
        for (size_t point = 0; point < trajectories.size(); point++) {
            point_cloud.push_back(trajectories[point][timestep]);
        }
        distance_matrices.push_back(compute_distance_matrix(point_cloud, metric));
    }
    return distance_matrices;
}



std::vector<Simplex<input_t>> calculate_birth_grades2(std::vector<std::vector<std::vector<input_t>>> trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices) {
    std::vector<Simplex<input_t>> new_simplices;
    for(size_t simplex_index=0; simplex_index<simplices.size(); simplex_index++){
        input_t prev_value = 0;
        for (size_t i=1; i< simplices[simplex_index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                prev_value = trajectory_dms[0][i][j] < prev_value ? prev_value : trajectory_dms[0][i][j];
            }
        }
        input_t curr_value = 0;
        for (size_t i=1; i< simplices[simplex_index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                curr_value = trajectory_dms[1][i][j] < curr_value ? curr_value : trajectory_dms[1][i][j];
            }
        }
        if ( prev_value <= curr_value && prev_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(0);
            grade.push_back(prev_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[simplex_index].vertices), std::vector<input_t>(grade), simplex_index));
        }
        
        for (size_t time_simplex_index=2; time_simplex_index < trajectory_dms.size(); time_simplex_index++){
            input_t next_value = 0;
            for (size_t i=1; i< simplices[simplex_index].vertices.size(); i++) {
                for (size_t j=0; j<i ; j++) {
                    next_value = trajectory_dms[time_simplex_index][i][j] < next_value ? next_value : trajectory_dms[time_simplex_index][i][j];
                }
            }
            
            if ( curr_value <= prev_value && curr_value <= next_value && curr_value < max_metric_value ) {
                std::vector<input_t> grade;
                grade.push_back(trajectory_dms.size()-time_simplex_index);
                grade.push_back(time_simplex_index-1);
                grade.push_back(curr_value);
                new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[simplex_index].vertices), std::vector<input_t>(grade), simplex_index));
            }
            
            prev_value = curr_value;
            curr_value = next_value;
        }
        
        if ( curr_value <= prev_value && curr_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(0);
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(curr_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[simplex_index].vertices), std::vector<input_t>(grade), simplex_index));
        }
        
    }
    return new_simplices;
}



input_t interlevel_rips_radius(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, std::vector<size_t>& vertices, size_t start_time, size_t end_time) {

    input_t radius = 0.0;

    for (size_t i=1; i< vertices.size(); i++) {
        for (size_t j=0; j<i ; j++) {
            
            input_t edge_length = trajectory_dms[start_time][vertices[i]][vertices[j]];

            for (size_t t = start_time+1; t <= end_time; t++){
                input_t new_edge_length = trajectory_dms[t][vertices[i]][vertices[j]];
                edge_length = new_edge_length < edge_length ? new_edge_length : edge_length;
            }
            
            radius = edge_length > radius ? edge_length : radius;
        }
    }
    
    return radius;
}



std::vector<Simplex<input_t>> calculate_birth_grades(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices) {
    std::vector<Simplex<input_t>> new_simplices;

    size_t num_time_steps = trajectory_dms.size();

    // radius_map[t0, t1] records the radius of the current simplex under the interlevel-rips filtration at interval [t0,t1]
    std::vector<std::vector<input_t>> radius_map(num_time_steps, std::vector<input_t>(num_time_steps, -1));
    
    int previous_length = 0;

    // We'll go through each simplex one by one and add all the relevant grades in lex order
    // The final result needs to be in reverse lex, so we will reverse it
    for(size_t simplex_index=0; simplex_index<simplices.size(); simplex_index++){
        

        std::vector<size_t> vertices = simplices[simplex_index].vertices;

        // inverted_start_time = num_time_steps - 1 - start_time
        // so start_time = num_time_steps - 1 - inverted_start_time
        for (size_t inverted_start_time=0; inverted_start_time < num_time_steps; inverted_start_time++){
            
            size_t start_time = num_time_steps - 1 - inverted_start_time;

            for (size_t end_time=start_time; end_time < num_time_steps; end_time++) {

                input_t radius = interlevel_rips_radius(trajectory_dms, vertices, start_time, end_time);
                
                radius_map[start_time][end_time] = radius;
                
                // Check the two previous values in the filtration if they exist
                bool inserted_already = (start_time < end_time) && ( (radius_map[start_time + 1][end_time] <= radius) ||  (radius_map[start_time][end_time-1] <= radius)) ;
                
                
                // std::cout << "Processing grade" << start_time << end_time << radius << std::endl;

                // if (inserted_already) {
                //     std::cout << "INSERTED ALREADY" << std::endl;
                // }
                    
                if ((!inserted_already) && (radius < max_metric_value)){
                    // std::cout << "Inserting at grade" << start_time << end_time << radius << std::endl;

                    std::vector<input_t> grade;
                    grade.push_back(inverted_start_time);
                    grade.push_back(end_time);
                    grade.push_back(radius);
                    new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(vertices), std::vector<input_t>(grade), simplex_index));
                }
                
            } 

        }

         // Reset radius_map
        for (size_t i = 0; i < num_time_steps; i++) {
            std::fill(radius_map[i].begin(), radius_map[i].end(), -1);
        }
        
        // Reverse the newly added simplices
        std::reverse(new_simplices.begin() + previous_length, new_simplices.end());

        previous_length = new_simplices.size();
       
    }
    return new_simplices;
}




std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>> compute_boundary_matrices_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric, input_t max_metric_value, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     trajectories {std::vector<std::vector<std::vector<input_t>>>} -- a time-varying point cloud, i.e. a list of trajectories
     metric {Metric} -- a metric on the points that should be used to contruct the complex
     max_metric_value {input_t} -- the maximum allowed distance for the metric
     hom_dim {int} -- the dimension of the homology to be computed
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    // The discretisation map on the parameters
    std::vector<std::vector<input_t>> index_value_lists;

    if(trajectories.size() > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        // Distance matrix of the point cloud at each timestep
        std::vector<std::vector<std::vector<input_t>>> trajectory_dms = get_trajectory_dms(trajectories, metric);
        
        for (size_t vertex = 0; vertex < trajectories.size(); vertex++) {
            std::vector<size_t> vertices;
            vertices.push_back(vertex);
            
            std::vector<input_t> curr_values;
            curr_values.push_back(0);
            
            
            build_VR_subcomplex_spatiotemporal(trajectory_dms, metric, vertices, curr_values, max_metric_value, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        
        std::vector<Simplex<input_t>> high_simplices_extended = calculate_birth_grades(trajectory_dms, max_metric_value, high_simplices);
        
        std::vector<Simplex<input_t>> mid_simplices_extended = calculate_birth_grades(trajectory_dms, max_metric_value, mid_simplices);
        
        std::pair<std::vector<hash_map<input_t, index_t>>,std::vector<std::vector<input_t>>> grade_map = compute_grade_map<input_t>(high_simplices_extended, mid_simplices_extended, low_simplices);
        add_boundary_columns_multicritical<input_t>(high_simplices_extended, mid_simplices_extended, high_matrix, grade_map.first);
        add_boundary_columns<input_t>(mid_simplices_extended, low_simplices, low_matrix, grade_map.first);
        index_value_lists = grade_map.second;
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>(high_matrix, low_matrix, index_value_lists);
}



std::pair<Matrix, Matrix> compute_boundary_matrices_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& grades, grade_t& max_grade, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    if(points.size() > 0){
        std::vector<Simplex<index_t>> low_simplices;
        std::vector<Simplex<index_t>> mid_simplices;
        std::vector<Simplex<index_t>> high_simplices;
        
        for (size_t i = 0; i < points.size(); i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            grade_t grade = grades[i];
            
            std::vector<std::vector<input_t>> points_bak = points;
            
            build_VR_subcomplex_grades(points, grades, vertices, grade, max_grade, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<index_t, index_t>> grade_map = compute_grade_map<index_t>(high_simplices, mid_simplices, low_simplices).first;
        add_boundary_columns<index_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<index_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}


std::pair<Matrix, Matrix> compute_boundary_matrices_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     distance_matrices {std::vector<std::vector<std::vector<input_t>>>} -- a list of distance matrices.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     filters {std::vector<std::vector<input_t>>} -- a list of evaluations of filter functions on the points that will be used to contruct the complex.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    size_t n_points = distance_matrices[0].size();
    
    if(n_points > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        for (size_t i = 0; i < n_points; i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            std::vector<input_t> curr_values;
            for(size_t metric_index=0; metric_index < distance_matrices.size(); metric_index++){
                curr_values.push_back(0);
            }
            for(size_t filter_index=0; filter_index < filters.size(); filter_index++){
                curr_values.push_back(filters[filter_index][i]);
            }
            
            build_VR_subcomplex_dm(distance_matrices, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices, mid_simplices, low_simplices).first;
        add_boundary_columns<input_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}



void verify_kernel(Matrix& kernel, Matrix& map){
    std::cout << "Starting to verify kernel of map...";
    for(size_t i=0; i<kernel.size(); i++){
        SignatureColumn working_column(kernel[i].grade, i, kernel[i]);
        SignatureColumn result(kernel[i].grade, i);
        while(true){
            index_t pivot = working_column.pop_pivot().get_index();
            if(pivot != -1){
                result.plus(map[pivot]);
            }else{
                break;
            }
        }
        index_t pivot = result.get_pivot().get_index();
        if(pivot!=-1){
            std::cerr << "Column is not in kernel of map";
            throw std::runtime_error("Kernel verification failed.");
        }
    }
    std::cout << "Finished verifying kernel of map.";
    std::cout << "Success.";
}


std::vector<std::vector<input_t>> read_point_cloud(std::istream& input_stream, int precision) {
    /** Reads a point cloud from input.
     
     The input should be formatted as: p_11 p_12 ... p_1m
                                       p_21 p_22 ... p_2m etc
     
     Arguments:
            std::istream& input_stream -- input stream of the point cloud
            int precision -- the number of bits of precision TODO: refactor the precision.
     
     Returns:
            std::vector<std::vector<input_t>> -- a list of points represented as vectors with entries of type 'input_t'.
     **/
    std::cout << "Starting to read point cloud...";
    std::vector<std::vector<input_t>> points;
    std::string line;
    input_t value;
    while (std::getline(input_stream, line)) {
        std::vector<input_t> point;
        std::istringstream s(line);
        while (s >> value) {
            point.push_back(input_t(value*pow(10, precision)));
            s.ignore();
        }
        if (!point.empty()){
            points.push_back(point);
        }
    }
    std::cout << "Finished reading points: " << points.size();
    return points;
}



void apply_diagonal_grade_transform(Matrix& columns) {
    for (size_t col_index=0; col_index<columns.size(); col_index++) {
        for (size_t i=0; i<columns[col_index].grade.size()-1; i++) {
            index_t updated_value = columns[col_index].grade[i] - columns[col_index].grade[columns[col_index].grade.size()-1];
            columns[col_index].grade[i] = updated_value >= 0 ? updated_value : 0;
        }
    }
}


void write_trajectories_to_file(std::vector<std::vector<std::vector<input_t>>> trajectories, std::ofstream& output_stream) {
    
    for ( auto& trajectory : trajectories ) {
        for ( auto& point : trajectory ) {
            for ( auto& coord : point ) {
                output_stream << coord << ",";
            }
            output_stream << "\n";
        }
        output_stream << "::\n";
    }
    
}