#include "spatiotemporal.h"

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


                if ((!inserted_already) && (radius < max_metric_value)){

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
