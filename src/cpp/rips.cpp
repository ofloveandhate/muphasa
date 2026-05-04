#include "rips.h"




Metric* parse_metric(int metric_index){
    switch(metric_index){
        case 0:
            return new SquaredEuclideanMetric();
        default:
            throw std::runtime_error("Failed to parse metric");
    }
}

Filter* parse_filter(int filter_index){
    if(filter_index >= 0){
        return new XFilter(filter_index);
    }
    switch(filter_index){
        case -1:
            return new XFilter(0); //TODO: implement more filters.
        default:
            throw std::runtime_error("Failed to parse metric");
    }
}

void build_VR_subcomplex(PointCloud& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.

     points {PointCloud} -- a list of points in space.
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
                input_t d = (*filters[filter_index]).eval(points[j]);
                curr_values[metrics.size()+filter_index] = d < curr_values[metrics.size()+filter_index] ? curr_values[metrics.size()+filter_index] : d;
            }
            vertices.push_back(j);
            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}




void build_VR_subcomplex_grades(PointCloud& points, std::vector<grade_t>& point_grades, std::vector<size_t>& vertices, grade_t& prev_grade, grade_t& max_grade, std::vector<Simplex<index_t>>& low_simplices, std::vector<Simplex<index_t>>& mid_simplices, std::vector<Simplex<index_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.

     points {PointCloud} -- a list of points in space.
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



void build_VR_subcomplex_dm(DistanceMatrices& distance_matrices, std::vector<std::vector<input_t>>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.

     distance_matrices {DistanceMatrices} -- a list of distance matrices.
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


std::pair<Matrix, Matrix> compute_boundary_matrices(PointCloud& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<input_t>& max_metric_values, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.

     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.

     */
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

            PointCloud points_bak = points;

            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices, mid_simplices, low_simplices).first;
        add_boundary_columns<input_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}


DistanceMatrix compute_distance_matrix(const PointCloud& point_cloud, Metric* metric) {
    /* Computes the lower-triangular pairwise-distance matrix for a point cloud under the given
     metric. Row i has length i+1: entries [i][j] for j <= i.

     point_cloud {PointCloud} -- list of points; each point is a vector of spatial coordinates.
     metric {Metric*} -- metric used to compute pairwise distances.

     Returns:
     DistanceMatrix -- lower-triangular distance matrix, indexed [i][j] with j <= i.
     */
    DistanceMatrix values;
    for (size_t i = 0; i < point_cloud.size(); i++) {
        std::vector<input_t> row;
        for (size_t j = 0; j <= i; j++) {
            row.push_back(metric->eval(point_cloud[i], point_cloud[j]));
        }
        values.push_back(row);
    }
    return values;
}




std::pair<Matrix, Matrix> compute_boundary_matrices_grades(PointCloud& points, std::vector<grade_t>& grades, grade_t& max_grade, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.

     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.

     */
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

            PointCloud points_bak = points;

            build_VR_subcomplex_grades(points, grades, vertices, grade, max_grade, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<index_t, index_t>> grade_map = compute_grade_map<index_t>(high_simplices, mid_simplices, low_simplices).first;
        add_boundary_columns<index_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<index_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}


std::pair<Matrix, Matrix> compute_boundary_matrices_dm(DistanceMatrices& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.

     distance_matrices {DistanceMatrices} -- a list of distance matrices.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     filters {std::vector<std::vector<input_t>>} -- a list of evaluations of filter functions on the points that will be used to contruct the complex.
     hom_dim {int} -- the dimension of the homology to be computed.

     */
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
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}



void verify_kernel(Matrix& kernel, Matrix& map){
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
            throw std::runtime_error("Kernel verification failed.");
        }
    }
}


PointCloud read_point_cloud(std::istream& input_stream, int precision) {
    /** Reads a point cloud from input.

     The input should be formatted as: p_11 p_12 ... p_1m
                                       p_21 p_22 ... p_2m etc

     Arguments:
            std::istream& input_stream -- input stream of the point cloud
            int precision -- the number of bits of precision TODO: refactor the precision.

     Returns:
            PointCloud -- a list of points represented as vectors with entries of type 'input_t'.
     **/
    PointCloud points;
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
