#include "bindings.h"

#include <utility>
#include <vector>

#include "grade.h"
#include "groebner.h"
#include "landscapes.h"
#include "matrix.h"
#include "presentation.h"
#include "rips.h"
#include "signatureColumn.h"
#include "spatiotemporal.h"
#include "utils.h"

Matrix translateInputMatrix(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& column_grades){
    Matrix M;
    for(size_t i=0; i<matrix.size(); i++){
        grade_t grade;
        for(size_t j=0; j<column_grades[i].size(); j++){
            grade.push_back(column_grades[i][j]);
        }
        SignatureColumn column(grade, i);
        for(size_t j=0; j<matrix[i].size(); j++){
            if(matrix[i][j] != 0){
                column.push(column_entry_t(1, j)); // TODO: implement support for other modulus than Z/2Z.
            }
        }
        column.syzygy.push(column_entry_t(1, i));
        M.push_back(column);
    }
    return M;
}

GradedMatrix translateOutputMatrix(Matrix& matrix){
    GradedMatrix graded_matrix;
    for(SignatureColumn column : matrix){
        std::vector<std::pair<int, int>> sparse_column;
        while(!column.empty()){
            column_entry_t entry = column.pop_pivot();
            if(entry.get_index() != -1){
                sparse_column.push_back(std::pair<int, int>(entry.get_value(), entry.get_index()));
            }
        }
        graded_matrix.matrix.push_back(sparse_column);
        std::vector<int> column_grade;
        for(size_t grade_index=0; grade_index<column.grade.size(); grade_index++){
            column_grade.push_back((int)column.grade[grade_index]);
        }
        graded_matrix.column_grades.push_back(column_grade);
    }
    return graded_matrix;
}

GradedMatrix compute_kernel(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    Matrix kernel;
    if(column_grades.size()>0 && column_grades[0].size()==2){
        kernel = computeKernel_2p(M);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
        kernel = gbs.first;
    }
    GradedMatrix Ker = translateOutputMatrix(kernel);
    Ker.row_grades = column_grades;
    return Ker;
}

std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
    GradedMatrix Im = translateOutputMatrix(gbs.first);
    GradedMatrix Ker = translateOutputMatrix(gbs.second);
    Ker.row_grades = column_grades;
    Im.row_grades = row_grades;
    return std::pair<GradedMatrix, GradedMatrix>(Im, Ker);
}

GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l){
    Matrix M_h = translateInputMatrix(high_matrix, column_grades_h);
    Matrix M_l = translateInputMatrix(low_matrix, column_grades_l);
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation_output;
    if(M_l.size()>0 && M_l[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(M_l);
        presentation_output = compute_presentation_2p(M_h, kernel);
    }else{
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(M_h, M_l);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation_output.first);
    for(std::pair<size_t, grade_t> entry : presentation_output.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation_dm(DistanceMatrices& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices_dm(distance_matrices, max_metric_values, filters, hom_dim);
    verify_kernel(boundary_matrices.first, boundary_matrices.second);
    for(size_t i=0; i<boundary_matrices.second.size(); i++){
        boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
    }
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases(boundary_matrices.second);
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(boundary_matrices.first, gbs.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation(PointCloud& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim){
    std::vector<Metric*> metrics;
    for(size_t i=0; i<_metrics.size(); i++){
        metrics.push_back(parse_metric(_metrics[i]));
    }
    std::vector<Filter*> filters;
    for(size_t i=0; i<_filters.size(); i++){
        filters.push_back(parse_filter(_filters[i]));
    }
    std::vector<input_t> max_metric_values;
    for(size_t i=0; i<_metrics.size(); i++){
        max_metric_values.push_back(_max_metric_values[i]);
    }
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(_points, metrics, filters, max_metric_values, hom_dim);
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
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(boundary_matrices.first, boundary_matrices.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}


namespace {

std::vector<int> grade_indices_to_int_vector(const grade_t& grade) {
    std::vector<int> out;
    for ( auto& e : grade ) {
        out.push_back((int)e);
    }
    return out;
}

std::pair<std::vector<int>, std::vector<std::vector<int>>> translate_compressed_slice(const std::pair<grade_t, std::vector<grade_t>>& slice) {
    std::vector<std::vector<int>> syzygies;
    for ( auto& s : slice.second ) {
        syzygies.push_back(grade_indices_to_int_vector(s));
    }
    return {grade_indices_to_int_vector(slice.first), syzygies};
}


PythonCompressedLandscape boundary_matrices_to_compressed_landscape(
        Matrix& high_matrix,
        Matrix& low_matrix,
        std::vector<std::vector<input_t>>& index_value_lists) {
    /* Shared tail of the spatiotemporal entry points: minimal presentation -> compressed landscape
     -> Python-facing struct. index_value_lists translates integer critical values to underlying
     filtration values and is forwarded as-is. */
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation =
        computeMinimalPresentation_3p(high_matrix, low_matrix);
    CompressedLandscape compressed_landscape =
        computeCompressedLandscape(minimal_presentation.first, minimal_presentation.second);

    PythonCompressedLandscape result;
    for (auto& slice : compressed_landscape) {
        result.pairings.push_back(translate_compressed_slice(slice));
    }
    result.index_value_lists = index_value_lists;
    return result;
}

}  // namespace

PythonCompressedLandscape landscapes_spatiotemporal(Trajectories& trajectories, input_t max_metric_value, int hom_dim){
    /* Computes the spatiotemporal compressed landscape for a time-varying point cloud.

     trajectories {Trajectories} -- a time-varying point cloud of shape n_agents x n_timepoints x n_spatial_dim
     max_metric_value {input_t} -- the maximum allowed metric value for inclusion in the Vietoris-Rips complex.
     hom_dim {int} -- homology dimesion

     Returns:
     PythonCompressedLandscape -- a list of (birth, syzygies) pairs whose grade components are integer indices,
     plus index_value_lists, the per-axis lookup table that maps those indices back to real filtration values.
     */
    SquaredEuclideanMetric metric;
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) =
        compute_boundary_matrices_spatiotemporal(trajectories, &metric, max_metric_value, hom_dim);
    return boundary_matrices_to_compressed_landscape(high_matrix, low_matrix, index_value_lists);
}

PythonCompressedLandscape landscapes_spatiotemporal_tree(PositionsPerTime& positions_per_t,
                                                         std::vector<std::vector<int>>& parents_per_t,
                                                         input_t max_metric_value,
                                                         int hom_dim){
    /* Computes the spatiotemporal compressed landscape for a tree-structured time-varying dataset
     (e.g. a branching/mitosis process), where the vertex set is the agents alive at the final
     timestep.

     positions_per_t -- per-timestep positions of agents alive at each timestep, indexed
        [timestep][agent_at_timestep][spatial_dim].
     parents_per_t -- per-timestep parent indices into positions_per_t[t-1], or -1 if the agent's
        lineage starts at this timestep. parents_per_t[0] is unused.
     max_metric_value -- the maximum allowed metric value for inclusion in the Vietoris-Rips
        complex (also used as the "edge missing" sentinel for pairs whose lineages weren't both
        alive at a given time).
     hom_dim -- homology dimension.

     Returns:
     PythonCompressedLandscape -- same shape as landscapes_spatiotemporal output.
     */
    SquaredEuclideanMetric metric;
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) =
        compute_boundary_matrices_spatiotemporal_tree(positions_per_t, parents_per_t, &metric,
                                                      max_metric_value, hom_dim);
    return boundary_matrices_to_compressed_landscape(high_matrix, low_matrix, index_value_lists);
}

PythonLandscape landscapes_spatiotemporal_naive(Trajectories& trajectories, input_t& max_metric_value, int hom_dim, int landscape_dim){
    Metric* m = new SquaredEuclideanMetric();
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) = compute_boundary_matrices_spatiotemporal(trajectories, m, max_metric_value, hom_dim);

    grade_t v_max = low_matrix[0].grade;
    for ( auto& column : low_matrix ) {
        v_max = column.grade.join(v_max);
    }
    for ( auto& column : high_matrix ) {
        v_max = column.grade.join(v_max);
    }

    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(high_matrix, low_matrix);

    Landscape landscape = computeLandscapeNaive(minimal_presentation.first, minimal_presentation.second, v_max, landscape_dim);

    PythonLandscape python_landscape;
    for ( auto& slice : landscape ) {
        std::vector<int> g = grade_indices_to_int_vector(slice.first);
        int f = (int)slice.second;
        python_landscape.landscape.push_back(std::pair<std::vector<int>, int>(g, f));
    }
    return python_landscape;
}
