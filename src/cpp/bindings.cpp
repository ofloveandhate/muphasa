#include "bindings.h"

#include <algorithm>
#include <limits>
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
    /* Computes the spatiotemporal compressed landscape for a tree structured time-varying dataset
     (e.g. a branching/mitosis process). The vertex set is every leaf of the lineage tree (every
     node with no children, whether it dies at the final timestep or earlier).

     positions_per_t: per-timestep positions of agents alive at each timestep, indexed
        [timestep][agent_at_timestep][spatial_dim].
     parents_per_t: per-timestep parent indices into positions_per_t[t-1], or -1 if the agent's
        lineage starts at this timestep. parents_per_t[0] is unused.
     max_metric_value: the maximum allowed metric value for inclusion in the Vietoris-Rips
        complex (also used as the "edge missing" sentinel for pairs whose lineages weren't both
        alive at a given time, including times before birth or after death).
     hom_dim: homology dimension.

     Returns:
     PythonCompressedLandscape: same shape as landscapes_spatiotemporal output.
     */
    SquaredEuclideanMetric metric;
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) =
        compute_boundary_matrices_spatiotemporal_tree(positions_per_t, parents_per_t, &metric,
                                                      max_metric_value, hom_dim);
    return boundary_matrices_to_compressed_landscape(high_matrix, low_matrix, index_value_lists);
}

namespace {

bool float_grade_leq_poset(const std::vector<input_t>& a, const std::vector<input_t>& b) {
    /* Componentwise <= on full grades — matches Grade.__le__ in _mph.py. */
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

bool float_grade_prefix_leq_poset(const std::vector<input_t>& a, const std::vector<input_t>& b,
                                  size_t prefix_len) {
    /* Componentwise <= on the first prefix_len components — matches Python's syzygy[:-1] <= grade[:-1]
     when prefix_len == grade.size() - 1. */
    for (size_t i = 0; i < prefix_len; ++i) {
        if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

bool float_grade_lt_colex(const std::vector<input_t>& a, const std::vector<input_t>& b) {
    /* Strict colex order: compare components from last to first. Mirrors Grade.colex_compare
     in _mph.py. */
    for (size_t i = a.size(); i-- > 0; ) {
        if (a[i] < b[i]) {
            return true;
        }
        if (a[i] > b[i]) {
            return false;
        }
    }
    return false;
}

bool int_grade_leq_poset(const std::vector<int>& a, const std::vector<int>& b) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

bool int_grade_prefix_leq_poset(const std::vector<int>& a, const std::vector<int>& b, size_t prefix_len) {
    for (size_t i = 0; i < prefix_len; ++i) {
        if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

bool int_grade_lt_colex(const std::vector<int>& a, const std::vector<int>& b) {
    for (size_t i = a.size(); i-- > 0; ) {
        if (a[i] < b[i]) {
            return true;
        }
        if (a[i] > b[i]) {
            return false;
        }
    }
    return false;
}

}  // namespace

std::vector<std::vector<input_t>> eval_sparse_landscape_batch(
        std::vector<std::pair<std::vector<input_t>, std::vector<std::vector<input_t>>>>& pairings,
        std::vector<std::vector<input_t>>& grades,
        int k_max) {
    /* Batch-evaluate a SparseLandscape on N grades for layers 1..k_max, mirroring the
     Python SparseLandscape.eval algorithm but in C++. Pairings are passed in with float
     grades (axis-2 already sqrt'd into Euclidean distance by the Python caller).

     pairings -- (birth, syzygies) entries; syzygies are sorted by colex in-place once
        per batch call so the per-grade inner loop can break on the first prefix-fitting one.
     grades -- N input grades, shape (N, d). All grades and births must share dimension d.
     k_max -- top landscape layer to evaluate; layer k goes to out[i][k-1].

     Returns:
     out -- shape (N, k_max). out[i][k-1] is layer k of the landscape evaluated at grades[i],
        or 0 if fewer than k bars are alive there. */

    if (k_max <= 0) {
        return std::vector<std::vector<input_t>>(grades.size(), std::vector<input_t>());
    }
    std::vector<std::vector<input_t>> out(grades.size(), std::vector<input_t>(k_max, 0.0));
    if (grades.empty() || pairings.empty()) {
        return out;
    }

    for (auto& p : pairings) {
        std::sort(p.second.begin(), p.second.end(), float_grade_lt_colex);
    }

    std::vector<input_t> diffs;
    for (size_t gi = 0; gi < grades.size(); ++gi) {
        const std::vector<input_t>& grade = grades[gi];
        const size_t r = grade.size() - 1;
        diffs.clear();

        for (const auto& entry : pairings) {
            const std::vector<input_t>& birth = entry.first;
            if (!float_grade_leq_poset(birth, grade)) {
                continue;
            }
            const input_t low_distance = grade[r] - birth[r];
            input_t high_distance = std::numeric_limits<input_t>::infinity();
            for (const auto& syzygy : entry.second) {
                if (float_grade_prefix_leq_poset(syzygy, grade, r)) {
                    high_distance = syzygy[r] - grade[r];
                    break;
                }
            }
            input_t contribution = std::min(low_distance, high_distance);
            if (contribution < 0) {
                contribution = 0;
            }
            diffs.push_back(contribution);
        }

        // Top-k descending. Full sort is fine: diffs is at most O(#pairings) which is small,
        // and we read out the first k_max entries.
        std::sort(diffs.begin(), diffs.end(), std::greater<input_t>());
        const size_t fill = std::min(diffs.size(), static_cast<size_t>(k_max));
        for (size_t k = 0; k < fill; ++k) {
            out[gi][k] = diffs[k];
        }
    }

    return out;
}

std::vector<int> eval_sparse_landscape_batch_int(
        std::vector<std::pair<std::vector<int>, std::vector<std::vector<int>>>>& pairings,
        std::vector<std::vector<int>>& grades,
        int k) {
    /* Integer-grade analogue of eval_sparse_landscape_batch, evaluating one layer per call.
     Pairings carry integer grades (as produced directly by the spatiotemporal compressed
     landscape pipeline before any index-to-value mapping); query grades are also integer.

     pairings -- (birth, syzygies); syzygies are sorted in-place by colex, once per call.
     grades -- N integer query grades, shape (N, d).
     k -- layer; out[i] = layer-k height at grades[i], 0 if fewer than k bars are alive there.

     Returns: vector<int> of length N. */
    std::vector<int> out(grades.size(), 0);
    if (k <= 0 || grades.empty() || pairings.empty()) {
        return out;
    }
    for (auto& p : pairings) {
        std::sort(p.second.begin(), p.second.end(), int_grade_lt_colex);
    }
    std::vector<int> diffs;
    const int int_inf = std::numeric_limits<int>::max();
    for (size_t gi = 0; gi < grades.size(); ++gi) {
        const std::vector<int>& grade = grades[gi];
        const size_t r = grade.size() - 1;
        diffs.clear();
        for (const auto& entry : pairings) {
            const std::vector<int>& birth = entry.first;
            if (!int_grade_leq_poset(birth, grade)) {
                continue;
            }
            const int low_distance = grade[r] - birth[r];
            int high_distance = int_inf;
            for (const auto& syzygy : entry.second) {
                if (int_grade_prefix_leq_poset(syzygy, grade, r)) {
                    high_distance = syzygy[r] - grade[r];
                    break;
                }
            }
            int contribution = std::min(low_distance, high_distance);
            if (contribution < 0) {
                contribution = 0;
            }
            diffs.push_back(contribution);
        }
        std::sort(diffs.begin(), diffs.end(), std::greater<int>());
        if ((size_t)k <= diffs.size()) {
            out[gi] = diffs[k - 1];
        }
    }
    return out;
}

namespace {

std::vector<int> column_to_int_indices(SignatureColumn column) {
    /* Pop a SignatureColumn into the list of nonzero row indices (Z/2Z). Order is whatever
     pop_pivot produces; translateInputMatrix re-pushes via column.push so order doesn't matter
     for round-trip correctness. */
    std::vector<int> out;
    while (!column.empty()) {
        column_entry_t entry = column.pop_pivot();
        if (entry.get_index() != -1) {
            out.push_back((int)entry.get_index());
        }
    }
    return out;
}

std::vector<grade_t> int_grades_to_grade_t(std::vector<std::vector<int>>& grades) {
    std::vector<grade_t> out;
    out.reserve(grades.size());
    for (auto& g : grades) {
        grade_t gt;
        for (int v : g) {
            gt.push_back((index_t)v);
        }
        out.push_back(gt);
    }
    return out;
}

Matrix translateInputMatrixSparse(std::vector<std::vector<int>>& matrix_sparse,
                                  std::vector<std::vector<int>>& column_grades) {
    /* Sparse counterpart to translateInputMatrix: each column comes in as the list of nonzero
     row indices (Z/2Z) instead of a dense 0/1 row vector. Used by the naive-bench fixture
     where dense encoding would balloon memory for typical minimal presentations. */
    Matrix M;
    for (size_t i = 0; i < matrix_sparse.size(); i++) {
        grade_t grade;
        for (int g : column_grades[i]) {
            grade.push_back(g);
        }
        SignatureColumn column(grade, i);
        for (int j : matrix_sparse[i]) {
            column.push(column_entry_t(1, j));
        }
        column.syzygy.push(column_entry_t(1, i));
        M.push_back(column);
    }
    return M;
}

}  // namespace

PythonNaiveBenchFixture naive_bench_setup_spatiotemporal(Trajectories& trajectories,
                                                        input_t max_metric_value,
                                                        int hom_dim) {
    /* Setup entry point for the naive-vs-eval2 bench. Computes boundary matrices, the minimal
     presentation, and the compressed landscape from a single trajectories input, packaging
     everything into a fixture so subsequent eval calls operate on identical data. The
     presentation is encoded as a sparse integer matrix (column = list of row indices) so it
     can round-trip through Cython into Python. */
    SquaredEuclideanMetric metric;
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) =
        compute_boundary_matrices_spatiotemporal(trajectories, &metric, max_metric_value, hom_dim);

    std::pair<Matrix, std::vector<grade_t>> minimal_presentation =
        computeMinimalPresentation_3p(high_matrix, low_matrix);
    Matrix& presentation = minimal_presentation.first;
    std::vector<grade_t>& row_grades = minimal_presentation.second;

    Matrix presentation_for_compress;
    for (auto& column : presentation) {
        presentation_for_compress.push_back(column);
    }
    CompressedLandscape compressed_landscape =
        computeCompressedLandscape(presentation_for_compress, row_grades);

    PythonNaiveBenchFixture fixture;
    for (auto& slice : compressed_landscape) {
        fixture.pairings.push_back(translate_compressed_slice(slice));
    }
    fixture.index_value_lists = index_value_lists;

    grade_t v_max;
    if (!row_grades.empty()) {
        v_max = row_grades[0];
    } else if (!presentation.empty()) {
        v_max = presentation[0].grade;
    }
    for (auto& column : presentation) {
        fixture.presentation_matrix.push_back(column_to_int_indices(column));
        fixture.presentation_column_grades.push_back(grade_indices_to_int_vector(column.grade));
        v_max = v_max.join(column.grade);
    }
    for (auto& g : row_grades) {
        fixture.presentation_row_grades.push_back(grade_indices_to_int_vector(g));
        v_max = v_max.join(g);
    }
    fixture.v_max = grade_indices_to_int_vector(v_max);
    return fixture;
}

std::vector<int> eval_naive_rank_int(std::vector<std::vector<int>>& presentation_matrix,
                                     std::vector<std::vector<int>>& presentation_column_grades,
                                     std::vector<std::vector<int>>& presentation_row_grades,
                                     std::vector<std::vector<int>>& query_grades,
                                     int landscape_dim) {
    /* Reconstructs the C++ Matrix from the integer presentation passed in by the Python caller,
     then runs computeLandscapeNaiveRankAtGrades for each query grade at layer landscape_dim.
     Reconstruction is one pass over the presentation; the per-grade rank reduction dominates. */
    Matrix M = translateInputMatrixSparse(presentation_matrix, presentation_column_grades);
    std::vector<grade_t> row_grades = int_grades_to_grade_t(presentation_row_grades);
    std::vector<grade_t> queries = int_grades_to_grade_t(query_grades);
    std::vector<size_t> heights = computeLandscapeNaiveRankAtGrades(M, row_grades, queries, landscape_dim);
    std::vector<int> out;
    out.reserve(heights.size());
    for (auto h : heights) {
        out.push_back((int)h);
    }
    return out;
}

std::vector<int> eval_naive_cached_int(std::vector<std::vector<int>>& presentation_matrix,
                                       std::vector<std::vector<int>>& presentation_column_grades,
                                       std::vector<std::vector<int>>& presentation_row_grades,
                                       std::vector<std::vector<int>>& query_grades,
                                       int landscape_dim) {
    /* Same as eval_naive_rank_int but calling the H-basis-cached variant. The cache lives only
     for the duration of this call (built fresh inside computeLandscapeNaiveCachedAtGrades) so
     repeat-bench invocations don't share state. */
    Matrix M = translateInputMatrixSparse(presentation_matrix, presentation_column_grades);
    std::vector<grade_t> row_grades = int_grades_to_grade_t(presentation_row_grades);
    std::vector<grade_t> queries = int_grades_to_grade_t(query_grades);
    std::vector<size_t> heights = computeLandscapeNaiveCachedAtGrades(M, row_grades, queries, landscape_dim);
    std::vector<int> out;
    out.reserve(heights.size());
    for (auto h : heights) {
        out.push_back((int)h);
    }
    return out;
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
