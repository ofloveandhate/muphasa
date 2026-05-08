#ifndef MPH_BINDINGS_INCLUDED
#define MPH_BINDINGS_INCLUDED

#include <utility>
#include <vector>

#include "rips.h"
#include "spatiotemporal.h"
#include "utils.h"

struct GradedMatrix{
    std::vector< std::vector<std::pair<int, int>> > matrix;
    std::vector< std::vector<int> > column_grades;
    std::vector< std::vector<int> > row_grades;
};

struct PythonLandscape{
    std::vector<std::pair<std::vector<int>, size_t>> landscape;
};

struct PythonCompressedLandscape{
    std::vector<std::pair<std::vector<int>, std::vector<std::vector<int>>>> pairings;
    std::vector<std::vector<input_t>> index_value_lists;
};

struct PythonNaiveBenchFixture{
    /* Shared fixture for the eval2-vs-naive bench. From a single trajectories input we build
     the minimal presentation once and expose it alongside the compressed-landscape integer
     pairings, so both algorithms run on identical data. The presentation is stored as a sparse
     integer matrix (each column = list of row indices, Z/2Z coefficients) so it can round-trip
     through Cython. */
    std::vector<std::pair<std::vector<int>, std::vector<std::vector<int>>>> pairings;
    std::vector<std::vector<input_t>> index_value_lists;
    std::vector<std::vector<int>> presentation_matrix;
    std::vector<std::vector<int>> presentation_column_grades;
    std::vector<std::vector<int>> presentation_row_grades;
    std::vector<int> v_max;
};

GradedMatrix presentation(PointCloud& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim);
GradedMatrix presentation_dm(DistanceMatrices& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim);
GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l);
std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades);
PythonCompressedLandscape landscapes_spatiotemporal(Trajectories& trajectories, input_t max_metric_value, int hom_dim);
PythonCompressedLandscape landscapes_spatiotemporal_tree(PositionsPerTime& positions_per_t, std::vector<std::vector<int>>& parents_per_t, input_t max_metric_value, int hom_dim);

std::vector<std::vector<input_t>> eval_sparse_landscape_batch(
    std::vector<std::pair<std::vector<input_t>, std::vector<std::vector<input_t>>>>& pairings,
    std::vector<std::vector<input_t>>& grades,
    int k_max);

std::vector<int> eval_sparse_landscape_batch_int(
    std::vector<std::pair<std::vector<int>, std::vector<std::vector<int>>>>& pairings,
    std::vector<std::vector<int>>& grades,
    int k);

PythonNaiveBenchFixture naive_bench_setup_spatiotemporal(Trajectories& trajectories,
                                                        input_t max_metric_value,
                                                        int hom_dim);

std::vector<int> eval_naive_rank_int(std::vector<std::vector<int>>& presentation_matrix,
                                     std::vector<std::vector<int>>& presentation_column_grades,
                                     std::vector<std::vector<int>>& presentation_row_grades,
                                     std::vector<std::vector<int>>& query_grades,
                                     int landscape_dim);

std::vector<int> eval_naive_cached_int(std::vector<std::vector<int>>& presentation_matrix,
                                       std::vector<std::vector<int>>& presentation_column_grades,
                                       std::vector<std::vector<int>>& presentation_row_grades,
                                       std::vector<std::vector<int>>& query_grades,
                                       int landscape_dim);

#endif // MPH_BINDINGS_INCLUDED
