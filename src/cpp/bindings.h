#ifndef MPH_BINDINGS_INCLUDED
#define MPH_BINDINGS_INCLUDED

#include <utility>
#include <vector>

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

GradedMatrix presentation(std::vector<std::vector<input_t>>& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim);
GradedMatrix presentation_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim);
GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l);
std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades);
PythonCompressedLandscape landscapes_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, input_t max_metric_value, int hom_dim);

#endif // MPH_BINDINGS_INCLUDED
