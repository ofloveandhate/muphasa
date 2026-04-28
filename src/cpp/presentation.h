#ifndef MPH_PRESENTATION_INCLUDED
#define MPH_PRESENTATION_INCLUDED

#include <set>
#include <utility>
#include <vector>

#include "grade.h"
#include "groebner.h"
#include "matrix.h"
#include "signatureColumn.h"
#include "utils.h"

std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation(Matrix& image_columns, Matrix& columns, bool debug=true);
std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation_3p(Matrix& image_columns, Matrix& columns, bool debug=true);

hash_map<size_t, size_t> compute_local_pairs(Matrix& columns, hash_map<size_t, grade_t>& row_grades);
Matrix compute_global_columns(Matrix& columns, hash_map<size_t, size_t>& positive_pairs, std::set<size_t>& negative_rows);
std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_2p(Matrix& image_generators, Matrix& generating_set_kernel);

std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_schreyer(Matrix& image_generators, Matrix& kernel_generators, bool debug=false);

#endif // MPH_PRESENTATION_INCLUDED
