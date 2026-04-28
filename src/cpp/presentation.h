#ifndef MPH_PRESENTATION_INCLUDED
#define MPH_PRESENTATION_INCLUDED

#include <utility>
#include <vector>

#include "grade.h"
#include "matrix.h"
#include "signatureColumn.h"

bool cmpID(const signature_t& lhs, const signature_t& rhs);
void insert_sorted(std::vector<signature_t>& cont, signature_t value);

Matrix compute_minimal_generating_set(Matrix& generators);

std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation_3p(Matrix& image_columns, Matrix& columns, bool debug=true);

#endif // MPH_PRESENTATION_INCLUDED
