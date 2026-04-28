#ifndef MPH_GROEBNER_INCLUDED
#define MPH_GROEBNER_INCLUDED

#include <utility>
#include <vector>

#include "grade.h"
#include "matrix.h"
#include "signatureColumn.h"

bool cmpID(const signature_t& lhs, const signature_t& rhs);
void insert_sorted(std::vector<signature_t>& cont, signature_t value);

Matrix compute_minimal_generating_set(Matrix& generators);
std::pair<Matrix, Matrix> compute_minimal_generating_set2(Matrix generators);
Matrix compute_syzygy_module(Matrix groebner_basis);

std::pair<Matrix, Matrix> computeGroebnerBases(std::vector<SignatureColumn>& columns);
std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt_min(std::vector<SignatureColumn>& columns);
std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt(std::vector<SignatureColumn>& columns);
Matrix computekernel_gradeopt(std::vector<SignatureColumn>& columns);
Matrix ImageGB(std::vector<SignatureColumn>& columns);
std::pair<Matrix, std::vector<std::pair<grade_t, std::vector<size_t>>>> ImageGB_gradeopt_min_sig_basis(std::vector<SignatureColumn>& columns, MultigradedBasis ker_basis);
Matrix buchberger(Matrix& columns);

Matrix computeKernel_2p(std::vector<SignatureColumn>& columns);

#endif // MPH_GROEBNER_INCLUDED
