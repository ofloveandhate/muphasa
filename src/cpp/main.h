#ifndef MPH_MAIN_INCLUDED
#define MPH_MAIN_INCLUDED



int c_presentation(float* points, int N, int D, int hom_dim);

#include "utils.h"
#include "grade.h"
#include "column.h"
#include "signatureColumn.h"
#include "matrix.h"
#include "IO.h"
// #include "examples.h"
#include "landscapes.h"





/* Python interface */

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


std::pair<Matrix, Matrix> computeGroebnerBases(std::vector<SignatureColumn>& columns);
std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt_min(std::vector<SignatureColumn>& columns);
std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt(std::vector<SignatureColumn>& columns);
Matrix computekernel_gradeopt(std::vector<SignatureColumn>& columns);
Matrix ImageGB(std::vector<SignatureColumn>& columns);
std::pair<Matrix, std::vector<std::pair<grade_t, std::vector<size_t>>>> ImageGB_gradeopt_min_sig_basis(std::vector<SignatureColumn>& columns, MultigradedBasis ker_basis);
CompressedLandscape computeCompressedLandscape(std::vector<SignatureColumn>& presentation, std::vector<grade_t> row_grades);
DiagLandscape computeSparseDiagLandscape(Matrix& presentation, std::vector<grade_t> row_grades);
Matrix buchberger(Matrix& columns);
Matrix computeKernel_2p(std::vector<SignatureColumn>& columns);
hash_map<size_t, size_t> compute_local_pairs(Matrix& columns, hash_map<size_t, grade_t>& row_grades);
Matrix compute_global_columns(Matrix& columns, hash_map<size_t, size_t>& positive_pairs, std::set<size_t>& negative_rows);
Matrix compute_syzygy_module(Matrix groebner_basis);
std::pair<Matrix, Matrix> compute_minimal_generating_set2(Matrix generators);
std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation(Matrix& image_columns, Matrix& columns, bool debug=true);



// these here are certainly necessary
GradedMatrix presentation(std::vector<std::vector<input_t>>& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim);
GradedMatrix presentation_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim);
GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l);
std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades);
PythonCompressedLandscape landscapes_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, input_t& max_metric_value, int hom_dim);



#endif // MPH_MAIN_INCLUDED
