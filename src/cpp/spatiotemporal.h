#ifndef MPH_SPATIOTEMPORAL_INCLUDED
#define MPH_SPATIOTEMPORAL_INCLUDED

#include <fstream>
#include <tuple>
#include <vector>

#include "matrix.h"
#include "rips.h"
#include "utils.h"

void build_VR_subcomplex_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, Metric* metric, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, input_t& max_metric_value, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim);

std::vector<std::vector<std::vector<input_t>>> get_trajectory_dms(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric);

std::vector<Simplex<input_t>> calculate_birth_grades2(std::vector<std::vector<std::vector<input_t>>> trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices);

input_t interlevel_rips_radius(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, std::vector<size_t>& vertices, size_t start_time, size_t end_time);

std::vector<Simplex<input_t>> calculate_birth_grades(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices);

std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>> compute_boundary_matrices_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric, input_t max_metric_value, int hom_dim);

void write_trajectories_to_file(std::vector<std::vector<std::vector<input_t>>> trajectories, std::ofstream& output_stream);

#endif // MPH_SPATIOTEMPORAL_INCLUDED
