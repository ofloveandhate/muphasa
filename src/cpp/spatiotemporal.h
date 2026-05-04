#ifndef MPH_SPATIOTEMPORAL_INCLUDED
#define MPH_SPATIOTEMPORAL_INCLUDED

#include <tuple>
#include <vector>

#include "matrix.h"
#include "rips.h"
#include "utils.h"


using Trajectories     = std::vector<std::vector<std::vector<input_t>>>;  // [point][t][dim]
using PositionsPerTime = std::vector<std::vector<std::vector<input_t>>>;  // [t][agent][dim]


void build_VR_subcomplex_spatiotemporal(DistanceMatrices& trajectory_dms,
                                        std::vector<size_t>& vertices,
                                        input_t prev_max_dist,
                                        input_t max_metric_value,
                                        std::vector<Simplex<input_t>>& low_simplices,
                                        std::vector<Simplex<input_t>>& mid_simplices,
                                        std::vector<Simplex<input_t>>& high_simplices,
                                        int hom_dim);

DistanceMatrices get_trajectory_dms(Trajectories& trajectories, Metric* metric);

input_t interlevel_rips_radius(const DistanceMatrices& trajectory_dms,
                               const std::vector<size_t>& vertices,
                               size_t start_time,
                               size_t end_time);

std::vector<Simplex<input_t>> calculate_birth_grades(const DistanceMatrices& trajectory_dms,
                                                     input_t max_metric_value,
                                                     const std::vector<Simplex<input_t>>& simplices);

std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal(Trajectories& trajectories,
                                         Metric* metric,
                                         input_t max_metric_value,
                                         int hom_dim);

std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal_dm(DistanceMatrices& trajectory_dms,
                                            input_t max_metric_value,
                                            int hom_dim);

std::vector<std::vector<int>>
build_ancestor_lookup(const std::vector<std::vector<int>>& parents_per_t, size_t n_leaves);

DistanceMatrices tree_positions_to_dms(const PositionsPerTime& positions_per_t,
                                       const std::vector<std::vector<int>>& ancestor_lookup,
                                       Metric* metric,
                                       input_t max_metric_value);

std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal_tree(PositionsPerTime& positions_per_t,
                                              std::vector<std::vector<int>>& parents_per_t,
                                              Metric* metric,
                                              input_t max_metric_value,
                                              int hom_dim);

#endif // MPH_SPATIOTEMPORAL_INCLUDED
