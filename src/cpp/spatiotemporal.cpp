#include "spatiotemporal.h"

#include <algorithm>

// Sentinel for "no parent / lineage absent" in parent-index arrays and ancestor-lookup tables.
static constexpr int NO_PARENT = -1;

static input_t try_extend_at_timestep(const DistanceMatrix& dm_at_t,
                                      size_t new_vertex,
                                      const std::vector<size_t>& existing,
                                      input_t max_metric_value,
                                      input_t prev_max_dist) {
    input_t curr_max_dist = prev_max_dist;
    for (size_t old_vertex : existing) {
        input_t d = dm_at_t[new_vertex][old_vertex];
        if (d >= max_metric_value) return -1;
        curr_max_dist = std::max(curr_max_dist, d);
    }
    return curr_max_dist;
}

void build_VR_subcomplex_spatiotemporal(DistanceMatrices& trajectory_dms,
                                        std::vector<size_t>& vertices,
                                        input_t prev_max_dist,
                                        input_t max_metric_value,
                                        std::vector<Simplex<input_t>>& low_simplices,
                                        std::vector<Simplex<input_t>>& mid_simplices,
                                        std::vector<Simplex<input_t>>& high_simplices,
                                        int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex for a dynamic metric space
     in the sense of Memoli/Kim (2020).

     trajectory_dms {DistanceMatrices} -- per-timestep distance matrices for the trajectories,
        indexed [timestep][i][j].
     vertices {std::vector<size_t>} -- vertex list of the simplex currently being built
     prev_max_dist {input_t} -- the largest pairwise distance seen so far across all edges of the
        current simplex (i.e. the simplex's Rips radius).
     max_metric_value {input_t} -- maximum allowed distance
     low_simplices, mid_simplices, high_simplices -- output buckets for simplices of dimension
        hom_dim-1, hom_dim, and hom_dim+1 respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     */

    // Push back the current simplex
    if (vertices.size() == hom_dim) {                  // dim hom_dim - 1
        low_simplices.push_back(Simplex<input_t>(vertices, std::vector<input_t>{prev_max_dist}));
    } else if (vertices.size() == hom_dim + 1) {       // dim hom_dim
        mid_simplices.push_back(Simplex<input_t>(vertices, std::vector<input_t>{prev_max_dist}));
    } else if (vertices.size() == hom_dim + 2) {       // dim hom_dim + 1
        high_simplices.push_back(Simplex<input_t>(vertices, std::vector<input_t>{prev_max_dist}));
        return;
    }

    // Recursively try new simplices: for each candidate new vertex, accept the
    // first timestep at which every new edge fits below max_metric_value.
    for (size_t new_vertex = vertices.back() + 1; new_vertex < trajectory_dms[0].size(); new_vertex++) {
        for (size_t time_index = 0; time_index < trajectory_dms.size(); time_index++) {
            input_t curr_max_dist = try_extend_at_timestep(trajectory_dms[time_index], new_vertex,
                                                           vertices, max_metric_value, prev_max_dist);
            if (curr_max_dist >= 0) {
                vertices.push_back(new_vertex);
                build_VR_subcomplex_spatiotemporal(trajectory_dms, vertices, curr_max_dist, max_metric_value,
                                                   low_simplices, mid_simplices, high_simplices, hom_dim);
                vertices.pop_back();
                break;
            }
        }
    }
}

DistanceMatrices get_trajectory_dms(Trajectories& trajectories, Metric* metric){
    /* Computes the pairwise-distance matrix at each timestep of a time-varying point cloud.

     trajectories {Trajectories} -- the time-varying point cloud, indexed [point][timestep][spatial_dim].
     metric {Metric*} -- metric used to compute pairwise distances.

     Returns:
     DistanceMatrices -- per-timestep distance matrices, indexed [timestep][i][j].
     */
    DistanceMatrices distance_matrices;

    if (trajectories.size() == 0) {
        return distance_matrices;
    }

    for (size_t timestep = 0; timestep < trajectories[0].size(); timestep++) {
        PointCloud point_cloud;
        for (size_t point = 0; point < trajectories.size(); point++) {
            point_cloud.push_back(trajectories[point][timestep]);
        }
        distance_matrices.push_back(compute_distance_matrix(point_cloud, metric));
    }
    return distance_matrices;
}


std::vector<std::vector<int>>
build_ancestor_lookup(const std::vector<std::vector<int>>& parents_per_t, size_t n_leaves) {
    /* For each leaf and each timestep, the index of that leaf's ancestor in positions_per_t[t]
     (or -1 if the lineage hadn't started yet at that timestep).

     parents_per_t {std::vector<std::vector<int>>} -- per-timestep parent indices, indexed
        [timestep][agent_at_timestep]. Each entry is an index into positions_per_t[t-1], or -1 if
        the agent's lineage starts at this timestep. parents_per_t[0] is unused (no t-1 to point
        into) but is expected to be present so the outer dimension matches positions_per_t.
     n_leaves {size_t} -- number of agents at the final timestep, which form the vertex set.

     Returns:
     std::vector<std::vector<int>> -- ancestor lookup table indexed [leaf][timestep]; entry is the
        ancestor's index in positions_per_t[t], or -1 if the lineage didn't exist yet at t.
     */
    size_t T = parents_per_t.size();
    std::vector<std::vector<int>> lookup(n_leaves, std::vector<int>(T, NO_PARENT));
    if (T == 0) {
        return lookup;
    }
    for (size_t l = 0; l < n_leaves; l++) {
        lookup[l][T - 1] = (int)l;
        for (size_t t_plus_one = T - 1; t_plus_one > 0; t_plus_one--) {
            int current = lookup[l][t_plus_one];
            if (current == NO_PARENT) {
                break;
            }
            lookup[l][t_plus_one - 1] = parents_per_t[t_plus_one][current];
        }
    }
    return lookup;
}


DistanceMatrices tree_positions_to_dms(const PositionsPerTime& positions_per_t,
                                       const std::vector<std::vector<int>>& ancestor_lookup,
                                       Metric* metric,
                                       input_t max_metric_value) {
    /* Builds per-timestep N_leaves x N_leaves distance matrices for tree-structured input.

     positions_per_t {PositionsPerTime} -- per-timestep positions of agents alive at each timestep,
        indexed [timestep][agent_at_timestep][spatial_dim].
     ancestor_lookup {std::vector<std::vector<int>>} -- output of build_ancestor_lookup, indexed
        [leaf][timestep].
     metric {Metric*} -- metric used to compute pairwise distances between ancestor positions.
     max_metric_value {input_t} -- value used to mark "no edge at this timestep" (when at least
        one of the two leaves' lineages didn't exist yet).

     Returns:
     DistanceMatrices -- per-timestep distance matrices in the same lower-triangular layout as
        compute_distance_matrix: dm[t][i] has length i+1, valid for accesses dm[t][i][j] with
        j <= i. Entries where either lineage was missing at that timestep are set to
        max_metric_value (which the engine treats as "edge missing").
     */
    size_t T = positions_per_t.size();
    size_t n_leaves = ancestor_lookup.size();

    // Lower-triangular layout (matching compute_distance_matrix / get_trajectory_dms): the engine
    // only ever reads dm[bigger][smaller], so the upper triangle is never needed. Halves memory
    // and zero-init traffic vs a full square.
    DistanceMatrices dms(T);
    for (size_t t = 0; t < T; t++) {
        dms[t].resize(n_leaves);
        for (size_t i = 0; i < n_leaves; i++) {
            dms[t][i].resize(i + 1);
            int a_i = ancestor_lookup[i][t];
            for (size_t j = 0; j < i; j++) {
                int a_j = ancestor_lookup[j][t];
                if (a_i == NO_PARENT || a_j == NO_PARENT) {
                    dms[t][i][j] = max_metric_value;
                } else {
                    dms[t][i][j] = metric->eval(positions_per_t[t][a_i], positions_per_t[t][a_j]);
                }
            }
        }
    }
    return dms;
}


input_t interlevel_rips_radius(const DistanceMatrices& trajectory_dms,
                               const std::vector<size_t>& vertices,
                               size_t start_time,
                               size_t end_time) {
    /* The Rips radius of a simplex over the time window [start_time, end_time] of the
     interlevel-rips filtration.

     trajectory_dms {DistanceMatrices} -- per-timestep distance matrices for the trajectories,
        indexed [timestep][i][j].
     vertices {std::vector<size_t>} -- vertex list of the simplex.
     start_time {size_t} -- inclusive start of the time window.
     end_time {size_t} -- inclusive end of the time window.

     Returns:
     input_t -- the simplex's Rips radius over [start_time, end_time].
     */
    input_t radius = 0.0;
    for (size_t i = 1; i < vertices.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            input_t edge_length = trajectory_dms[start_time][vertices[i]][vertices[j]];
            for (size_t t = start_time + 1; t <= end_time; t++) {
                edge_length = std::min(edge_length, trajectory_dms[t][vertices[i]][vertices[j]]);
            }
            radius = std::max(radius, edge_length);
        }
    }
    return radius;
}



std::vector<Simplex<input_t>> calculate_birth_grades(const DistanceMatrices& trajectory_dms,
                                                     input_t max_metric_value,
                                                     const std::vector<Simplex<input_t>>& simplices) {
    /*

     trajectory_dms {DistanceMatrices} -- per-timestep distance matrices for the trajectories,
        indexed [timestep][i][j].
     max_metric_value {input_t} -- maximum allowed Rips radius; grades whose radius is at or
        above this are dropped.
     simplices {std::vector<Simplex<input_t>>} -- input simplices, each with a 1-element grade
        carrying the running max edge distance from build_VR_subcomplex_spatiotemporal.

     Returns:
     std::vector<Simplex<input_t>> -- new simplices, each carrying a 3-element grade
        (inverted_start_time, end_time, radius). One input simplex may contribute multiple
        output simplices (one per minimal window).
     */
    std::vector<Simplex<input_t>> new_simplices;

    size_t num_time_steps = trajectory_dms.size();

    // radius_map[t0][t1] caches the radius for window [t0, t1] so the
    // poset-minimality check below is a constant-time lookup.
    std::vector<std::vector<input_t>> radius_map(num_time_steps, std::vector<input_t>(num_time_steps, -1));

    size_t previous_length = 0;

    for (size_t simplex_index = 0; simplex_index < simplices.size(); simplex_index++) {
        const std::vector<size_t>& vertices = simplices[simplex_index].vertices;

        for (size_t inverted_start_time = 0; inverted_start_time < num_time_steps; inverted_start_time++) {
            size_t start_time = num_time_steps - 1 - inverted_start_time;

            for (size_t end_time = start_time; end_time < num_time_steps; end_time++) {
                input_t radius = interlevel_rips_radius(trajectory_dms, vertices, start_time, end_time);
                radius_map[start_time][end_time] = radius;

                // Skip if a strictly-smaller window (one step shorter at either
                // end) already births the simplex at radius <= ours — that
                // smaller window's grade dominates ours in the poset.
                bool inserted_already = (start_time < end_time)
                    && ((radius_map[start_time + 1][end_time] <= radius)
                        || (radius_map[start_time][end_time - 1] <= radius));

                if (!inserted_already && radius < max_metric_value) {
                    std::vector<input_t> grade{(input_t)inverted_start_time, (input_t)end_time, radius};
                    new_simplices.push_back(Simplex<input_t>(vertices, grade, simplex_index));
                }
            }
        }

        // Reset the cache for the next simplex.
        for (size_t i = 0; i < num_time_steps; i++) {
            std::fill(radius_map[i].begin(), radius_map[i].end(), -1);
        }

        // Reverse the grades just emitted for this simplex.
        std::reverse(new_simplices.begin() + previous_length, new_simplices.end());
        previous_length = new_simplices.size();
    }
    return new_simplices;
}




std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal_dm(DistanceMatrices& trajectory_dms,
                                            input_t max_metric_value,
                                            int hom_dim){
    /* Computes the two boundary matrices of the dynamic Vietoris-Rips complex needed to compute the
     homology of dimension 'hom_dim', from pre-computed per-timestep distance matrices.

     trajectory_dms {DistanceMatrices} -- per-timestep distance matrices, indexed [timestep][i][j].
        Entries at or above max_metric_value are treated as "edge missing" at that timestep, which
        lets the caller encode unrealised vertex pairs (e.g. two leaves of a lineage tree whose
        ancestors weren't both alive yet).
     max_metric_value {input_t} -- maximum allowed distance for inclusion in the complex.
     hom_dim {int} -- the dimension of the homology to be computed.

     Returns:
     std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>> -- (high boundary matrix, low
        boundary matrix, index_value_lists). The matrices' columns carry integer-index grades;
        index_value_lists is the per-axis lookup table that maps those indices back to the real
        filtration values.
     */
    Matrix high_matrix;
    Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;

    if (trajectory_dms.size() == 0 || trajectory_dms[0].size() == 0) {
        return std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>(high_matrix, low_matrix,
                                                                             index_value_lists);
    }

    std::vector<Simplex<input_t>> low_simplices;
    std::vector<Simplex<input_t>> mid_simplices;
    std::vector<Simplex<input_t>> high_simplices;

    size_t n_agents = trajectory_dms[0].size();
    for (size_t vertex = 0; vertex < n_agents; vertex++) {
        std::vector<size_t> vertices;
        vertices.push_back(vertex);
        build_VR_subcomplex_spatiotemporal(trajectory_dms, vertices, 0, max_metric_value,
                                           low_simplices, mid_simplices, high_simplices, hom_dim);
    }

    std::vector<Simplex<input_t>> high_simplices_extended =
        calculate_birth_grades(trajectory_dms, max_metric_value, high_simplices);
    std::vector<Simplex<input_t>> mid_simplices_extended =
        calculate_birth_grades(trajectory_dms, max_metric_value, mid_simplices);

    std::pair<std::vector<hash_map<input_t, index_t>>, std::vector<std::vector<input_t>>> grade_map =
        compute_grade_map<input_t>(high_simplices_extended, mid_simplices_extended, low_simplices);
    add_boundary_columns_multicritical<input_t>(high_simplices_extended, mid_simplices_extended,
                                                high_matrix, grade_map.first);
    add_boundary_columns<input_t>(mid_simplices_extended, low_simplices, low_matrix, grade_map.first);
    index_value_lists = grade_map.second;

    return std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>(high_matrix, low_matrix,
                                                                         index_value_lists);
}


std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal(Trajectories& trajectories,
                                         Metric* metric,
                                         input_t max_metric_value,
                                         int hom_dim){
    /* Convenience wrapper around compute_boundary_matrices_spatiotemporal_dm: builds per-timestep
     distance matrices from a time-varying point cloud, then forwards.

     trajectories {Trajectories} -- a time-varying point cloud, indexed [point][timestep][spatial_dim].
     metric {Metric*} -- metric used to construct the per-timestep distance matrices.
     max_metric_value {input_t} -- maximum allowed distance for inclusion in the complex.
     hom_dim {int} -- the dimension of the homology to be computed.
     */
    DistanceMatrices trajectory_dms = get_trajectory_dms(trajectories, metric);
    return compute_boundary_matrices_spatiotemporal_dm(trajectory_dms, max_metric_value, hom_dim);
}


std::tuple<Matrix, Matrix, std::vector<std::vector<input_t>>>
compute_boundary_matrices_spatiotemporal_tree(PositionsPerTime& positions_per_t,
                                              std::vector<std::vector<int>>& parents_per_t,
                                              Metric* metric,
                                              input_t max_metric_value,
                                              int hom_dim){
    /* Convenience wrapper around compute_boundary_matrices_spatiotemporal_dm: takes tree-
     structured input, builds per-timestep distance matrices on the leaves (= agents at the final
     timestep), and forwards. Pairs of leaves whose lineages weren't both alive at a given time
     get distance >= max_metric_value at that time, which the engine treats as "edge missing".

     positions_per_t {PositionsPerTime} -- per-timestep positions of agents alive at each timestep,
        indexed [timestep][agent_at_timestep][spatial_dim].
     parents_per_t {std::vector<std::vector<int>>} -- per-timestep parent indices, indexed
        [timestep][agent_at_timestep]. Each entry is an index into positions_per_t[t-1], or -1
        if the agent's lineage starts at this timestep. parents_per_t[0] is unused.
     metric {Metric*} -- metric used to compute pairwise distances between ancestor positions.
     max_metric_value {input_t} -- maximum allowed distance for inclusion in the complex.
     hom_dim {int} -- the dimension of the homology to be computed.
     */
    if (positions_per_t.size() == 0) {
        DistanceMatrices empty_dms;
        return compute_boundary_matrices_spatiotemporal_dm(empty_dms, max_metric_value, hom_dim);
    }
    size_t n_leaves = positions_per_t.back().size();
    std::vector<std::vector<int>> ancestor_lookup = build_ancestor_lookup(parents_per_t, n_leaves);
    DistanceMatrices trajectory_dms =
        tree_positions_to_dms(positions_per_t, ancestor_lookup, metric, max_metric_value);
    return compute_boundary_matrices_spatiotemporal_dm(trajectory_dms, max_metric_value, hom_dim);
}
