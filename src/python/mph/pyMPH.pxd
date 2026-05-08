

from libcpp.vector cimport vector
from libcpp.queue cimport priority_queue
from libcpp.pair cimport pair

cdef extern from "bindings.h":
	ctypedef struct GradedMatrix:
		vector[vector[pair[int, int]]] matrix
		vector[vector[int]] row_grades
		vector[vector[int]] column_grades

	ctypedef double input_t

	ctypedef struct PythonCompressedLandscape:
		vector[pair[vector[int], vector[vector[int]] ]] pairings
		vector[vector[input_t]] index_value_lists

	ctypedef struct PythonNaiveBenchFixture:
		vector[pair[vector[int], vector[vector[int]] ]] pairings
		vector[vector[input_t]] index_value_lists
		vector[vector[int]] presentation_matrix
		vector[vector[int]] presentation_column_grades
		vector[vector[int]] presentation_row_grades
		vector[int] v_max

	# ctypedef struct PythonLandscape:
	# 	vector[pair[vector[int], size_t]] landscape

		
	GradedMatrix presentation(vector[vector[input_t]]& _points, vector[int]& _metrics, vector[input_t]& _max_metric_values, vector[int]& _filters, int hom_dim) except +
	
	GradedMatrix presentation_dm(vector[vector[vector[input_t]]]& distance_matrices, vector[input_t]& max_metric_values, vector[vector[input_t]]& filters, int hom_dim) except +

	GradedMatrix presentation_FIrep(vector[vector[int]]& high_matrix, vector[vector[int]]& column_grades_h, vector[vector[int]]& low_matrix, vector[vector[int]]& column_grades_l) except +

	pair[GradedMatrix, GradedMatrix] groebner_bases(vector[vector[int]]& matrix, vector[vector[int]]& row_grades, vector[vector[int]]& column_grades) except +

	PythonCompressedLandscape landscapes_spatiotemporal(vector[vector[vector[input_t]]]& trajectories, input_t max_metric_value, int hom_dim) except +

	PythonCompressedLandscape landscapes_spatiotemporal_tree(vector[vector[vector[input_t]]]& positions_per_t, vector[vector[int]]& parents_per_t, input_t max_metric_value, int hom_dim) except +

	vector[vector[input_t]] eval_sparse_landscape_batch(vector[pair[vector[input_t], vector[vector[input_t]]]]& pairings, vector[vector[input_t]]& grades, int k_max) except +

	vector[int] eval_sparse_landscape_batch_int(vector[pair[vector[int], vector[vector[int]]]]& pairings, vector[vector[int]]& grades, int k) except +

	PythonNaiveBenchFixture naive_bench_setup_spatiotemporal(vector[vector[vector[input_t]]]& trajectories, input_t max_metric_value, int hom_dim) except +

	vector[int] eval_naive_rank_int(vector[vector[int]]& presentation_matrix, vector[vector[int]]& presentation_column_grades, vector[vector[int]]& presentation_row_grades, vector[vector[int]]& query_grades, int landscape_dim) except +

	vector[int] eval_naive_cached_int(vector[vector[int]]& presentation_matrix, vector[vector[int]]& presentation_column_grades, vector[vector[int]]& presentation_row_grades, vector[vector[int]]& query_grades, int landscape_dim) except +

	# PythonLandscape landscapes_spatiotemporal_naive(vector[vector[vector[input_t]]] trajectories, input_t max_metric_value, int hom_dim, int landscape_dim) except +
