import math, functools
from typing import List, Tuple

import numpy as np

from pyMPH import doPresentation, doPresentationFIrep, doPresentationDm, computeGroebnerBases, computeSpatiotemporalLandscapesSparse


METRIC_DICT = {
    "euclidean" : 0,
    "codensity" : 1
}
FILTER_DICT = {
    "codensity" : -1
}

class GradedMatrix():
    def __init__(self, matrix: List[List[Tuple[int, int]]], column_grades: List[List[int]], row_grades: List[List[int]]):
        self.matrix = matrix
        self.column_grades = column_grades
        self.row_grades = row_grades

    def asnparray(self):
        _matrix = np.zeros(shape=(len(self.row_grades), len(self.column_grades)), dtype=np.int32)
        for index, entry in enumerate(self.matrix):
            _matrix[entry[1], index] =  entry[0]
        return _matrix

    def __len__(self):
        return len(self.column_grades)

    def __str__(self):
        matrix = self.asnparray()
        s = [[" "]+[str(grade) for grade in self.column_grades]]+[[str(e) for e in [self.row_grades[i]]+matrix[i, :].tolist()] for i in range(matrix.shape[0])]
        lens = [max(map(len, col)) for col in zip(*s)]
        fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
        table = [fmt.format(*row) for row in s]
        return '\n'.join(table)

def parse_log_level(log_level: str) -> int:
    if log_level == "info":
        return 1
    if log_level == "debug":
        return 2
    return 0


class Grade(tuple):
    
    def __getitem__(self, i):
        item = super().__getitem__(self, i)
        return Grade(item) if isinstance(i, slice) else item

    def __lt__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(f"Cannot compare grade of size {len(self)} with grade of size {len(v)}")
        return all(self.__getitem__(i) < v[i] for i in range(len(self)))

    def __le__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(f"Cannot compare grade of size {len(self)} with grade of size {len(v)}")
        return all(self.__getitem__(i) <= v[i] for i in range(len(self)))
   
    def __gt__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(f"Cannot compare grade of size {len(self)} with grade of size {len(v)}")
        return all(self.__getitem__(i) > v[i] for i in range(len(self)))
    
    def __ge__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(f"Cannot compare grade of size {len(self)} with grade of size {len(v)}")
        return all(self.__getitem__(i) >= v[i] for i in range(len(self)))


    @classmethod
    def colex_compare(cls, v, w):
        if len(v) != len(w):
            raise ValueError(f"Cannot compare grade of size {len(v)} with grade of size {len(w)}")
        for vi, wi in (zip(reversed(v), reversed(w))):
            if vi < wi:
                return -1
            if wi > vi:
                return 1
        return 0


    @classmethod
    def join(cls, v, w):
        return Grade(map(max, zip(v, w)))
    

class SparseLandscape():
    def __init__(self,
                 pairings: List[Tuple[Grade, List[Grade]]],
                 index_value_lists: List[List[np.float64]]
                ):
        
        def transform_grade(grade: Grade) -> Grade:
            return Grade(index_value_lists[parameter][v] for parameter, v in enumerate(grade))

        self.pairings = [ (transform_grade(v), list(map(transform_grade, syzygies)))
                                    for v, syzygies in pairings]
        
        self.index_value_lists = index_value_lists


    def eval(self, grade, k):
        """Evaluate the landscape layer k at the given grade"""
        if k <= 0:
            return 0
        diffs = []  # D
        for birth, syzygies in filter(lambda birth, _ : birth < grade, self.pairings):
            low_distance = grade[-1] - birth[-1]
            high_distance = math.inf
            for syzygy in sorted(syzygies, key=functools.cmp_to_key(Grade.colex_compare)):
                if syzygy[:-1] <= grade[:-1]:
                    high_distance = syzygy[-1] - grade[-1]
                    break
            diffs.append(max(min(low_distance, high_distance), 0))

        return sorted(diffs, reverse=True)[k-1] if k < len(diffs) else 0

def translate_metrics_filters(metrics: List[str], filters: List[str]):
    for metric in metrics: 
        assert metric in METRIC_DICT
    for _filter in filters:
        assert _filter in FILTER_DICT or _filter.isdigit()
    metric_trans = [METRIC_DICT[metric] for metric in metrics]
    filter_trans = [int(_filter) if _filter.isdigit() else FILTER_DICT[_filter] for _filter in filters]
    return np.ascontiguousarray(metric_trans, dtype=np.int32), np.ascontiguousarray(filter_trans, dtype=np.int32)

def presentation(points: np.ndarray, metrics: List[str]=["euclidean"], max_metric_values: List[float]=[100], filters: List[str]=["1"], hom_dim: int=1, log_level:str="silent") -> GradedMatrix:
    _metrics, _filters = translate_metrics_filters(metrics, filters)
    ret = doPresentation(np.ascontiguousarray(points, dtype=np.float64), _metrics, np.ascontiguousarray(max_metric_values, dtype=np.float64), _filters, hom_dim)
    return GradedMatrix(ret["matrix"], ret["column_grades"], ret["row_grades"])

def presentation_dm(distance_matrices: np.ndarray, max_metric_values: List[float]=None, filters: np.ndarray=np.asarray([[]]), hom_dim:int=1, log_level:str="silent") -> GradedMatrix:
    if not max_metric_values:
        max_metric_values = [ np.max(distance_matrices[i, :, :]) for i in range(distance_matrices.shape[0])]
    ret = doPresentationDm(np.ascontiguousarray(distance_matrices, dtype=np.float64), np.ascontiguousarray(max_metric_values, dtype=np.float64), np.ascontiguousarray(filters, dtype=np.float64), hom_dim)
    return GradedMatrix(ret["matrix"], ret["column_grades"], ret["row_grades"])

def presentation_FIrep(high_matrix: np.ndarray, column_grades_h: List[List[int]], low_matrix: np.ndarray, column_grades_l: List[List[int]], log_level:str="silent") -> GradedMatrix:
    ret = doPresentationFIrep(np.ascontiguousarray(high_matrix.T, dtype=np.int32), np.ascontiguousarray(column_grades_h, dtype=np.int32),
        np.ascontiguousarray(low_matrix.T, dtype=np.int32), np.ascontiguousarray(column_grades_l, dtype=np.int32))
    return GradedMatrix(ret["matrix"], ret["column_grades"], ret["row_grades"])

def groebner_bases(matrix: np.ndarray, column_grades: List[List[int]], row_grades: List[List[int]]=None, log_level:str="silent") -> Tuple[GradedMatrix, GradedMatrix]:
    if not row_grades:
        row_grades = [ [0]*len(column_grades[0]) ] * matrix.shape[0]
    # Transpose input matrix to simplify matrix construction in cpp
    ret = computeGroebnerBases(np.ascontiguousarray(matrix.T, dtype=np.int32), np.ascontiguousarray(row_grades, dtype=np.int32), np.ascontiguousarray(column_grades, dtype=np.int32))
    return GradedMatrix(ret[0]["matrix"], ret[0]["column_grades"], ret[0]["row_grades"]), GradedMatrix(ret[1]["matrix"], ret[1]["column_grades"], ret[1]["row_grades"])

def compute_spatiotemporal_landscapes_sparse(trajectories: np.ndarray, max_metric_value: float, hom_dim: int) -> SparseLandscape:
    ret = computeSpatiotemporalLandscapesSparse(np.ascontiguousarray(trajectories, dtype=np.float64), max_metric_value*max_metric_value, hom_dim)
    index_value_lists = ret["index_value_lists"]
    index_value_lists[2] = list(map(math.sqrt, index_value_lists[2]))
    return SparseLandscape(ret["pairings"], index_value_lists)

# def compute_spatiotemporal_landscapes_naive(trajectories: np.ndarray, max_metric_value: float, hom_dim: int, landscape_dim: int):
#     ret = computeSpatiotemporalLandscapesNaive(np.ascontiguousarray(trajectories, dtype=np.float64), max_metric_value*max_metric_value, hom_dim, landscape_dim)
#     return ret

__all__ = ["presentation", "presentation_dm", "presentation_FIrep", "groebner_bases", "GradedMatrix", "compute_spatiotemporal_landscapes_sparse"]
