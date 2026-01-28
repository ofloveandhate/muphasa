import math, functools
from typing import List, Tuple

import numpy as np

from .pyMPH import doPresentation, doPresentationFIrep, doPresentationDm, computeGroebnerBases, computeSpatiotemporalLandscapesSparse


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
    """
    Grade in `mph`.  A Grade can represent either a generator or a syzygy.  It consists of, essentially, a list of floating point values.

    This is derived from `tuple`, so it is an integer-indexed container starting at 0.  It provides comparitor methods <, >, <=, >=.
    """

    @classmethod
    def _comparison_error(cls, v, w):
        return f"Cannot compare grade of size {len(v)} with grade of size {len(w)}"

    def __getitem__(self, i):
        item = super().__getitem__(i)
        return Grade(item) if isinstance(i, slice) else item

    def __lt__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(Grade._comparison_error(self, v))
        found_difference = False
        for x, y in zip(self, v):
            if x > y:
                return False
            if x < y:
                found_difference = True
        return found_difference

    def __le__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(Grade._comparison_error(self, v))
        return all(self.__getitem__(i) <= v[i] for i in range(len(self)))
   
    def __gt__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(Grade._comparison_error(self, v))
        found_difference = False
        for x, y in zip(self, v):
            if x < y:
                return False
            if x > y:
                found_difference = True
        return found_difference

    
    def __ge__(self, v) -> bool:
        if len(self) != len(v):
            raise ValueError(Grade._comparison_error(self, v))
        return all(self.__getitem__(i) >= v[i] for i in range(len(self)))


    @classmethod
    def colex_compare(cls, v, w):
        if len(v) != len(w):
            raise ValueError(Grade._comparison_error(v, w))
        for vi, wi in (zip(reversed(v), reversed(w))):
            if vi < wi:
                return -1
            if wi > vi:
                return 1
        return 0


    @classmethod
    def join(cls, v, w):
        """
        Join two grades, by zipping them (as tuples) and applying the `max` function to each item.  This is essentially the component-wise max of two `Grade` s.
        """
        return Grade(map(max, zip(v, w)))
    

class SparseLandscape():
    """
    SparseLandscape in `mph`
    """

    def __init__(self, pairings: List[Tuple[Grade, List[Grade]]]):
        self.pairings = pairings
        

    def eval(self, grade, k):
        """
        Evaluate the landscape layer k at the given grade
        """

        grade = Grade(grade)
        if k <= 0:
            return 0
        diffs = []
        for birth, syzygies in filter(lambda pairing : pairing[0] <= grade, self.pairings):
            low_distance = grade[-1] - birth[-1]
            high_distance = math.inf
            for syzygy in sorted(syzygies, key=functools.cmp_to_key(Grade.colex_compare)):
                if syzygy[:-1] <= grade[:-1]:
                    high_distance = syzygy[-1] - grade[-1]
                    break
            diffs.append(max(min(low_distance, high_distance), 0))
        return sorted(diffs, reverse=True)[k-1] if k <= len(diffs) else 0

    def save(self, fname: str, precision=3):
        """
        Save the `SparseLandscape` to disk, under filename `fname`, with number-of-digits precision `precision`
        """
        
        def grade_to_str(grade):
            return ', '.join(f'{x:.{precision}f}' for x in grade)
        
        with open(fname, 'w', encoding="utf-8") as f:
            for v, syzygies in self.pairings:
                f.write('F ' + grade_to_str(v) + '\n')
                for syzygy in syzygies:
                    f.write('S ' + grade_to_str(syzygy) + '\n')
                    
    @classmethod
    def load(cls, fname: str):
        """
        Load a `SparseLandscape` from file `fname`.

        It reads a plaintext file, with lines.  Each line starts with a tag, either `F` or `S`.  
        Behind the tag is a ", "-separated list of floats, which represent a `Grade`.  

        If the tag is `F`, then the line represents a generator.
        If the tag is `S`, then the line represents a syzygy.

        The structure is an F, and then a bunch of S's.  This repeats.

        If an S appears before a generator is set via an F line, a exception is raised.  
        If a tag is invalid, an exception is raised.
        """

        pairings = []
        
        with open(fname, 'r', encoding='utf-8') as f:

            syzygies = []
            generator = None

            for lineno, line in enumerate(f.readlines()):
                tag = line[0]
                grade = list(map(float, line[1:].rstrip().lstrip().split(', ')))
                
                # New generator
                if tag == 'F':

                    # store, because found a new one.  (the only time this won't happen is the first generator)
                    if generator is not None:
                        pairings.append((generator, syzygies))

                    # set the current generator
                    generator = grade

                    # reset to empty for this generator.  subsequent lines will read the syz's
                    syzygies = []
                
                # New syzygy
                elif tag == 'S': 
                    if generator is None:
                        raise ValueError(f"Must provide a generator before a syzygy.  See line {lineno} in file {fname}")

                    syzygies.append(grade)

                else:
                    raise ValueError(f"The first character of a line in a SparseLandscape file must be a valid tag (F or S).  File {fname} has {tag} at line {lineno}")

            # wrap up by taking care of the last generator and its syzygyies
            if generator is not None:
                pairings.append((generator, syzygies)) 
        
        # finally, pass to the class constructor to make a SparseLandscape
        return cls(pairings)
                

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
    
    def transform_grade(grade: Grade) -> Grade:
        return Grade(index_value_lists[parameter][v] for parameter, v in enumerate(grade))
    
    pairings = [ (transform_grade(v), list(map(transform_grade, syzygies)))
                                    for v, syzygies in ret["pairings"]]

    return SparseLandscape(pairings)

# def compute_spatiotemporal_landscapes_naive(trajectories: np.ndarray, max_metric_value: float, hom_dim: int, landscape_dim: int):
#     ret = computeSpatiotemporalLandscapesNaive(np.ascontiguousarray(trajectories, dtype=np.float64), max_metric_value*max_metric_value, hom_dim, landscape_dim)
#     return ret

__all__ = ["presentation", "presentation_dm", "presentation_FIrep", "groebner_bases", "Grade", "GradedMatrix", "compute_spatiotemporal_landscapes_sparse", "SparseLandscape"] + ['List',
'Tuple']
