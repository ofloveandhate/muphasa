#include "grade.h"

std::vector<std::vector<grade_t>> get_grade_slices(const std::vector<std::vector<index_t>>& base_set)
{
    std::vector<std::vector<grade_t>> grade_slice_list;
    Iterator_lex iterator(base_set);
    while(iterator.has_next()){
        grade_t v = iterator.next();
        index_t sum = 0;
        for(auto& val : v){
            sum += val;
        }
        while(grade_slice_list.size() <= sum){
            grade_slice_list.push_back(std::vector<grade_t>());
        }
        grade_slice_list[sum].push_back(v);
    }
    return grade_slice_list;
}


void helper_grade_slice(grade_t& indices, size_t sum, size_t target_sum, size_t index, std::vector<grade_t*>& grades, std::vector<std::vector<index_t>>& base_set, size_t n_parameters){
    if(index+1==n_parameters){
        if(base_set[index].size()+sum > target_sum){
            std::vector<index_t> grade;
            for(size_t i=0; i<n_parameters-1; i++){
                grade.push_back(base_set[i][indices[i]]);
            }
            grade.push_back(base_set[index][target_sum-sum]);
            grades.push_back(new grade_t(grade));
        }
    }else{
        for(size_t i=0; i<base_set[index].size(); i++){
            if(sum+i <= target_sum){
                indices[index] = i;
                helper_grade_slice(indices, sum+i, target_sum, index+1, grades, base_set, n_parameters);
            }else{
                break;
            }
        }
    }
}

std::vector<grade_t*> get_grade_slice_rec(std::vector<std::vector<index_t>>& base_set, int val){
    std::vector<grade_t*> grade_slice_list;
    grade_t indices;
    for(size_t i=0; i<base_set.size(); i++){
        indices.push_back(0);
    }
    helper_grade_slice(indices, 0, val, 0, grade_slice_list, base_set, base_set.size());
    return grade_slice_list;
}


std::vector<grade_t> get_grade_slice(const std::vector<std::vector<index_t>>& base_set, int val){
    std::vector<grade_t> grade_slice_list;
    Iterator_lex iterator(base_set);
    while(iterator.has_next()){
        iterator.step();
        index_t sum = 0;
        for(auto& val : iterator.indices){
            sum += val;
        }
        if(sum==val){
            grade_slice_list.push_back(grade_t(iterator.values));
        }
    }
    return grade_slice_list;
}

