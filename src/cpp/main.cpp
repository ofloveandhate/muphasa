//
//  main.cpp
//  mph
//
//  Created by Oliver on 2020-04-09.
//  Copyright © 2020 Oliver. All rights reserved.
//

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <queue>
#include <iostream>
#include <unordered_map>
#include <set>
#include <time.h>


#include "main.h"
#include "presentation.h"



double res_memory=0;
double virt_memory=0;


CompressedLandscape computeCompressedLandscape(std::vector<SignatureColumn>& presentation, std::vector<grade_t> row_grades){
    /*
     The main function computing a compressed persistence landscape from a minimal presentation described by the list of columns 'presentation' and the row grades 'row_grades'.

     Arguments:
     presentation {std::vector<SignatureColumn>} -- columns describing a minimal presentation.
     row_grades {std::vector<grade_t>} grades of the rows in the matrix 'presentation'.

     Returns:
     CompressedLandscape -- a compressed representation of a multiparameter persistence landscape.
     */

    std::cout << "Starting to compute persistence landscape..." << std::endl;

    /* Sort columns colexicographically */
    sort(presentation.begin(), presentation.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<presentation.size(); i++){
        presentation[i].signature_index = i;
        index_map_high[presentation[i].grade[presentation[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;

    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    gb_columns.reserve(2*presentation.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(presentation.size());
    Syz.reserve(presentation.size());

    for(size_t i=0; i<presentation.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : presentation){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : presentation){
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(presentation.size());
    GradeHasher grade_hasher;
    for(auto& c : presentation){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    index_t pivot;

    int iter_index = 0;

    CompressedLandscape landscape;
    for(size_t i=0; i<row_grades.size(); i++){
        landscape.push_back(std::pair<grade_t, std::vector<grade_t>>(row_grades[i], std::vector<grade_t>()));
    }

    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }

        size_t grade_hash = grade_hasher(v);
        iter_index++;

        /* Initialize Macaulay matrix */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                bool in_syz = false;
                if(Syz[i].size()>0){
                    for(size_t j=Syz[i].size()-1; j < Syz[i].size(); j--){
                        if((Syz[i][j].get_grade()).leq_poset(v)){
                            in_syz = true;
                            break;
                        }
                    }
                }
                if(!in_syz){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(grade_hash == grade_hashes[i] && presentation[i].grade == v){
                    gb_columns.push_back(presentation[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        bool add_to_landscape = true;
                        for ( auto& w : landscape[pivot].second ) {
                            if ( w.leq_poset(v) ) {
                                add_to_landscape = false;
                                break;
                            }
                        }
                        if ( add_to_landscape ) {
                            landscape[pivot].second.push_back(v);
                        }
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_lists[pivot].size() == 0 || working_column.grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(working_column.grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grade_t g = working_column.grade.join(grade_lists[pivot][j]);
                                        if(visited_grades.find(g) == visited_grades.end()){
                                            grades.push(g);
                                            visited_grades.insert(g);
                                        }
                                    }
                                }
                                grade_lists[pivot].push_back(working_column.grade);
                            }
                        }
                    }else{
                        Syz[i].push_back(signature_t(working_column.grade, -1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        bool add_to_landscape = true;
                        for ( auto& w : landscape[pivot].second ) {
                            if ( w.leq_poset(v) ) {
                                add_to_landscape = false;
                                break;
                            }
                        }
                        if ( add_to_landscape ) {
                            landscape[pivot].second.push_back(v);
                        }
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(gb_columns[column_index].grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][j]);
                                        if(visited_grades.find(g) == visited_grades.end()){
                                            grades.push(g);
                                            visited_grades.insert(g);
                                        }
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "Finished computing compressed landscape." << std::endl;
    landscape.sort_colex();
    return landscape;
}

DiagLandscape computeSparseDiagLandscape(Matrix& presentation, std::vector<grade_t> row_grades){
    DiagLandscape landscape;
    for(size_t i=0; i<row_grades.size(); i++){
        landscape.push_back(std::pair<grade_t, std::vector<std::vector<grade_t>>>(row_grades[i], std::vector<std::vector<grade_t>>()));
    }
    size_t n_grade = row_grades[0].size();
    for(size_t i=0; i<n_grade; i++){
        Matrix presentation_copy(presentation);
        for(auto& column : presentation_copy){
            index_t val = column.grade[n_grade-1];
            column.grade[n_grade-1] = column.grade[i];
            column.grade[i] = val;
        }


        // Also switch index of rows and reorder, then keep track of reordering so it can be mapped back to the original index
        std::vector<std::pair<size_t, grade_t>> zipped_row_grades;
        size_t index = 0;
        for(auto& grade : row_grades){
            grade_t new_grade(grade);
            index_t val = new_grade[n_grade-1];
            new_grade[n_grade-1] = new_grade[i];
            new_grade[i] = val;
            zipped_row_grades.push_back(std::pair<size_t, grade_t>(index, new_grade));
            index++;
        }

        sort(zipped_row_grades.begin(), zipped_row_grades.end(), [ ](  std::pair<size_t, grade_t>& lhs,  std::pair<size_t, grade_t>& rhs )
             {
                 return lhs.second.lt_colex(rhs.second);
             });

        std::vector<grade_t> new_row_grades;
        hash_map<size_t, size_t> row_reorder_map;
        for(size_t j=0; j<zipped_row_grades.size(); j++){
            new_row_grades.push_back(zipped_row_grades[j].second);
            row_reorder_map[zipped_row_grades[j].first] = j;
        }

        Matrix reindex_presentation;
        for(size_t i=0; i<presentation_copy.size(); i++) {
            SignatureColumn pres_copy(presentation_copy[i]);
            SignatureColumn working_column(presentation_copy[i].grade, i);
            while(!pres_copy.empty()){
                column_entry_t entry = pres_copy.pop_pivot();
                working_column.push(column_entry_t(1, row_reorder_map[entry.get_index()]));
            }
            reindex_presentation.push_back(working_column);
        }

        CompressedLandscape c_landscape = computeCompressedLandscape(presentation_copy, new_row_grades);
        for(size_t j=0; j<row_grades.size(); j++){
            landscape[zipped_row_grades[j].first].second.push_back(c_landscape[j].second);
        }
    }
    return landscape;
}


Landscape compute_landscape_naive_rank(Matrix presentation, std::vector<grade_t> row_grades, grade_t v_max, int lanscape_dim) {
    Landscape lanscape;

    std::vector<std::vector<index_t>> base_set;
    for(size_t i=0; i<v_max.size(); i++){
        std::vector<index_t> row;
        for(size_t j=0; j<v_max[i]; j++){
            row.push_back(j);
        }
        base_set.push_back(row);
    }

    index_t max_pivot=0;
    for(auto& column : presentation){
        if(max_pivot < column.get_pivot().get_index()){
            max_pivot = column.get_pivot().get_index();
        }
    }

    Iterator_lex grade_iterator = Iterator_lex(base_set);

    grade_t vv;
    vv.push_back(4);vv.push_back(4);vv.push_back(7);

    while(grade_iterator.has_next()){
        grade_t v = grade_iterator.next();
        size_t max_rank_index = 0;
        if ( v == vv ) {
            std::cout << std::endl;
        }
        for ( size_t k=0; k<=v[v.size()-1]; k++ ) {
            // Calculate rank of map from v-i to v+i
            index_t pivot;
            std::vector<index_t> pivot_map(max_pivot+1, -1);
            grade_t w(v);
            w[w.size()-1] += k;
            Matrix reduced_basis;
            for( size_t i=0; i<presentation.size(); i++ ){
                if ( presentation[i].grade.leq_poset(w) ) {
                    SignatureColumn working_column = presentation[i];
                    pivot = working_column.get_pivot().get_index();
                    while(pivot != -1 && pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(reduced_basis[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        pivot_map[pivot] = reduced_basis.size();
                        reduced_basis.push_back(working_column);
                    }
                }
            }
            w[w.size()-1] = v[v.size()-1]-k;
            int current_rank = 0;
            for(size_t i=0; i<row_grades.size(); i++) {
                if ( row_grades[i].leq_poset(w) && pivot_map[i] == -1 ) {
                    current_rank++;
                }
            }
            if (current_rank < lanscape_dim) {
                break;
            } else {
                max_rank_index = k;
            }
        }
        lanscape.push_back(std::pair<grade_t, size_t>(v, max_rank_index));
    }
    return lanscape;
}

Landscape computeLandscapeNaive(Matrix presentation, std::vector<grade_t> row_grades, grade_t v_max, int lanscape_dim) {
    Landscape lanscape;

    std::vector<std::vector<index_t>> base_set;
    for(size_t i=0; i<v_max.size(); i++){
        std::vector<index_t> row;
        for(size_t j=0; j<v_max[i]; j++){
            row.push_back(j);
        }
        base_set.push_back(row);
    }

    index_t max_pivot=0;
    for(auto& column : presentation){
        if(max_pivot < column.get_pivot().get_index()){
            max_pivot = column.get_pivot().get_index();
        }
    }

    Iterator_lex grade_iterator = Iterator_lex(base_set);

    std::unordered_map<grade_t, std::pair<Matrix, Matrix>, GradeHasher> H_basis_map;

    while(grade_iterator.has_next()){
        grade_t v = grade_iterator.next();
        size_t max_rank_index = 0;
        for ( size_t k=0; k<=v[v.size()-1]; k++ ) {
            // Calculate rank of map from v-i to v+i
            grade_t w(v);
            w[w.size()-1] += k;

            if ( H_basis_map.find(w) == H_basis_map.end() ) {
                index_t pivot;
                std::vector<index_t> pivot_map(max_pivot+1, -1);

                Matrix reduced_basis;
                for( size_t i=0; i<presentation.size(); i++ ){
                    if ( presentation[i].grade.leq_poset(w) ) {
                        SignatureColumn working_column = presentation[i];
                        pivot = working_column.get_pivot().get_index();
                        while(pivot != -1 && pivot_map[pivot] > -1){
                            working_column.plus(reduced_basis[pivot_map[pivot]]);
                            pivot = working_column.get_pivot().get_index();
                        }
                        if(pivot != -1){
                            pivot_map[pivot] = reduced_basis.size();
                            reduced_basis.push_back(working_column);
                        }
                    }
                }
                size_t im_basis_index = reduced_basis.size();
                // Construct H-basis
                for ( size_t i=0; i<row_grades.size(); i++ ) {
                    if (row_grades[i].leq_poset(w)){
                        SignatureColumn working_column(row_grades[i], reduced_basis.size());
                        working_column.push(column_entry_t(1, i));
                        pivot = working_column.get_pivot().get_index();
                        while(pivot != -1 && pivot_map[pivot] > -1){
                            working_column.plus(reduced_basis[pivot_map[pivot]]);
                            pivot = working_column.get_pivot().get_index();
                        }
                        if(pivot != -1){
                            pivot_map[pivot] = reduced_basis.size();
                            reduced_basis.push_back(working_column);
                        }
                    }
                }
                hash_map<size_t, size_t> reindex_map;
                size_t pivot_index = 0;
                for(size_t i=0; i<row_grades.size(); i++){
                    if( row_grades[i].leq_poset(w) ) {
                        reindex_map[i] = pivot_index;
                        pivot_index++;
                    }
                }

                if ( pivot_index != reduced_basis.size() ) {
                    throw std::runtime_error("not basis");
                }

                // Reindex reduced_basis
                Matrix H_basis;
                Matrix reindexed_reduced_basis;
                for(size_t i=0; i<reduced_basis.size(); i++) {
                    SignatureColumn sig_copy(reduced_basis[i]);
                    SignatureColumn working_column(reduced_basis[i].grade, reduced_basis[i].signature_index);
                    while(!sig_copy.empty()){
                        column_entry_t entry = sig_copy.pop_pivot();
                        working_column.push(column_entry_t(1, reindex_map[entry.get_index()]));
                    }
                    reindexed_reduced_basis.push_back(working_column);
                    if ( i>=im_basis_index ) {
                        H_basis.push_back(working_column);
                    }
                }

                Matrix reordered_basis;
                for(size_t i=reindexed_reduced_basis.size()-H_basis.size(); i<reindexed_reduced_basis.size(); i++) {
                    reordered_basis.push_back(reindexed_reduced_basis[i]);
                }
                for(size_t i=0; i<reindexed_reduced_basis.size()-H_basis.size(); i++) {
                    reordered_basis.push_back(reindexed_reduced_basis[i]);
                }
                Matrix reduced_basis_inverse = inverse(reordered_basis);

                Matrix H_basis_inv = mat_v_cut(reduced_basis_inverse, H_basis.size());

                H_basis_map[w] = std::pair<Matrix, Matrix>(H_basis, H_basis_inv);
            }
            grade_t z(v);
            z[z.size()-1] = v[v.size()-1]-k;
            if ( H_basis_map.find(z) == H_basis_map.end() ) {
                break;
            } else {
                Matrix H_basis = H_basis_map[z].first;
                Matrix H_basis_inv = H_basis_map[w].second;
                hash_map<size_t, size_t> basis_order_map;
                for (size_t i=0; i<H_basis_map[w].first.size(); i++) {
                    basis_order_map[H_basis_map[w].first[i].signature_index] = i;
                }
                Matrix A;
                for (auto& column : H_basis) {
                    SignatureColumn working_column(column.grade, A.size());
                    working_column.push(column_entry_t(1, basis_order_map[column.signature_index]));
                    A.push_back(working_column);
                }


                Matrix inter = matmul(A, H_basis);
                Matrix B = matmul(H_basis_inv, inter);
                std::vector<index_t> pivot_map(row_grades.size()+1, -1);

                for( size_t i=0; i<B.size(); i++ ){
                    index_t pivot = B[i].get_pivot().get_index();
                    while(pivot != -1 && pivot_map[pivot] > -1){
                        B[i].plus(B[pivot_map[pivot]]);
                        pivot = B[i].get_pivot().get_index();
                    }
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }

                size_t rank = 0;
                for(size_t i=0; i<row_grades.size(); i++) {
                    if ( pivot_map[i] != -1 ) {
                        rank++;
                    }
                }
                if (rank < lanscape_dim) {
                    break;
                } else {
                    max_rank_index = k;
                }
            }
        }
        lanscape.push_back(std::pair<grade_t, size_t>(v, max_rank_index));
    }
    return lanscape;
}

Landscape computeLandscapeNaive_diag(Matrix presentation, std::vector<grade_t> row_grades, grade_t v_max, int lanscape_dim) {
    Landscape lanscape;

    std::vector<std::vector<index_t>> base_set;
    for(size_t i=0; i<v_max.size(); i++){
        std::vector<index_t> row;
        for(size_t j=0; j<v_max[i]; j++){
            row.push_back(j);
        }
        base_set.push_back(row);
    }

    index_t max_pivot=0;
    for(auto& column : presentation){
        if(max_pivot < column.get_pivot().get_index()){
            max_pivot = column.get_pivot().get_index();
        }
    }

    Iterator_lex grade_iterator = Iterator_lex(base_set);

    std::unordered_map<grade_t, std::pair<Matrix, Matrix>, GradeHasher> H_basis_map;

    while(grade_iterator.has_next()){
        grade_t v = grade_iterator.next();
        size_t vr_max = 0;
        for (size_t i=0; i<v.size(); i++){
            vr_max = v[i] > vr_max ? v[i] : vr_max;
        }
        size_t max_rank_index = 0;
        for ( size_t k=0; k<=vr_max; k++ ) {
            // Calculate rank of map from v-i to v+i
            grade_t w(v);
            for (size_t i=0; i<w.size(); i++) {
                w[i] += k;
            }
            if ( H_basis_map.find(w) == H_basis_map.end() ) {
                index_t pivot;
                std::vector<index_t> pivot_map(max_pivot+1, -1);

                Matrix reduced_basis;
                for( size_t i=0; i<presentation.size(); i++ ){
                    if ( presentation[i].grade.leq_poset(w) ) {
                        SignatureColumn working_column = presentation[i];
                        pivot = working_column.get_pivot().get_index();
                        while(pivot != -1 && pivot_map[pivot] > -1){
                            working_column.plus(reduced_basis[pivot_map[pivot]]);
                            pivot = working_column.get_pivot().get_index();
                        }
                        if(pivot != -1){
                            pivot_map[pivot] = reduced_basis.size();
                            reduced_basis.push_back(working_column);
                        }
                    }
                }
                size_t im_basis_index = reduced_basis.size();
                // Construct H-basis
                for ( size_t i=0; i<row_grades.size(); i++ ) {
                    if (row_grades[i].leq_poset(w)){
                        SignatureColumn working_column(row_grades[i], i);
                        working_column.push(column_entry_t(1, i));
                        pivot = working_column.get_pivot().get_index();
                        while(pivot != -1 && pivot_map[pivot] > -1){
                            working_column.plus(reduced_basis[pivot_map[pivot]]);
                            pivot = working_column.get_pivot().get_index();
                        }
                        if(pivot != -1){
                            pivot_map[pivot] = reduced_basis.size();
                            reduced_basis.push_back(working_column);
                        }
                    }
                }
                hash_map<size_t, size_t> reindex_map;
                size_t pivot_index = 0;
                for(size_t i=0; i<row_grades.size(); i++){
                    if( row_grades[i].leq_poset(w) ) {
                        reindex_map[i] = pivot_index;
                        pivot_index++;
                    }
                }

                if ( pivot_index != reduced_basis.size() ) {
                    throw std::runtime_error("not basis");
                }

                // Reindex reduced_basis
                Matrix H_basis;
                Matrix reindexed_reduced_basis;
                for(size_t i=0; i<reduced_basis.size(); i++) {
                    SignatureColumn sig_copy(reduced_basis[i]);
                    SignatureColumn working_column(reduced_basis[i].grade, reduced_basis[i].signature_index);
                    while(!sig_copy.empty()){
                        column_entry_t entry = sig_copy.pop_pivot();
                        working_column.push(column_entry_t(1, reindex_map[entry.get_index()]));
                    }
                    reindexed_reduced_basis.push_back(working_column);
                    if ( i>=im_basis_index ) {
                        H_basis.push_back(working_column);
                    }
                }

                Matrix reordered_basis;
                for(size_t i=reindexed_reduced_basis.size()-H_basis.size(); i<reindexed_reduced_basis.size(); i++) {
                    reordered_basis.push_back(reindexed_reduced_basis[i]);
                }
                for(size_t i=0; i<reindexed_reduced_basis.size()-H_basis.size(); i++) {
                    reordered_basis.push_back(reindexed_reduced_basis[i]);
                }
                Matrix reduced_basis_inverse = inverse(reordered_basis);

                Matrix I = matmul(reordered_basis, reduced_basis_inverse);

                for ( size_t i=0; i<I.size(); i++ ) {
                    while(!I[i].empty()) {
                        column_entry_t p = I[i].pop_pivot();
                        if ( p.get_index() != i ) {
                            throw std::runtime_error("Non-identity matrix");
                        }
                    }
                }

                Matrix H_basis_inv = mat_v_cut(reduced_basis_inverse, H_basis.size());

                H_basis_map[w] = std::pair<Matrix, Matrix>(H_basis, H_basis_inv);
            }


            grade_t z(v);
            for (size_t i=0; i<z.size(); i++) {
                z[i] = v[i]-k;
            }
            if ( H_basis_map.find(z) == H_basis_map.end() ) {
                break;
            } else {
                Matrix H_basis = H_basis_map[z].first;
                Matrix H_basis_inv = H_basis_map[w].second;
                hash_map<size_t, size_t> basis_order_map;
                for (size_t i=0; i<H_basis_map[w].first.size(); i++) {
                    basis_order_map[H_basis_map[w].first[i].signature_index] = i;
                }
                Matrix A;

                for (auto& column : H_basis) {
                    SignatureColumn working_column(column.grade, A.size());
                    working_column.push(column_entry_t(1, basis_order_map[column.signature_index]));
                    A.push_back(working_column);
                }

                Matrix inter = matmul(A, H_basis);
                Matrix B = matmul(H_basis_inv, inter);
                std::vector<index_t> pivot_map(row_grades.size()+1, -1);

                for( size_t i=0; i<B.size(); i++ ){
                    index_t pivot = B[i].get_pivot().get_index();
                    while(pivot != -1 && pivot_map[pivot] > -1){
                        B[i].plus(B[pivot_map[pivot]]);
                        pivot = B[i].get_pivot().get_index();
                    }
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }

                size_t rank = 0;
                for(size_t i=0; i<row_grades.size(); i++) {
                    if ( pivot_map[i] != -1 ) {
                        rank++;
                    }
                }
                if (rank < lanscape_dim) {
                    break;
                } else {
                    max_rank_index = k;
                }
            }
        }
        lanscape.push_back(std::pair<grade_t, size_t>(v, max_rank_index));
    }
    return lanscape;
}


Landscape compute_landscape_peaks(MultigradedBasis ker_basis, MultigradedBasis p_basis) {
    Landscape lanscape;


    return lanscape;
}

Landscape compute_landscape_naive(MultigradedBasis ker_basis, MultigradedBasis p_basis, grade_t v_max) {
    Landscape lanscape;

    size_t v_r = v_max[v_max.size()-1];

    std::vector<std::vector<index_t>> base_set;
    for(size_t i=0; i<v_max.size()-1; i++){
        std::vector<index_t> row;
        for(size_t j=0; j<v_max[i]; j++){
            row.push_back(j);
        }
        base_set.push_back(row);
    }

    Iterator_lex grade_iterator = Iterator_lex(base_set);

    while(grade_iterator.has_next()){
        grade_t v = grade_iterator.next();
        for ( size_t i=0; i<v_r; i++ ) {

        }
    }

    return lanscape;
}



Metric* parse_metric(int metric_index){
    switch(metric_index){
        case 0:
            return new SquaredEuclideanMetric();
        default:
            throw std::runtime_error("Failed to parse metric");
    }
}

Filter* parse_filter(int filter_index){
    if(filter_index >= 0){
        return new XFilter(filter_index);
    }
    switch(filter_index){
        case -1:
            return new XFilter(0); //TODO: implement more filters.
        default:
            throw std::runtime_error("Failed to parse metric");
    }
}

Matrix translateInputMatrix(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& column_grades){
    Matrix M;
    for(size_t i=0; i<matrix.size(); i++){
        grade_t grade;
        for(size_t j=0; j<column_grades[i].size(); j++){
            grade.push_back(column_grades[i][j]);
        }
        SignatureColumn column(grade, i);
        for(size_t j=0; j<matrix[i].size(); j++){
            if(matrix[i][j] != 0){
                column.push(column_entry_t(1, j)); // TODO: implement support for other modulus than Z/2Z.
            }
        }
        column.syzygy.push(column_entry_t(1, i));
        M.push_back(column);
    }
    return M;
}

GradedMatrix translateOutputMatrix(Matrix& matrix){
    GradedMatrix graded_matrix;
    for(SignatureColumn column : matrix){
        std::vector<std::pair<int, int>> sparse_column;
        while(!column.empty()){
            column_entry_t entry = column.pop_pivot();
            if(entry.get_index() != -1){
                sparse_column.push_back(std::pair<int, int>(entry.get_value(), entry.get_index()));
            }
        }
        graded_matrix.matrix.push_back(sparse_column);
        std::vector<int> column_grade;
        for(size_t grade_index=0; grade_index<column.grade.size(); grade_index++){
            column_grade.push_back((int)column.grade[grade_index]);
        }
        graded_matrix.column_grades.push_back(column_grade);
    }
    return graded_matrix;
}

GradedMatrix compute_kernel(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    Matrix kernel;
    if(column_grades.size()>0 && column_grades[0].size()==2){
        kernel = computeKernel_2p(M);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
        kernel = gbs.first;
    }
    GradedMatrix Ker = translateOutputMatrix(kernel);
    Ker.row_grades = column_grades;
    return Ker;
}

std::pair<GradedMatrix, GradedMatrix> groebner_bases(std::vector<std::vector<int>>& matrix, std::vector<std::vector<int>>& row_grades, std::vector<std::vector<int>>& column_grades){
    Matrix M = translateInputMatrix(matrix, column_grades);
    std::pair<Matrix, Matrix> gbs = computeGroebnerBases_gradeopt_min(M);
    GradedMatrix Im = translateOutputMatrix(gbs.first);
    GradedMatrix Ker = translateOutputMatrix(gbs.second);
    Ker.row_grades = column_grades;
    Im.row_grades = row_grades;
    return std::pair<GradedMatrix, GradedMatrix>(Im, Ker);
}

GradedMatrix presentation_FIrep(std::vector<std::vector<int>>& high_matrix, std::vector<std::vector<int>>& column_grades_h, std::vector<std::vector<int>>& low_matrix, std::vector<std::vector<int>>& column_grades_l){
    Matrix M_h = translateInputMatrix(high_matrix, column_grades_h);
    Matrix M_l = translateInputMatrix(low_matrix, column_grades_l);
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation_output;
    if(M_l.size()>0 && M_l[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(M_l);
        presentation_output = compute_presentation_2p(M_h, kernel);
    }else{
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(M_h, M_l);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "Time elapsed: " << duration.count() << " seconds" << std::endl;
    GradedMatrix graded_matrix = translateOutputMatrix(presentation_output.first);
    for(std::pair<size_t, grade_t> entry : presentation_output.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices_dm(distance_matrices, max_metric_values, filters, hom_dim);
    verify_kernel(boundary_matrices.first, boundary_matrices.second);
    for(size_t i=0; i<boundary_matrices.second.size(); i++){
        boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
    }
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        std::pair<Matrix, Matrix> gbs = computeGroebnerBases(boundary_matrices.second);
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(boundary_matrices.first, gbs.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

GradedMatrix presentation(std::vector<std::vector<input_t>>& _points, std::vector<int>& _metrics, std::vector<input_t>& _max_metric_values, std::vector<int>& _filters, int hom_dim){
    std::vector<Metric*> metrics;
    for(size_t i=0; i<_metrics.size(); i++){
        metrics.push_back(parse_metric(_metrics[i]));
    }
    std::vector<Filter*> filters;
    for(size_t i=0; i<_filters.size(); i++){
        filters.push_back(parse_filter(_filters[i]));
    }
    std::vector<input_t> max_metric_values;
    for(size_t i=0; i<_metrics.size(); i++){
        max_metric_values.push_back(_max_metric_values[i]);
    }
    std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(_points, metrics, filters, max_metric_values, hom_dim);
    verify_kernel(boundary_matrices.first, boundary_matrices.second);
    for(auto& metric : metrics){
        delete metric;
    }
    for(auto& filter : filters){
        delete filter;
    }
    for(size_t i=0; i<boundary_matrices.second.size(); i++){
        boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
    }
    Matrix input_copy;
    for(auto& column : boundary_matrices.second){
        input_copy.push_back(SignatureColumn(column.grade, 0, column));
    }
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    if(boundary_matrices.second.size()>0 && boundary_matrices.second[0].grade.size()==2){
        Matrix kernel = computeKernel_2p(boundary_matrices.second);
        presentation = compute_presentation_2p(boundary_matrices.first, kernel);
    }else{
        std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(boundary_matrices.first, boundary_matrices.second);
    }
    GradedMatrix graded_matrix = translateOutputMatrix(presentation.first);
    for(std::pair<size_t, grade_t> entry : presentation.second){
        std::vector<int> row_grade;
        for(size_t grade_index=0; grade_index<entry.second.size(); grade_index++){
            row_grade.push_back((int)entry.second[grade_index]);
        }
        graded_matrix.row_grades.push_back(row_grade);
    }
    return graded_matrix;
}

std::vector<int> translate_grade(grade_t _grade) {
    std::vector<int> grade;
    for ( auto& e : _grade ) {
        grade.push_back((int)e);
    }
    return grade;
}

PythonCompressedLandscape landscapes_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, input_t& max_metric_value, int hom_dim){
    Metric* m = new SquaredEuclideanMetric();
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) = compute_boundary_matrices_spatiotemporal(trajectories, m, max_metric_value, hom_dim);

    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(high_matrix, low_matrix);

    CompressedLandscape compressed_landscape = computeCompressedLandscape(minimal_presentation.first, minimal_presentation.second);

    PythonCompressedLandscape python_compressed_landscape;
    for ( auto& slice : compressed_landscape ) {
        std::pair<std::vector<int>, std::vector<std::vector<int>>> translated_slice;
        std::vector<int> g = translate_grade(slice.first);
        std::vector<std::vector<int>> f;
        for ( auto& l : slice.second ) {
            f.push_back(translate_grade(l));
        }
        python_compressed_landscape.pairings.push_back(std::pair<std::vector<int>, std::vector<std::vector<int>>>(g, f));
    }
    python_compressed_landscape.index_value_lists = index_value_lists;
    return python_compressed_landscape;
}

PythonLandscape landscapes_spatiotemporal_naive(std::vector<std::vector<std::vector<input_t>>>& trajectories, input_t& max_metric_value, int hom_dim, int landscape_dim){
    Metric* m = new SquaredEuclideanMetric();
    Matrix high_matrix; Matrix low_matrix;
    std::vector<std::vector<input_t>> index_value_lists;
    std::tie(high_matrix, low_matrix, index_value_lists) = compute_boundary_matrices_spatiotemporal(trajectories, m, max_metric_value, hom_dim);

    grade_t v_max = low_matrix[0].grade;
    for ( auto& column : low_matrix ) {
        v_max = column.grade.join(v_max);
    }
    for ( auto& column : high_matrix ) {
        v_max = column.grade.join(v_max);
    }

    std::pair<Matrix, hash_map<size_t, grade_t>> presentation;
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(high_matrix, low_matrix);

    Landscape landscape = computeLandscapeNaive(minimal_presentation.first, minimal_presentation.second, v_max, landscape_dim);

    PythonLandscape python_landscape;
    for ( auto& slice : landscape ) {
        std::vector<int> g = translate_grade(slice.first);
        int f = (int)slice.second;
        python_landscape.landscape.push_back(std::pair<std::vector<int>, int>(g, f));
    }
    return python_landscape;
}

/*
 Testing
 */

std::vector<double> run_gb_singatures_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, std::vector<grade_t>> minimal_presentation = computeMinimalPresentation_3p(image, input_matrix, debug);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << minimal_presentation.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = minimal_presentation.first.size();
    return std::vector<double>{time, res_memory, size};
}

std::vector<double> run_gb_shreyer_pres(Matrix& input_matrix, Matrix& image, bool debug=false){
    auto start = std::chrono::high_resolution_clock::now();
    Matrix kernel = computeGroebnerBases_gradeopt(input_matrix).second;
    std::pair<Matrix, hash_map<size_t, grade_t>> presentation = compute_presentation_schreyer(image, kernel);
    auto stop = std::chrono::high_resolution_clock::now();
    auto b_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " Presentation size: " << presentation.first.size() << std::endl;//gbs2.size() << std::endl;
    double time = b_time.count();
    double size = presentation.first.size();
    return std::vector<double>{time, res_memory, size};
}


/* Input */

void print_usage_and_exit(int exit_code) {
    std::cerr
    << "Usage: "
    << "mph "
    << "[options] [filename]" << std::endl
    << std::endl
    << "Options:" << std::endl
    << std::endl
    << "  --help           print this screen" << std::endl
    << "  --dim <k>        compute presentation matrix of the k-th persistent homology module" << std::endl
    << "  --firep <k>      input file with rivet <firep> file format" << std::endl
    << std::endl;
    exit(exit_code);
}



int main(int argc, char** argv) {
    const char* filename = nullptr;

    int dim_max = 1;
    bool firep = false;

    if(argc==1){
        print_usage_and_exit(0);
    }

    for (index_t i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage_and_exit(0);
        } else if (arg == "--dim") {
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            dim_max = (int)std::stol(parameter, &next_pos);
            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        } else if (arg == "--firep"){
            firep = true;
        } else {
            if (filename) { print_usage_and_exit(-1); }
            filename = argv[i];
        }
    }

    std::pair<Matrix, std::vector<grade_t>> presentation;
    if(firep){
        Matrix high_matrix, low_matrix;
        std::ifstream file_stream(filename);
        if (file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }
        read_input_file<SyzColumn, SyzColumn>(file_stream, high_matrix, low_matrix);
        presentation = computeMinimalPresentation_3p(high_matrix, low_matrix);
    } else{
        std::ifstream file_stream(filename);
        if (filename && file_stream.fail()) {
            std::cerr << "couldn't open file " << filename << std::endl;
            exit(-1);
        }

        // TODO: how should choice of metrics and filters be specified in input?

        std::vector<Metric*> metrics;
        metrics.push_back(new SquaredEuclideanMetric());
        std::vector<Filter*> filters;
        filters.push_back(new Filter());
        std::vector<input_t> max_metric_values;
        max_metric_values.push_back(10);
        std::vector<std::vector<input_t>> points = read_point_cloud(file_stream, 0);
        std::pair<Matrix, Matrix> boundary_matrices = compute_boundary_matrices(points, metrics, filters, max_metric_values, dim_max);
        verify_kernel(boundary_matrices.first, boundary_matrices.second);
        for(auto& metric : metrics){
            delete metric;
        }
        for(auto& filter : filters){
            delete filter;
        }
        for(size_t i=0; i<boundary_matrices.second.size(); i++){
            boundary_matrices.second[i].syzygy.push(column_entry_t(1, i));
        }
        presentation = computeMinimalPresentation_3p(boundary_matrices.first, boundary_matrices.second);
    }


    exit(0);
}
