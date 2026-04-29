#include "landscapes.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_set>

#include "matrix.h"
#include "signatureColumn.h"
#include "utils.h"

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
