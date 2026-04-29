#include "presentation.h"

#include <queue>
#include <set>
#include <algorithm>


hash_map<size_t, size_t> compute_local_pairs(Matrix& columns, hash_map<size_t, grade_t>& row_grades){
    /** Computes the local positive and negative pairs of columns and pivots.

     Arguments:
     columns {Matrix} -- the matrix with columns to be labeled local positive, negative or global.
     row_grades {std::vector<grade_t>} -- a list of the multigrade of each row.

     Returns:
     hash_map<size_t, size_t> -- a map that sends a local positive row to a local negative column.

     */
    hash_map<size_t, size_t> pairs;
    for(size_t column_index=0; column_index < columns.size(); column_index++){
        if(columns[column_index].local != 1 && columns[column_index].get_pivot().get_index() != -1 && row_grades[columns[column_index].get_pivot().get_index()] == columns[column_index].grade){

            if(pairs.find(columns[column_index].get_pivot().get_index()) != pairs.end()){
                SignatureColumn working_column(columns[column_index].grade, -1);
                size_t column_to_add = column_index;
                while(true){
                    working_column.plus(columns[column_to_add]);
                    index_t pivot = working_column.get_pivot().get_index();
                    if(pivot != -1 && row_grades[pivot] == working_column.grade){
                        if (pairs.find(pivot) != pairs.end()) {
                            column_to_add = pairs[pivot];
                        } else {
                            pairs[pivot] = column_index;
                            columns[column_index].local = -1;
                            break;
                        }
                    } else{
                        break;
                    }
                }
            }else{
                pairs[columns[column_index].get_pivot().get_index()] = column_index;
                columns[column_index].local = -1;
            }
        }
    }
    return pairs;
}

Matrix compute_global_columns(Matrix& columns, hash_map<size_t, size_t>& positive_pairs, std::set<size_t>& negative_rows){
    /** Performs the presentation minimization step by removing local positive and negative columns.

     Arguments:
     columns {Matrix} -- a matrix with columns labeled local positive, negative or global.
     positive_pairs {hash_map<size_t, size_t>} -- a map sending local positive rows to local negative columns.
     negative_rows {std::set<size_t>} -- a set containing the indices of local negative rows.

     Returns:
     Matrix -- a minimized presentation matrix.

     */

    Matrix global_columns;
    for(size_t index_column_to_reduce = 0; index_column_to_reduce<columns.size();index_column_to_reduce++) {
        if(columns[index_column_to_reduce].local != 0)
            continue;
        SignatureColumn working_boundary(columns[index_column_to_reduce].grade, index_column_to_reduce);
        SignatureColumn global_column(columns[index_column_to_reduce].grade, index_column_to_reduce);
        working_boundary.plus(columns[index_column_to_reduce]);
        while(true) {
            index_t pivot = working_boundary.get_pivot().get_index();
            if (pivot != -1) {
                if (negative_rows.find(pivot) != negative_rows.end()) {
                    working_boundary.pop_pivot();
                }else if(positive_pairs.find(pivot) != positive_pairs.end()) {
                    working_boundary.plus(columns[positive_pairs[pivot]]);
                }else{
                    global_column.push(working_boundary.pop_pivot());
                }
            }else{
                break;
            }
        }
        if(global_column.get_pivot().get_index() != -1){
            global_columns.push_back(global_column);
        }
    }
    return global_columns;
}


std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_2p(Matrix& image_generators, Matrix& generating_set_kernel){
    /** Computes a minimal presentation of the ith homology module using generators of the kernel and image of
     the boundary matrices \Delta_{i-1} and \Delta_i respectively.

     Arguments:
     image_generators {Matrix} -- a generating set of the image of the boundary map \Delta_i.
     kernel_generators {Matrix} -- a Groebner basis of the kernel of the boundary matrix \Delta_{i-1}.

     Returns:
     Matrix -- a matrix describing a minimal presentation of the ith homology module.

     */

    Matrix generating_set_image = compute_minimal_generating_set(image_generators);

    //Sort image and kernel gens in colex and lex respectively
   sort(generating_set_kernel.begin(), generating_set_kernel.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade) ;
         });
     sort(generating_set_image.begin(), generating_set_image.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });

    hash_map<size_t, grade_t> row_grade_map;
    for(size_t j=0; j<generating_set_kernel.size(); j++){
        row_grade_map[j] = grade_t(generating_set_kernel[j].grade);
    }
    for(size_t i=0; i<generating_set_kernel.size(); i++){
        generating_set_kernel[i].syzygy = SyzColumn();
        generating_set_kernel[i].syzygy.push(column_entry_t(1, i));
    }


    hash_map<size_t, size_t> index_map_low, index_map_high;
    for(size_t i=0; i<generating_set_kernel.size(); i++){
        generating_set_kernel[i].signature_index = i;
        index_map_high[generating_set_kernel[i].grade.back()] = i;
        if(index_map_low.find(generating_set_kernel[i].grade.back()) == index_map_low.end()){
            index_map_low[generating_set_kernel[i].grade.back()] = i;
        }
    }


    /* Main algorithm that iterates through the index set */
    index_t max_pivot=0;
    for(auto& generator : generating_set_kernel){
        if(max_pivot < generator.get_pivot().get_index()){
            max_pivot = generator.get_pivot().get_index();
        }
    }
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    Matrix presentation_matrix;
    presentation_matrix.reserve(generating_set_image.size());
    for(auto& column : generating_set_image){
        grade_t& v = column.grade;

        /* Initialize Macaulay matrix */
        for( size_t i=0; i<=index_map_high[v.back()]; i++ ){
            if(generating_set_kernel[i].size()>0 && generating_set_kernel[i].grade[0] <= v[0]){
                pivot = generating_set_kernel[i].get_pivot().get_index();
                if(pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                    SignatureColumn& working_column = generating_set_kernel[i];
                    while(pivot != -1 && pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                        working_column.plus(generating_set_kernel[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        pivot_map[pivot] = i;
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }
            }
        }


        // Reduce the image column with the pivot set.
        column.syzygy = SyzColumn();
        while(true){
            pivot = column.get_pivot().get_index();
            if(pivot != -1){
                if(pivot_map[pivot] == -1){
                    throw std::runtime_error("Failed to express image column in terms of kernel generating set.");
                }
                column.plus(generating_set_kernel[pivot_map[pivot]]);
            }else{
                break;
            }
        }
        presentation_matrix.push_back(SignatureColumn(column.grade, column.signature_index, column.syzygy));
        presentation_matrix.back().refresh();
    }


    hash_map<size_t, size_t> pairs = compute_local_pairs(presentation_matrix, row_grade_map);

    std::set<size_t> negative_columns;
    presentation_matrix = compute_global_columns(presentation_matrix, pairs, negative_columns);

    for(auto& entry : pairs){
        if(row_grade_map.find(entry.first) != row_grade_map.end()){
            row_grade_map.erase(entry.first);
        }
    }



    // Reindex the rows of the presentation matrix to account for deleted rows

    std::vector<index_t> rows;
    for(auto& entry : row_grade_map){
        rows.push_back(entry.first);
    }
    sort(rows.begin(), rows.end());
    hash_map<size_t, size_t> row_index_map;
    for(size_t i=0; i<rows.size(); i++){
        row_index_map[rows[i]] = i;
    }

    Matrix minimized_presentation;
    for(size_t i=0;i<presentation_matrix.size(); i++){
        SignatureColumn column(presentation_matrix[i].grade, presentation_matrix[i].signature_index);
        for(auto& entry : presentation_matrix[i]){
            if(row_grade_map.find(entry) != row_grade_map.end()){
                column.push(column_entry_t(entry, row_index_map[entry]));
            }
        }
        minimized_presentation.push_back(column);
    }

    return std::pair<Matrix, hash_map<size_t, grade_t>>(minimized_presentation, row_grade_map);
}



std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation_3p(Matrix& image_columns, Matrix& columns, bool debug){
    /*
     Optimized for 3-parameter modules. Uses the fact the module has projective dimension 3 => syzygies of Z are free. A consequence is that a syzygy is either immideately reduced to zero or never reduced to zero within each 2-dim slice.

     The main function computing a minimal presentation of the homology of the chain complex

                                    C ->^F A ->^G B

     where the columns of F are given by 'image_columns' and the columns of G by 'columns'.

     Arguments:
     image_columns {std::vector<SignatureColumn>} -- a minimal generating set of the image map.
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded modules.

     Returns:
     Matrix -- the minimal presentation matrix.
     hash_map<size_t, grade_t> -- a map of the row grades.
     */


    /* Sort columns colexicographically */
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
    }
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    hash_map<size_t, size_t> kernel_reorder_map;
    for(size_t j=0; j<columns.size(); j++){
        kernel_reorder_map[columns[j].signature_index] = j;
    }

    // Reindex image columns
    for(size_t i=0; i<image_columns.size(); i++) {
        SignatureColumn image_copy(image_columns[i]);
        SignatureColumn working_column(image_columns[i].grade, i);
        while(!image_copy.empty()){
            column_entry_t entry = image_copy.pop_pivot();
            working_column.push(column_entry_t(1, kernel_reorder_map[entry.get_index()]));
        }
        image_columns[i] = working_column;
    }
    sort(image_columns.begin(), image_columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade<rhs.grade;
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        columns[i].syzygy = SyzColumn();
        columns[i].syzygy.push(column_entry_t(1, i));
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::vector<std::vector<grade_t>> grade_listsH;

    /* Vectors to store the columns of the GBs */
    Matrix gb_columns, gb_columnsH, syzygiesH;


    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB, GBH;
    std::vector<std::vector<signature_t>> Syz, SyzH;

    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        grades.push(column.grade);
        grade_listsH.push_back(std::vector<grade_t>());
    }
    for(auto& column : image_columns){
        grades.push(column.grade);
    }
    index_t pivot;

    std::vector<index_t> pivot_map(max_pivot+1, -1);
    std::vector<index_t> pivot_mapH(columns.size(), -1);

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    std::vector<size_t> L_F;
    std::vector<signature_t> reorder_Z_map;

    size_t image_index=0;

    int iter_index = 0;
    while(!grades.empty()){
        grade_t v = grades.top();
        while(v == grades.top()){
            grades.pop();
        }

        iter_index++;
        size_t grade_hash = grade_hasher(v);

        // B^h_z
        for( auto& signature : reorder_Z_map ){
            size_t i = signature.get_index();
            column_index = -1;
            //bool in_syz = false;
            //size_t syz_size = SyzH[i].size();
            if(SyzH[i].size()==0 && GBH[i][0].get_grade().leq_poset(v)){
                size_t gb_size = GBH[i].size();
                column_index = (int)GBH[i][0].get_index();
                for(size_t j=gb_size-1; j >= 1; j--){
                    if(GBH[i][j].get_grade().leq_poset(v)){
                        column_index = (int)GBH[i][j].get_index();
                        break;
                    }
                }
                /*for(size_t j=syz_size-1; j < syz_size; j--){
                    if((SyzH[i][j].get_grade()).leq_poset(v)){
                        in_syz = true;
                        break;
                    }
                }*/
            }

            /*if(!in_syz){
                size_t gb_size = GBH[i].size();
                for(size_t j=gb_size-1; j < gb_size; j--){
                    if(GBH[i][j].get_grade().leq_poset(v)){
                        column_index = (int)GBH[i][j].get_index();
                        break;
                    }
                }
            }*/


            if(column_index > -1){
                pivot = gb_columnsH[column_index].get_pivot_index();
                if(pivot > -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columnsH[column_index];
                    while(pivot != -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        gb_columnsH.push_back(working_column);
                        GBH[i].push_back(signature_t(grade, gb_columnsH.size()-1));
                        pivot_mapH[pivot] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(grade == v){
                            if(grade_listsH[pivot].size() == 0 || grade != grade_listsH[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                                    grade_t m_ji = grade_listsH[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_listsH[pivot][j]));
                                    }
                                }
                                grade_listsH[pivot].push_back(grade);
                            }
                        }
                    }else{
                        working_column.syzygy.refresh();
                            syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
                            SyzH[i].push_back(signature_t(working_column.grade, syzygiesH.size()-1));
                            syzygiesH[syzygiesH.size()-1].last_updated = iter_index;
                      /*  }else{
                            SyzH[i].push_back(signature_t(working_column.grade, syzygiesH.size()-1));
                        }*/
                    }
                } else{
                    if(pivot != -1){
                        pivot_mapH[pivot] = column_index;
                        gb_columnsH[column_index].last_updated = iter_index;
                    }
                }
            }
        }

        /* Reduce image columns */
        while(image_index < image_columns.size() && image_columns[image_index].grade == v){
            SignatureColumn working_column(image_columns[image_index]);
            working_column.syzygy = SyzColumn();
            pivot = working_column.get_pivot_index();
            while(pivot != -1 && pivot_mapH[pivot] != -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                pivot = working_column.get_pivot_index();
            }
            working_column.syzygy.refresh();
            if(pivot != -1){
                working_column.refresh();
                working_column.signature_index = GBH.size();
                L_F.push_back(GBH.size()); // The set keeping track of which kernel columns to remove.
                GBH.push_back(std::vector<signature_t>());
                SyzH.push_back(std::vector<signature_t>());
                gb_columnsH.push_back(working_column);
                GBH[GBH.size()-1].push_back(signature_t(working_column.grade, gb_columnsH.size()-1));
                insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                pivot_mapH[pivot] = gb_columnsH.size()-1;
                gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                if(grade_listsH[pivot].size() == 0 || working_column.grade != grade_listsH[pivot].back()){
                    std::vector<grade_t> minimal_elements_tmp;
                    for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                        grade_t m_ji = grade_listsH[pivot][j].m_ji(working_column.grade);
                        bool is_minimal = true;
                        for(auto& el : minimal_elements_tmp){
                            if(el.leq_poset_m(m_ji)){
                                is_minimal = false;
                                break;
                            }
                        }
                        if(is_minimal){
                            minimal_elements_tmp.push_back(m_ji);
                            grades.push(working_column.grade.join(grade_listsH[pivot][j]));
                        }
                    }
                    grade_listsH[pivot].push_back(working_column.grade);
                }
            }else{
                syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
            }
            image_index++;
        }

        /* Reduce columns of map G */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                if(pivot_mapH[i] == -1 || gb_columnsH[pivot_mapH[i]].last_updated != iter_index){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(pivot_mapH[i] == -1 && grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        if(is_new){
                            gb_columns[column_index] = working_column;
                            GB[i].push_back(signature_t(grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }else{
                            /*if(GB[i].size()==1){
                                gb_columns.push_back(working_column);
                                GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            }else{
                                if(pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] == GB[i][1].get_index()){
                                    pivot_map[gb_columns[GB[i][1].get_index()].get_pivot_index()] = -1;
                                }
                                gb_columns[GB[i][1].get_index()] = working_column;
                                GB[i][1].first = grade;
                            }*/
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }
                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        SignatureColumn syz(working_column.grade, GBH.size(), working_column.syzygy);
                        syz.syzygy.push(column_entry_t(1, GBH.size()));
                        GBH.push_back(std::vector<signature_t>());
                        SyzH.push_back(std::vector<signature_t>());
                        gb_columnsH.push_back(syz);
                        GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                        insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                        pivot_mapH[i] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[i].size(); j++){
                                    grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(working_column.grade.join(grade_listsH[i][j]));
                                    }
                                }
                                grade_listsH[i].push_back(working_column.grade);
                            }
                        }
                    }
                } else{
                    if(pivot != -1){
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(gb_columns[column_index].grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(gb_columns[column_index].grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }else{
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }
                    } else {
                        if(is_new){
                            SignatureColumn working_column = gb_columns[column_index];
                            SyzColumn syz_column = SyzColumn();
                            syz_column.push(column_entry_t(1, column_index));
                            SignatureColumn syz(working_column.grade, GBH.size(), syz_column);
                            syz.syzygy.push(column_entry_t(1, GBH.size()));
                            GBH.push_back(std::vector<signature_t>());
                            SyzH.push_back(std::vector<signature_t>());
                            gb_columnsH.push_back(syz);
                            GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                            insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                            pivot_mapH[i] = gb_columnsH.size()-1;
                            gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                            if(working_column.grade == v){
                                if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                    std::vector<grade_t> minimal_elements_tmp;
                                    for(size_t j=0; j<grade_listsH[i].size(); j++){
                                        grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                        bool is_minimal = true;
                                        for(auto& el : minimal_elements_tmp){
                                            if(el.leq_poset_m(m_ji)){
                                                is_minimal = false;
                                                break;
                                            }
                                        }
                                        if(is_minimal){
                                            minimal_elements_tmp.push_back(m_ji);
                                            grades.push(working_column.grade.join(grade_listsH[i][j]));
                                        }
                                    }
                                    grade_listsH[i].push_back(working_column.grade);
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    Matrix filtered_gb_columnsH;
    size_t j = 0;
    for(size_t i=0; i<GBH.size(); i++){
        if(j<L_F.size() && i == L_F[j]){
            j++;
        } else {
            filtered_gb_columnsH.push_back(gb_columnsH[GBH[i][0].get_index()]);
        }
    }

    sort(filtered_gb_columnsH.begin(), filtered_gb_columnsH.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });

    hash_map<size_t, size_t> signature_map;
    for(size_t j=0; j<filtered_gb_columnsH.size(); j++){
        signature_map[filtered_gb_columnsH[j].signature_index] = j;
    }

    // Reindex syzygiesH
    Matrix reindexed_syzygies;
    for(size_t i=0; i<syzygiesH.size(); i++) {
        SignatureColumn sig_copy(syzygiesH[i]);
        SignatureColumn working_column(syzygiesH[i].grade, i);
        while(!sig_copy.empty()){
            column_entry_t entry = sig_copy.pop_pivot();
            working_column.push(column_entry_t(1, signature_map[entry.get_index()]));
        }
        reindexed_syzygies.push_back(working_column);
    }

    reindexed_syzygies = compute_minimal_generating_set(reindexed_syzygies);

    std::vector<grade_t> row_grades;
    for(auto& column : filtered_gb_columnsH) {
        row_grades.push_back(column.grade);
    }

    return std::pair<Matrix, std::vector<grade_t>>(reindexed_syzygies, row_grades);
}

std::pair<Matrix, hash_map<size_t, grade_t>> compute_presentation_schreyer(Matrix& image_generators, Matrix& kernel_generators, bool debug){
    /** Input a minimal set of generators for the image and a Groebner basis for the kernel.

     */

    sort(kernel_generators.begin(), kernel_generators.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });

    for(size_t i=0; i<kernel_generators.size(); i++){
        kernel_generators[i].signature_index = i;
        kernel_generators[i].syzygy = SyzColumn();
        kernel_generators[i].syzygy.push(column_entry_t(1, i));
    }

    /* Translate image matrix */
    hash_map<index_t, std::vector<size_t>> pivots_indices;
    pivots_indices.reserve(kernel_generators.size());
    for(auto& column : kernel_generators){
        if(pivots_indices.find(column.get_pivot_index()) == pivots_indices.end()){
            pivots_indices[column.get_pivot_index()] = std::vector<size_t>();
        }
        pivots_indices[column.get_pivot_index()].push_back(column.signature_index);
    }

    Matrix translated_image_columns;
    for(auto& column : image_generators){
        SignatureColumn working_column(column);
        working_column.syzygy = SyzColumn();
        index_t pivot = working_column.get_pivot_index();
        while(pivot != -1){
            bool has_eliminated = false;
            std::vector<SignatureColumn> pivot_columns;
            for(auto& pivot_index : pivots_indices[pivot]){
                pivot_columns.push_back(kernel_generators[pivot_index]);
                if(kernel_generators[pivot_index].grade.leq_poset(working_column.grade)){
                    working_column.plus(kernel_generators[pivot_index]);
                    has_eliminated = true;
                    break;
                }
            }
            if(!has_eliminated){
                throw std::runtime_error("Cannot express image column in terms of kernel.");
            }
            pivot = working_column.get_pivot_index();
        }
        translated_image_columns.push_back(SignatureColumn(working_column.grade, working_column.signature_index, working_column.syzygy));
    }

    Matrix syzygy_module = compute_syzygy_module(kernel_generators);
    for(auto& column : translated_image_columns){
        syzygy_module.push_back(column);
    }

    std::pair<Matrix, Matrix> m_pair = compute_minimal_generating_set2(kernel_generators);
    Matrix& generating_set_kernel = m_pair.first;
    Matrix& change_of_basis_map = m_pair.second;


    hash_map<size_t, grade_t> row_grade_map;
    for(size_t j=0; j<generating_set_kernel.size(); j++){
        row_grade_map[j] = generating_set_kernel[j].grade;
    }

    Matrix presentation;
    for(auto& syz_column : syzygy_module){
        SignatureColumn column(syz_column.grade, syz_column.signature_index);
        index_t pivot = syz_column.get_pivot().get_index();
        while(pivot != -1){
            column.plus(change_of_basis_map[pivot]);
            syz_column.pop_pivot();
            pivot = syz_column.get_pivot().get_index();
        }
        if(debug){
            if(syz_column.grade != column.grade){
                throw std::runtime_error("Grade of columns has changed.");
            }
            SignatureColumn ker_column(syz_column.grade, syz_column.signature_index);
            SignatureColumn copy(column);
            pivot = copy.get_pivot().get_index();
            while(pivot != -1){
                ker_column.plus(generating_set_kernel[pivot]);
                copy.pop_pivot();
                pivot = copy.get_pivot().get_index();
            }
            if(ker_column.get_pivot().get_index() != -1){
                throw std::runtime_error("Column non-zero");
            }
        }
        if(column.get_pivot().get_index() != -1){
            presentation.push_back(column);
        }
    }


    // Compute minimal generating set from syz-module and image columns
    sort(presentation.begin(), presentation.end(),[ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade < rhs.grade ;
         });

    hash_map<size_t, size_t> pairs = compute_local_pairs(presentation, row_grade_map);

    std::set<size_t> negative_columns;
    presentation = compute_global_columns(presentation, pairs, negative_columns);

    for(auto& entry : pairs){
        if(row_grade_map.find(entry.first) != row_grade_map.end()){
            row_grade_map.erase(entry.first);
        }
    }

    std::vector<index_t> rows;
    for(auto& entry : row_grade_map){
        rows.push_back(entry.first);
    }
    sort(rows.begin(), rows.end());
    hash_map<size_t, size_t> row_index_map;
    for(size_t i=0; i<rows.size(); i++){
        row_index_map[rows[i]] = i;
    }

    for(size_t i=0;i<presentation.size(); i++){
        SignatureColumn column(presentation[i].grade, presentation[i].signature_index);
        for(auto& entry : presentation[i]){
            if(row_grade_map.find(entry) != row_grade_map.end()){
                column.push(column_entry_t(entry, row_index_map[entry]));
            }
        }
        presentation[i].swap(column);
    }
    presentation = compute_minimal_generating_set(presentation);
    return std::pair<Matrix, hash_map<size_t, grade_t>>(presentation, row_grade_map);
}

std::pair<Matrix, std::vector<grade_t>> computeMinimalPresentation(Matrix& image_columns, Matrix& columns, bool debug){
    /*
     The main function computing a minimal presentation of the homology of the chain complex

                                    C ->^F A ->^G B

     where the columns of F are given by 'image_columns' and the columns of G by 'columns'.

     Arguments:
     image_columns {std::vector<SignatureColumn>} -- a minimal generating set of the image map.
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded modules.

     Returns:
     Matrix -- the minimal presentation matrix.
     hash_map<size_t, grade_t> -- a map of the row grades.
     */


    /* Sort columns colexicographically */
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
    }
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    hash_map<size_t, size_t> kernel_reorder_map;
    for(size_t j=0; j<columns.size(); j++){
        kernel_reorder_map[columns[j].signature_index] = j;
    }

    // Reindex image columns
    for(size_t i=0; i<image_columns.size(); i++) {
        SignatureColumn image_copy(image_columns[i]);
        SignatureColumn working_column(image_columns[i].grade, i);
        while(!image_copy.empty()){
            column_entry_t entry = image_copy.pop_pivot();
            working_column.push(column_entry_t(1, kernel_reorder_map[entry.get_index()]));
        }
        image_columns[i] = working_column;
    }
    sort(image_columns.begin(), image_columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade<rhs.grade;
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        columns[i].syzygy = SyzColumn();
        columns[i].syzygy.push(column_entry_t(1, i));
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::vector<std::vector<grade_t>> grade_listsH;

    /* Vectors to store the columns of the GBs */
    Matrix gb_columns, gb_columnsH, syzygiesH;


    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB, GBH;
    std::vector<std::vector<signature_t>> Syz, SyzH;

    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot_index()){
            max_pivot = column.get_pivot_index();
        }
    }
    for(size_t i=0; i<=max_pivot; i++){
        grade_lists.push_back(std::vector<grade_t>());
    }
    for(auto& column : columns){
        grades.push(column.grade);
        grade_listsH.push_back(std::vector<grade_t>());
    }
    for(auto& column : image_columns){
        grades.push(column.grade);
    }
    index_t pivot;

    std::vector<index_t> pivot_map(max_pivot+1, -1);
    std::vector<index_t> pivot_mapH(columns.size(), -1);

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    std::vector<size_t> L_F;
    std::vector<signature_t> reorder_Z_map;

    size_t image_index=0;
    int iter_index = 0;
    while(!grades.empty()){
        grade_t v = grades.top();
        while(v == grades.top()){
            grades.pop();
        }

        iter_index++;
        size_t grade_hash = grade_hasher(v);

        // B^h_z
        for( auto& signature : reorder_Z_map ){
            size_t i = signature.get_index();
            column_index = -1;
            bool in_syz = false;
            size_t syz_size = SyzH[i].size();
            if(syz_size>0 && GBH[i][0].get_grade().leq_poset(v)){
                for(size_t j=syz_size-1; j < syz_size; j--){
                    if((SyzH[i][j].get_grade()).leq_poset(v)){
                        in_syz = true;
                        break;
                    }
                }
            }

            if(!in_syz){
                size_t gb_size = GBH[i].size();
                for(size_t j=gb_size-1; j < gb_size; j--){
                    if(GBH[i][j].get_grade().leq_poset(v)){
                        column_index = (int)GBH[i][j].get_index();
                        break;
                    }
                }
            }


            if(column_index > -1){
                pivot = gb_columnsH[column_index].get_pivot_index();
                if(pivot > -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columnsH[column_index];
                    while(pivot != -1 && pivot_mapH[pivot] > -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        gb_columnsH.push_back(working_column);
                        GBH[i].push_back(signature_t(grade, gb_columnsH.size()-1));
                        pivot_mapH[pivot] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(grade == v){
                            if(grade_listsH[pivot].size() == 0 || grade != grade_listsH[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                                    grade_t m_ji = grade_listsH[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_listsH[pivot][j]));
                                    }
                                }
                                grade_listsH[pivot].push_back(grade);
                            }
                        }
                    }else{
                        working_column.syzygy.refresh();
                        syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
                        SyzH[i].push_back(signature_t(working_column.grade, syzygiesH.size()-1));
                        syzygiesH[syzygiesH.size()-1].last_updated = iter_index;
                    }
                } else{
                    if(pivot != -1){
                        pivot_mapH[pivot] = column_index;
                        gb_columnsH[column_index].last_updated = iter_index;
                    }
                }
            }
        }

        /* Reduce image columns */
        while(image_index < image_columns.size() && image_columns[image_index].grade == v){
            SignatureColumn working_column(image_columns[image_index]);
            working_column.syzygy = SyzColumn();
            pivot = working_column.get_pivot_index();
            while(pivot != -1 && pivot_mapH[pivot] != -1 && gb_columnsH[pivot_mapH[pivot]].last_updated == iter_index){
                working_column.plus(gb_columnsH[pivot_mapH[pivot]]);
                pivot = working_column.get_pivot_index();
            }
            working_column.syzygy.refresh();
            if(pivot != -1){
                working_column.refresh();
                working_column.signature_index = GBH.size();
                L_F.push_back(GBH.size()); // The set keeping track of which kernel columns to remove.
                GBH.push_back(std::vector<signature_t>());
                SyzH.push_back(std::vector<signature_t>());
                gb_columnsH.push_back(working_column);
                GBH[GBH.size()-1].push_back(signature_t(working_column.grade, gb_columnsH.size()-1));
                insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                pivot_mapH[pivot] = gb_columnsH.size()-1;
                gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                if(grade_listsH[pivot].size() == 0 || working_column.grade != grade_listsH[pivot].back()){
                    std::vector<grade_t> minimal_elements_tmp;
                    for(size_t j=0; j<grade_listsH[pivot].size(); j++){
                        grade_t m_ji = grade_listsH[pivot][j].m_ji(working_column.grade);
                        bool is_minimal = true;
                        for(auto& el : minimal_elements_tmp){
                            if(el.leq_poset_m(m_ji)){
                                is_minimal = false;
                                break;
                            }
                        }
                        if(is_minimal){
                            minimal_elements_tmp.push_back(m_ji);
                            grades.push(working_column.grade.join(grade_listsH[pivot][j]));
                        }
                    }
                    grade_listsH[pivot].push_back(working_column.grade);
                }
            }else{
                syzygiesH.push_back(SignatureColumn(working_column.grade, syzygiesH.size(), working_column.syzygy));
            }
            image_index++;
        }

        /* Reduce columns of map G */
        size_t& index_bound = index_map_high[v[v.size()-1]];
        for( size_t i=0; i<=index_bound; i++ ){
            column_index = -1;
            bool is_new=false;
            if(GB[i].size() > 0){
                if(pivot_mapH[i] == -1 || gb_columnsH[pivot_mapH[i]].last_updated != iter_index){
                    for(size_t j=GB[i].size()-1; j < GB[i].size(); j--){
                        if(GB[i][j].get_grade().leq_poset(v)){
                            column_index = (int)GB[i][j].get_index();
                            break;
                        }
                    }
                }
            } else{
                if(pivot_mapH[i] == -1 && grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        grade_t grade = working_column.grade;
                        if(is_new){
                            gb_columns[column_index] = working_column;
                            GB[i].push_back(signature_t(grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }else{
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }
                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        SignatureColumn syz(working_column.grade, GBH.size(), working_column.syzygy);
                        syz.syzygy.push(column_entry_t(1, GBH.size()));
                        GBH.push_back(std::vector<signature_t>());
                        SyzH.push_back(std::vector<signature_t>());
                        gb_columnsH.push_back(syz);
                        GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                        insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                        pivot_mapH[i] = gb_columnsH.size()-1;
                        gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_listsH[i].size(); j++){
                                    grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(working_column.grade.join(grade_listsH[i][j]));
                                    }
                                }
                                grade_listsH[i].push_back(working_column.grade);
                            }
                        }
                    }
                } else{
                    if(pivot != -1){
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                std::vector<grade_t> minimal_elements_tmp;
                                for(size_t j=0; j<grade_lists[pivot].size(); j++){
                                    grade_t m_ji = grade_lists[pivot][j].m_ji(gb_columns[column_index].grade);
                                    bool is_minimal = true;
                                    for(auto& el : minimal_elements_tmp){
                                        if(el.leq_poset_m(m_ji)){
                                            is_minimal = false;
                                            break;
                                        }
                                    }
                                    if(is_minimal){
                                        minimal_elements_tmp.push_back(m_ji);
                                        grades.push(gb_columns[column_index].grade.join(grade_lists[pivot][j]));
                                    }
                                }
                                grade_lists[pivot].push_back(gb_columns[column_index].grade);
                            }
                        }else{
                            pivot_map[pivot] = column_index;
                            gb_columns[column_index].last_updated = iter_index;
                        }
                    } else {
                        if(is_new){
                            SignatureColumn working_column = gb_columns[column_index];
                            SyzColumn syz_column = SyzColumn();
                            syz_column.push(column_entry_t(1, column_index));
                            SignatureColumn syz(working_column.grade, GBH.size(), syz_column);
                            syz.syzygy.push(column_entry_t(1, GBH.size()));
                            GBH.push_back(std::vector<signature_t>());
                            SyzH.push_back(std::vector<signature_t>());
                            gb_columnsH.push_back(syz);
                            GBH[GBH.size()-1].push_back(signature_t(syz.grade, gb_columnsH.size()-1));
                            insert_sorted(reorder_Z_map, signature_t(working_column.grade, GBH.size()-1));
                            pivot_mapH[i] = gb_columnsH.size()-1;
                            gb_columnsH[gb_columnsH.size()-1].last_updated = iter_index;
                            if(working_column.grade == v){
                                if(grade_listsH[i].size() == 0 || working_column.grade != grade_listsH[i].back()){
                                    std::vector<grade_t> minimal_elements_tmp;
                                    for(size_t j=0; j<grade_listsH[i].size(); j++){
                                        grade_t m_ji = grade_listsH[i][j].m_ji(working_column.grade);
                                        bool is_minimal = true;
                                        for(auto& el : minimal_elements_tmp){
                                            if(el.leq_poset_m(m_ji)){
                                                is_minimal = false;
                                                break;
                                            }
                                        }
                                        if(is_minimal){
                                            minimal_elements_tmp.push_back(m_ji);
                                            grades.push(working_column.grade.join(grade_listsH[i][j]));
                                        }
                                    }
                                    grade_listsH[i].push_back(working_column.grade);
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    Matrix filtered_gb_columnsH;
    size_t j = 0;
    for(size_t i=0; i<GBH.size(); i++){
        if(j<L_F.size() && i == L_F[j]){
            j++;
        } else {
            filtered_gb_columnsH.push_back(gb_columnsH[GBH[i][0].get_index()]);
        }
    }

    sort(filtered_gb_columnsH.begin(), filtered_gb_columnsH.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });

    hash_map<size_t, size_t> signature_map;
    for(size_t j=0; j<filtered_gb_columnsH.size(); j++){
        signature_map[filtered_gb_columnsH[j].signature_index] = j;
    }

    // Reindex syzygiesH
    Matrix reindexed_syzygies;
    for(size_t i=0; i<syzygiesH.size(); i++) {
        SignatureColumn sig_copy(syzygiesH[i]);
        SignatureColumn working_column(syzygiesH[i].grade, i);
        while(!sig_copy.empty()){
            column_entry_t entry = sig_copy.pop_pivot();
            working_column.push(column_entry_t(1, signature_map[entry.get_index()]));
        }
        reindexed_syzygies.push_back(working_column);
    }

    reindexed_syzygies = compute_minimal_generating_set(reindexed_syzygies);

    std::vector<grade_t> row_grades;
    for(auto& column : filtered_gb_columnsH) {
        row_grades.push_back(column.grade);
    }

    return std::pair<Matrix, std::vector<grade_t>>(reindexed_syzygies, row_grades);
}
