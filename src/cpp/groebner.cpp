#include "groebner.h"

#include <algorithm>
#include <iostream>
#include <queue>

bool cmpID(const signature_t& lhs, const signature_t& rhs)
{
    return lhs.first.lt_colex(rhs.first);
}

void insert_sorted( std::vector<signature_t>& cont, signature_t value ) {
    std::vector<signature_t>::iterator it = std::upper_bound( cont.begin(), cont.end(), value, cmpID); // find proper position in descending order
    cont.insert( it, value ); // insert before iterator it
}

std::pair<Matrix, Matrix> computeGroebnerBases(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::vector<std::vector<index_t>> grade_base_set = get_grade_base_set<SignatureColumn>(columns);
    int n_total_grades = 1;
    for(auto& grade_list : grade_base_set){
        n_total_grades *= grade_list.size();
    }
    Iterator_lex grade_iterator = Iterator_lex(grade_base_set);

    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<SignatureColumn> syzygies;

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;

    for(size_t i=0; i<columns.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t max_pivot=0;
    for(auto& column : columns){
        if(max_pivot < column.get_pivot().get_index()){
            max_pivot = column.get_pivot().get_index();
        }
    }
    index_t pivot;

    int iter_index = 0;

    while(grade_iterator.has_next()){
        grade_t& v = grade_iterator.next();
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
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
                if(columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(SignatureColumn(working_column.grade, working_column.signature_index, working_column.syzygy));
                        Syz[i].push_back(signature_t(syzygies.back().grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }

    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    Matrix gb_columns_output;
    gb_columns_output.reserve(gb_columns.size());
    for(size_t i=0; i<GB.size(); i++){
        for(size_t j=0; j<GB[i].size(); j++){
            gb_columns_output.push_back(gb_columns[GB[i][j].get_index()]);
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns_output, syzygies_output);
}

std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt_min(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;

    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<VectorColumn> syzygies;
    gb_columns.reserve(2*columns.size());
    syzygies.reserve(columns.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());

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
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    index_t pivot;

    int iter_index = 0;

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
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
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
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
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

    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    Matrix gb_columns_output;
    gb_columns_output.reserve(gb_columns.size());
    for(size_t i=0; i<GB.size(); i++){
        for(size_t j=0; j<GB[i].size(); j++){
            gb_columns_output.push_back(gb_columns[GB[i][j].get_index()]);
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns_output, syzygies_output);
}

std::pair<Matrix, Matrix> computeGroebnerBases_gradeopt(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;

    /* Vectors to store the columns of the GBs */
    Matrix gb_columns;
    Matrix syzygies;
    gb_columns.reserve(columns.size());
    syzygies.reserve(columns.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());

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
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;

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
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
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
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = grade.join(grade_lists[pivot][ig]);
                                    if(grade != g){
                                        grades.push(g);
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
                        syzygies.push_back(SignatureColumn(working_column.get_grade(), syzygies.size(), working_column.syzygy));
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][ig]);
                                    if(gb_columns[column_index].grade != g){
                                        grades.push(g);
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

    // get_mem_usage(virt_memory, res_memory);

    std::cout << "Finished computing Groebner bases." << std::endl;
    return std::pair<Matrix, Matrix>(gb_columns, syzygies);
}

Matrix computekernel_gradeopt(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;

    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<SyzColumn> syzygies;
    gb_columns.reserve(2*columns.size());
    syzygies.reserve(columns.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());

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
    }

    index_t pivot;

    int iter_index = 0;

    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
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
                if(columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1){// && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            gb_columns.push_back(working_column);
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        }else{
                            GB[i][1].first = working_column.grade;
                            gb_columns[GB[i][1].get_index()].swap(working_column);
                            pivot_map[pivot] = GB[i][1].get_index();
                            gb_columns[GB[i][1].get_index()].last_updated = iter_index;
                        }
                        grade_t& grade = GB[i][1].get_grade();
                        if(grade == v){
                            if(grade_lists[pivot].size() == 0 || grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    if(grade.join(grade_lists[pivot][ig]) != grade){
                                        grade_t g = grade_lists[pivot][ig];
                                        grades.push(grade.join(grade_lists[pivot][ig]));
                                    }
                                }
                                grade_lists[pivot].push_back(grade);
                            }
                        }

                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            gb_columns.push_back(gb_columns[column_index]);
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, gb_columns.size()-1));
                            pivot_map[pivot] = gb_columns.size()-1;
                            gb_columns[gb_columns.size()-1].last_updated = iter_index;
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    if(gb_columns[column_index].grade.join(grade_lists[pivot][ig]) != gb_columns[column_index].grade){
                                        grades.push(gb_columns[column_index].grade.join(grade_lists[pivot][ig]));
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

    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return syzygies_output;
}

Matrix ImageGB(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;

    /* Vectors to store the columns of the GBs */
    Matrix gb_columns;
    gb_columns.reserve(columns.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());

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
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;

    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        /*while(v == grades.top()){
            grades.pop();
        }*/
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
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot_index();
                if(pivot_map[pivot] != -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                        if(working_column.grade == v){
                            if(grade_lists[pivot].size() == 0 || working_column.grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = working_column.grade.join(grade_lists[pivot][ig]);
                                    if(visited_grades.find(g) == visited_grades.end()){
                                        grades.push(g);
                                        visited_grades.insert(g);
                                    }
                                }
                                grade_lists[pivot].push_back(working_column.grade);
                            }
                        }

                    }else{
                        if(is_new){
                            gb_columns.pop_back();
                        }
                        working_column.syzygy.refresh();
                        Syz[i].push_back(signature_t(working_column.grade, -1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                            if(grade_lists[pivot].size() == 0 || gb_columns[column_index].grade != grade_lists[pivot].back()){
                                for(size_t ig=0; ig<grade_lists[pivot].size(); ig++){
                                    grade_t g = gb_columns[column_index].grade.join(grade_lists[pivot][ig]);
                                    if(visited_grades.find(g) == visited_grades.end()){
                                        grades.push(g);
                                        visited_grades.insert(g);
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
    std::cout << "Finished computing Groebner bases." << std::endl;
    return gb_columns;
}

std::pair<Matrix, std::vector<std::pair<grade_t, std::vector<size_t>>>> ImageGB_gradeopt_min_sig_basis(std::vector<SignatureColumn>& columns, MultigradedBasis ker_basis){
    /*
     The main function computing a Groebner basis for the image and kernel of the map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.
     ker_basis {MultigradedBasis} assumed to be sorted in lex.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
    }

    /* Compute index set iterator */
    std::priority_queue<grade_t, std::vector<grade_t>, std::greater<grade_t>> grades;
    std::vector<std::vector<grade_t>> grade_lists;
    std::unordered_set<grade_t, GradeHasher> visited_grades;

    /* Vectors to store the columns of the GBs */
    std::vector<SignatureColumn> gb_columns;
    std::vector<VectorColumn> syzygies;
    gb_columns.reserve(2*columns.size());
    syzygies.reserve(columns.size());

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB;
    std::vector<std::vector<signature_t>> Syz;
    GB.reserve(columns.size());
    Syz.reserve(columns.size());

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
        if(visited_grades.find(column.grade) == visited_grades.end()){
            grades.push(column.grade);
            visited_grades.insert(column.grade);
        }
    }
    for(auto& entry : ker_basis){
        if(visited_grades.find(entry.first) == visited_grades.end()){
            grades.push(entry.first);
            visited_grades.insert(entry.first);
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(columns.size());
    GradeHasher grade_hasher;
    for(auto& c : columns){
        grade_hashes.push_back(grade_hasher(c.grade));
    }
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    index_t pivot;

    int iter_index = 0;

    std::vector<std::pair<grade_t, std::vector<size_t>>> sig_basis;

    size_t ker_basis_index = 0;

    while(!grades.empty()){
        grade_t v = grades.top();
        grades.pop();
        while(v == grades.top()){
            grades.pop();
        }

        size_t grade_hash = grade_hasher(v);
        iter_index++;

        std::vector<size_t> sig_basis_at_v;

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
                if(grade_hash == grade_hashes[i] && columns[i].grade == v){
                    gb_columns.push_back(columns[i]);
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
                        sig_basis_at_v.push_back(pivot);
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
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
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        sig_basis_at_v.push_back(pivot);
                        pivot_map[pivot] = column_index;
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
        sig_basis.push_back(std::pair<grade_t, std::vector<size_t>>(v, sig_basis_at_v));
    }
    Matrix gb_columns_output;
    gb_columns_output.reserve(gb_columns.size());
    for(size_t i=0; i<GB.size(); i++){
        for(size_t j=0; j<GB[i].size(); j++){
            gb_columns_output.push_back(gb_columns[GB[i][j].get_index()]);
        }
    }
    std::cout << "Finished computing Groebner basis." << std::endl;
    return std::pair<Matrix, std::vector<std::pair<grade_t, std::vector<size_t>>>>(gb_columns_output, sig_basis);
}

Matrix buchberger(Matrix& columns){
    std::cout << "Computing GB for the image using Buchbergers algorithm." << std::endl;
    size_t max_pivot = 0;
    for(size_t i=0; i<columns.size(); i++){
        if(columns[i].get_pivot_index() > max_pivot){
            max_pivot = columns[i].get_pivot_index();
        }
    }
    std::vector<std::vector<size_t>> pivot_map(max_pivot+1, std::vector<size_t>());
    Matrix G = columns;
    std::queue<std::pair<size_t, size_t>> M;
    for(size_t i=0; i<G.size(); i++){
        for(size_t j=0; j<i; j++){
            if(G[i].get_pivot_index() == G[j].get_pivot_index()){
                M.push(std::pair<size_t, size_t>(i, j));
            }
        }
    }
    while(M.size()>0){
        std::pair<size_t, size_t> p = M.front();
        M.pop();
        SignatureColumn c = G[p.first];
        c.plus(G[p.second]);
        index_t pivot = c.get_pivot_index();
        while(pivot != -1){
            bool has_reduced = false;
            for(size_t i=0; i<pivot_map[pivot].size(); i++){
                if(G[pivot_map[pivot][i]].grade.leq_poset(c.grade)){
                    c.plus(G[pivot_map[pivot][i]]);
                    c.refresh();
                    c.syzygy.refresh();
                    has_reduced = true;
                    break;
                }
            }
            if(has_reduced){
                pivot = c.get_pivot_index();
            }else{

                for(size_t i=0; i<G.size(); i++){
                    if(G[i].get_pivot_index() == pivot){
                        M.push(std::pair<size_t, size_t>(i, G.size()));
                    }
                }
                G.push_back(c);
                pivot_map[pivot].push_back(G.size()-1);
                break;
            }
        }

    }
    std::cout << "Finished computing a GB of size: " << G.size() << std::endl;
    return G;
}


Matrix computeKernel_2p(std::vector<SignatureColumn>& columns){
    /*
     The main function computing a Groebner basis for the kernel of the bigraded map described by the list of columns 'columns'.

     Arguments:
     columns {std::vector<SignatureColumn>} -- columns describing the matrix of a map between two free multigraded momdules.

     Returns:
     std::vector<SignatureColumn> -- a list of vectors decribing a minimal Groebner basis for the image of the map.
     std::vector<SignatureColumn> -- a list of vectors describing a minimal Groebner basis for the kernel of the map.
     */

    std::cout << "Starting to compute Groebner bases..." << std::endl;

    /* Sort columns colexicographically */
    sort(columns.begin(), columns.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });
    // The sorted columns should agree with the columns sorted by index of signature
    hash_map<size_t, size_t> index_map_low, index_map_high;
    for(size_t i=0; i<columns.size(); i++){
        columns[i].signature_index = i;
        index_map_high[columns[i].grade[columns[i].grade.size()-1]] = i;
        if(index_map_low.find(columns[i].grade[columns[i].grade.size()-1]) == index_map_low.end()){
            index_map_low[columns[i].grade[columns[i].grade.size()-1]] = i;
        }
    }

    /* Compute index set iterator */
    std::vector<std::vector<index_t>> grade_base_set = get_grade_base_set<SignatureColumn>(columns);
    int n_total_grades = 1;
    for(auto& grade_list : grade_base_set){
        n_total_grades *= grade_list.size();
    }
    Iterator_lex grade_iterator = Iterator_lex(grade_base_set);

    std::vector<SyzColumn> syzygies;
    std::vector<std::vector<signature_t>> Syz;

    for(size_t i=0; i<columns.size(); i++){
        Syz.push_back(std::vector<signature_t>());
    }

    /* Main algorithm that iterates through the index set */
    index_t pivot;
    std::vector<index_t> pivot_map(columns.size(), -1);

    while(grade_iterator.has_next()){
        grade_t& v = grade_iterator.next();

        /* Initialize Macaulay matrix */
        for( size_t i=index_map_low[v[v.size()-1]]; i<=index_map_high[v[v.size()-1]] && columns[i].grade[0] <= v[0]; i++ ){
            if(columns[i].size()>0){
                pivot = columns[i].get_pivot_index();
                if(pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                    SignatureColumn& working_column = columns[i];
                    while(pivot != -1 && pivot_map[pivot] > -1 && pivot_map[pivot] < i){
                        working_column.plus(columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        pivot_map[pivot] = i;
                    }else{
                        working_column.syzygy.refresh();
                        syzygies.push_back(working_column.syzygy);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = i;
                    }
                }
            }
        }
    }

    Matrix syzygies_output;
    syzygies_output.reserve(syzygies.size());
    for(size_t i=0; i<Syz.size(); i++){
        for(size_t j=0; j<Syz[i].size(); j++){
            syzygies_output.push_back(SignatureColumn(Syz[i][j].get_grade(), syzygies_output.size(), syzygies[Syz[i][j].get_index()]));
        }
    }
    std::cout << "Finished computing Groebner bases." << std::endl;
    return syzygies_output;
}


Matrix compute_minimal_generating_set(Matrix& generators){
    /** Computes a minimal generating set for the module described by the columns in 'generators'.

     Arguments:
     generators {Matrix} -- a matrix whose columns a generators of a module.

     Returns:
     Matrix -- a subset of the columns of 'generators' that constitute a minimal set of generators for the module.
     */

    /* Vectors to store the columns of the GBs */
    std::cout << "Starting to compute minimal generating set for columns of size: " << generators.size() << std::endl;

    /* Sort columns colexicographically */
    sort(generators.begin(), generators.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });

    Matrix gb_columns;

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB; //The columns of 'generators' with a non-empty signature list is a generator.
    std::vector<std::vector<signature_t>> Syz;

    for(size_t i=0; i<generators.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
        generators[i].signature_index = i;
    }

    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<generators.size(); i++){
        generators[i].signature_index = i;
        index_map_high[generators[i].grade[generators[i].grade.size()-1]] = i;
    }

    std::vector<grade_t> index_list;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        index_list.push_back(generators[column_index].grade);
    }
    sort(index_list.begin(), index_list.end()); // Sort grades lexicographically


    index_t max_pivot=0;
    for(auto& generator : generators){
        if(max_pivot < generator.get_pivot_index()){
            max_pivot = generator.get_pivot_index();
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(generators.size());
    GradeHasher grade_hasher;
    for(auto& c : generators){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t pivot;
    std::vector<index_t> pivot_map(max_pivot+1, -1);
    int iter_index = 0;

    for(size_t index=0 ; index<index_list.size(); index++){
        if(index > 0 && index_list[index]==index_list[index-1]){
            continue;
        }
        grade_t v = index_list[index];
        iter_index++;
        size_t grade_hash = grade_hasher(v);


        /* Initialize Macaulay matrix */
        for( size_t i=0; i<=index_map_high[v[v.size()-1]]; i++ ){
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
                if(grade_hash == grade_hashes[i] && generators[i].grade == v){
                    gb_columns.push_back(generators[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot > -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                    SignatureColumn working_column = gb_columns[column_index];
                    while(pivot != -1 && pivot_map[pivot] > -1 && gb_columns[pivot_map[pivot]].last_updated == iter_index){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(is_new){
                        gb_columns.pop_back();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        Syz[i].push_back(signature_t(working_column.grade, -1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }
    Matrix minimal_generating_set;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        if(GB[column_index].size()>0){
            minimal_generating_set.push_back(generators[column_index]);
        }
    }
    //get_mem_usage(virt_memory, res_memory);
    return minimal_generating_set;
}

std::pair<Matrix, Matrix> compute_minimal_generating_set2(Matrix generators){
    /** Computes a minimal generating set for the module described by the columns in 'generators'.

     Arguments:
     generators {Matrix} -- a matrix whose columns a generators of a module.

     Returns:
     Matrix -- a subset of the columns of 'generators' that constitute a minimal set of generators for the module.
     */

    /* Vectors to store the columns of the GBs */
    std::cout << "Starting to compute minimal generating set for columns of size: " << generators.size()  << std::endl;

    for(size_t i=0; i<generators.size(); i++){
        generators[i].signature_index = i;
    }

    /* Sort columns colexicographically */
    sort(generators.begin(), generators.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.grade.lt_colex(rhs.grade);
         });

    hash_map<size_t, size_t> reorder_map;
    for(size_t i=0; i<generators.size(); i++){
        reorder_map[i] = generators[i].signature_index;
    }

    Matrix gb_columns;
    std::vector<SyzColumn> syzygies;

    /* Index maps to keep track of signatures */
    std::vector<std::vector<signature_t>> GB; //The columns of 'generators' with a non-empty signature list is a generator.
    std::vector<std::vector<signature_t>> Syz;

    for(size_t i=0; i<generators.size(); i++){
        GB.push_back(std::vector<signature_t>());
        Syz.push_back(std::vector<signature_t>());
        generators[i].signature_index = i;
        generators[i].syzygy = SyzColumn();
        generators[i].syzygy.push(column_entry_t(1, i));
    }

    hash_map<size_t, size_t> index_map_high;
    for(size_t i=0; i<generators.size(); i++){
        index_map_high[generators[i].grade[generators[i].grade.size()-1]] = i;
    }

    std::vector<grade_t> index_list;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        index_list.push_back(generators[column_index].grade);
    }
    sort(index_list.begin(), index_list.end()); // Sort grades lexicographically

    index_t max_pivot=0;
    for(auto& generator : generators){
        if(max_pivot < generator.get_pivot().get_index()){
            max_pivot = generator.get_pivot().get_index();
        }
    }

    std::vector<size_t> grade_hashes;
    grade_hashes.reserve(generators.size());
    GradeHasher grade_hasher;
    for(auto& c : generators){
        grade_hashes.push_back(grade_hasher(c.grade));
    }

    /* Main algorithm that iterates through the index set */
    int column_index;
    index_t pivot;

    int iter_index = 0;

    for(size_t index=0 ; index<index_list.size(); index++){
       // if(index > 0 && index_list[index]==index_list[index-1]){
       //     continue;
       // }
        grade_t& v = index_list[index];
        iter_index++;
        std::vector<index_t> pivot_map(max_pivot+1, -1);
        size_t grade_hash = grade_hasher(v);
        /* Initialize Macaulay matrix */
        for( size_t i=0; i<generators.size(); i++ ){
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
                if(grade_hash == grade_hashes[i] && generators[i].grade == v){
                    gb_columns.push_back(generators[i]);
                    column_index = (int)gb_columns.size()-1;
                    is_new=true;
                }
            }
            if(column_index > -1){
                pivot = gb_columns[column_index].get_pivot().get_index();
                if(pivot > -1 && pivot_map[pivot] > -1){
                    SignatureColumn working_column(gb_columns[column_index]);
                    SyzColumn syz;
                    while(pivot != -1 && pivot_map[pivot] > -1){
                        working_column.plus(gb_columns[pivot_map[pivot]]);
                        syz.plus(gb_columns[pivot_map[pivot]].syzygy);
                        pivot = working_column.get_pivot().get_index();
                    }
                    if(pivot != -1){
                        working_column.refresh();
                        working_column.syzygy.refresh();
                        gb_columns.push_back(working_column);
                        GB[i].push_back(signature_t(working_column.grade, gb_columns.size()-1));
                        pivot_map[pivot] = gb_columns.size()-1;
                        gb_columns[gb_columns.size()-1].last_updated = iter_index;
                    }else{
                        working_column.syzygy.refresh();
                        syz.refresh();
                        syzygies.push_back(syz);
                        Syz[i].push_back(signature_t(working_column.grade, syzygies.size()-1));
                    }
                } else{
                    if(pivot != -1){
                        pivot_map[pivot] = column_index;
                        gb_columns[column_index].last_updated = iter_index;
                        if(is_new){
                            GB[i].push_back(signature_t(gb_columns[column_index].grade, column_index));
                        }
                    }
                }
            }
        }
    }
    Matrix minimal_generating_set;
    Matrix change_of_basis_map;
    hash_map<size_t, size_t> reindex_map;
    for(size_t column_index=0 ; column_index<generators.size(); column_index++){
        if(GB[column_index].size()>0){
            reindex_map[column_index] = minimal_generating_set.size();
            minimal_generating_set.push_back(generators[column_index]);
            SignatureColumn column(generators[column_index].grade, reorder_map[column_index]);
            column.push(column_entry_t(1, reindex_map[column_index]));
            change_of_basis_map.push_back(column);
        }else{
            SignatureColumn column(generators[column_index].grade, reorder_map[column_index]);
            index_t pivot = syzygies[Syz[column_index][0].get_index()].get_pivot_index();
            while(pivot != -1){
                column.push(column_entry_t(1, reindex_map[pivot]));
                syzygies[Syz[column_index][0].get_index()].pop_pivot();
                pivot = syzygies[Syz[column_index][0].get_index()].get_pivot_index();
            }
            change_of_basis_map.push_back(column);
        }
    }

    sort(change_of_basis_map.begin(), change_of_basis_map.end(), [ ](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.signature_index < rhs.signature_index;
         });
    // get_mem_usage(virt_memory, res_memory);
    std::cout << "Finished computing minimal generating set of size: " << minimal_generating_set.size();
    return std::pair<Matrix, Matrix>(minimal_generating_set, change_of_basis_map);
}

Matrix compute_syzygy_module(Matrix groebner_basis){
    /** Computes the syzygy module of a Groebner basis using Schreyer's algorithm.

    */

    for(size_t i=0; i<groebner_basis.size(); i++){
        groebner_basis[i].signature_index = i;
    }

    hash_map<size_t, std::vector<SignatureColumn>> pivot_partition;
    for(size_t i=0; i<groebner_basis.size(); i++){
        if(pivot_partition.find(groebner_basis[i].get_pivot().get_index()) == pivot_partition.end()){
            pivot_partition[groebner_basis[i].get_pivot().get_index()] = std::vector<SignatureColumn>();
        }
        pivot_partition[groebner_basis[i].get_pivot().get_index()].push_back(groebner_basis[i]);
    }
    Matrix syzygies;
    for(auto& entry : pivot_partition){
        std::vector<SignatureColumn>& columns = entry.second;
        for(size_t i=1; i<columns.size(); i++){
            std::vector<grade_t> minimal_elements_tmp;
            std::vector<size_t> minimal_indices_tmp;
            for(size_t j=0; j<i; j++){
                grade_t m_ji = columns[j].grade.m_ji(columns[i].grade);
                bool is_minimal = true;
                for(auto& el : minimal_elements_tmp){
                    if(el.leq_poset(m_ji)){
                        is_minimal = false;
                        break;
                    }
                }
                if(is_minimal){
                    minimal_elements_tmp.push_back(m_ji);
                    minimal_indices_tmp.push_back(j);
                }
            }
            std::vector<grade_t> minimal_elements;
            std::vector<size_t> minimal_indices;
            for(size_t j=0; j<minimal_elements_tmp.size(); j++){
                bool is_minimal = true;
                for(size_t k=0; k<minimal_elements_tmp.size(); k++){
                    if(k!=j && minimal_elements_tmp[k].leq_poset(minimal_elements_tmp[j])){
                        is_minimal = false;
                        break;
                    }
                }
                if(is_minimal){
                    minimal_elements.push_back(minimal_elements_tmp[j]);
                    minimal_indices.push_back(minimal_indices_tmp[j]);
                }
            }
            /*spdlog::info("Minimal elements...");
            for(auto& min_el : minimal_elements){
                min_el.print();
            }*/
            for(size_t j=0; j<minimal_elements.size(); j++){
                SignatureColumn column_to_reduce(columns[i]);
                column_to_reduce.plus(columns[minimal_indices[j]]);
                SignatureColumn syzygy(column_to_reduce.grade, -1);
                index_t pivot = column_to_reduce.get_pivot().get_index();
                while(pivot != -1){
                    if(pivot_partition.find(pivot) == pivot_partition.end()){
                        throw std::runtime_error("Unkown pivot.");
                    }
                    bool has_reduced = false;
                    for(auto& column : pivot_partition[pivot]){
                        if(column.grade.leq_poset(column_to_reduce.grade)){
                            column_to_reduce.plus(column);
                            syzygy.push(column_entry_t(1, column.signature_index));
                            has_reduced = true;
                            break;
                        }
                    }
                    if(!has_reduced){
                        throw std::runtime_error("Could not reduce column.");
                    }
                    pivot = column_to_reduce.get_pivot().get_index();
                }
                syzygy.push(column_entry_t(1, columns[i].signature_index));
                syzygy.push(column_entry_t(1, columns[minimal_indices[j]].signature_index));
                if(syzygy.get_pivot().get_index() != -1){
                    syzygies.push_back(syzygy);
                }

            }
        }
    }
    return syzygies;
}
