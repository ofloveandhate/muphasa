#include "presentation.h"

#include <iostream>
#include <queue>
#include <algorithm>


bool cmpID(const signature_t& lhs, const signature_t& rhs)
{
    return lhs.first.lt_colex(rhs.first);
}

void insert_sorted( std::vector<signature_t>& cont, signature_t value ) {
    std::vector<signature_t>::iterator it = std::upper_bound( cont.begin(), cont.end(), value, cmpID); // find proper position in descending order
    cont.insert( it, value ); // insert before iterator it
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

    std::cout << "Starting to presentation degree by degree..."  << std::endl;

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
    std::cout << "Nonminimal presentation of size "<< syzygiesH.size() << std::endl;

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

//get_mem_usage(virt_memory, res_memory);
    return std::pair<Matrix, std::vector<grade_t>>(reindexed_syzygies, row_grades);
}
