#include "matrix.h"


void reduce_columns(std::vector<SignatureColumn>& columns, std::vector<index_t>& indices, std::vector<SignatureColumn>& modified_columns, bool control_size) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    hash_map<index_t, index_t> pivot_map;
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(pivot_map.find(columns[indices[i]].get_pivot().get_index()) != pivot_map.end()){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            SignatureColumn working_column(columns[indices[i]].grade, columns[indices[i]].signature_index);
            working_column.plus(columns[indices[i]]);
            size_t input_size=working_column.size();

            while(!working_column.empty() && pivot_map.find(pivot) != pivot_map.end()){
                working_column.plus(columns[pivot_map[pivot]]);
                pivot = working_column.get_pivot().get_index();
                if(control_size && working_column.size() > 2*input_size){
                    working_column.refresh();
                    working_column.syzygy.refresh();
                    input_size = working_column.size();
                }
            }
            
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                pivot_map[columns[indices[i]].get_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}


void reduce_columns_custom_hash(Matrix& columns, std::vector<index_t>& indices, Matrix& modified_columns) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    CustomHashMap map(indices.size());
    
    
    for(size_t i=0; i<indices.size(); i++){
        if(map.containsKey(columns[indices[i]].get_pivot().get_index())){
            index_t pivot = columns[indices[i]].get_pivot().get_index();
            SignatureColumn working_column(columns[indices[i]].grade, columns[indices[i]].signature_index);
            working_column.plus(columns[indices[i]]);
            
            
            while(!working_column.empty() && map.containsKey(pivot)){
                working_column.plus(columns[map.get(pivot)]);
                pivot = working_column.get_pivot().get_index();
            }
            if(!working_column.empty()){
                map.insert(working_column.get_pivot().get_index(), indices[i]);
            }
            
            modified_columns.push_back(working_column);
        } else{
            if(!columns[indices[i]].empty()){
                map.insert(columns[indices[i]].get_pivot().get_index(),  indices[i]);
            }
        }
    }
    return;
}




void reduce_columns_par(Matrix& columns, Matrix& additional_columns, std::vector<index_t>& indices, std::vector<SignatureColumn*>& modified_columns) {
    /**  The function reduces a subset of the columns 'columns' described by the indices 'indices'. The function returns a pair containing a list of reduced columns and their syzygy-columns describing the operations performed on the columns. Only columns that were modified by the function will be returned in this list.
     
     Arguments:
     columns {std::vector<SignatureColumn>} -- a list of columns.
     indices {std::vector<index_t>} -- a list of indices that specify a submatrix from the list 'columns' that should be reduced.
     
     Returns:
     std::vector<SignatureColumn> -- a list of columns that were modified by the function in reduced form.
     */
    
    hash_map<index_t, index_t> pivot_map;
    
    
    
    for(size_t i=0; i<indices.size(); i++){
        bool base_list = indices[i] < columns.size();
        if(pivot_map.find((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index()) != pivot_map.end()){
            SignatureColumn working_column((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).grade, (base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).signature_index);
            working_column.plus((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]));
            index_t pivot = working_column.get_pivot().get_index();
            
            while(!working_column.empty() && (pivot_map.find(pivot) != pivot_map.end())){
                working_column.plus(pivot_map[pivot] < columns.size() ? columns[pivot_map[pivot]] : additional_columns[pivot_map[pivot]-columns.size()]);
                pivot = working_column.get_pivot().get_index();
            }
            if(!working_column.empty()){
                pivot_map[working_column.get_pivot().get_index()] = indices[i];
            }
            modified_columns.push_back(new SignatureColumn(working_column));
        } else{
            if((base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index() != -1){
                pivot_map[(base_list ? columns[indices[i]] : additional_columns[indices[i]-columns.size()]).get_cached_pivot().get_index()] = indices[i];
            }
        }
    }
    return;
}




Matrix inverse(Matrix& M) {
    Matrix M_inv;
    for ( size_t i=0; i<M.size(); i++ ) {
        M[i].syzygy = SyzColumn();
        M[i].syzygy.push(column_entry_t(1, i));
    }
    index_t pivot;
    std::vector<index_t> pivot_map(M.size()+1, -1);
    for( size_t i=0; i<M.size(); i++ ){
        SignatureColumn working_column = M[i];
        pivot = working_column.get_pivot().get_index();
        while(pivot != -1 && pivot_map[pivot] > -1){
            working_column.plus(M[pivot_map[pivot]]);
            pivot = working_column.get_pivot().get_index();
        }
        if(pivot != -1){
            pivot_map[pivot] = i;
        }
    }
    for( size_t i=0; i<M.size(); i++ ){
        SignatureColumn working_column = M[i];
        working_column.pop_pivot();
        pivot = working_column.get_pivot().get_index();
        while(pivot != -1 && pivot_map[pivot] > -1){
            working_column.plus(M[pivot_map[pivot]]);
            pivot = working_column.get_pivot().get_index();
        }
        if (pivot != -1) {
            throw std::runtime_error("Failed to compute inverse");
        }
        SignatureColumn inv_col(M[i].grade, M[i].get_pivot().get_index(), working_column.syzygy);
        M_inv.push_back(inv_col);
    }
    sort(M_inv.begin(), M_inv.end(), [](  SignatureColumn& lhs,  SignatureColumn& rhs )
         {
             return lhs.signature_index < rhs.signature_index ;
         });
    return M_inv;
}



Matrix mat_v_cut(Matrix& M, index_t index) {
    Matrix Mv;
    for ( size_t i=0; i<M.size(); i++ ) {
        SignatureColumn working_column = M[i];
        index_t pivot = working_column.get_pivot().get_index();
        while( pivot >= index ) {
            working_column.pop_pivot();
            pivot = working_column.get_pivot().get_index();
        }
        Mv.push_back(working_column);
    }
    return Mv;
}



Matrix matmul(Matrix& A, Matrix& B) {
    Matrix M;
    for ( size_t i=0; i<B.size(); i++ ) {
        SignatureColumn working_column = B[i];
        SignatureColumn ret(B[i].grade, i);
        while( !working_column.empty() ) {
            index_t pivot = working_column.get_pivot().get_index();
            ret.plus(A[pivot]);
            working_column.pop_pivot();
        }
        M.push_back(ret);
    }
    return M;
}

