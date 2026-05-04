#ifndef MPH_RIPS_INCLUDED
#define MPH_RIPS_INCLUDED

#include <iostream>
#include <stdexcept>
#include <vector>
#include "utils.h"
#include "grade.h"
#include "matrix.h"
#include "signatureColumn.h"


//TODO: fix the virtual/override of the metrics and filters.
class Metric{
    /** Abstract class representing a metric. **/
public:
    virtual input_t eval(std::vector<input_t> a, std::vector<input_t> b) {return 0;}
};

// TODO: implement more types of metrcs.
class SquaredEuclideanMetric : public Metric{
public:
    virtual input_t eval(std::vector<input_t> a, std::vector<input_t> b) override{
        assert(a.size()==b.size());
        input_t ret = 0;
        for(size_t i=0; i<a.size(); i++){
            ret += (a[i]-b[i])*(a[i]-b[i]);
        }
        return ret;
    }
};


class Filter{
    /** Abstract class representing a filter function.
     A filter function is a function f: R^N -> R taking a point to its filtration value. **/
public:
    virtual input_t eval(std::vector<input_t> a) {return 0;}
};

// TODO: Implement more types of filters.
class XFilter : public Filter{
public:
    int axis;
    input_t coeff;
    XFilter(int _axis, input_t _coeff = 1) : Filter(){
        axis = _axis;
        coeff = _coeff;
    }
    virtual input_t eval(std::vector<input_t> a) override{
        return a[this->axis]*coeff;
    }
};

class DensityFilter : public Filter{
public:
    std::vector<std::vector<input_t>> points;
    hash_map<std::vector<input_t>*, input_t> cache;
    DensityFilter(std::vector<std::vector<input_t>>& _points) : Filter(){
        points = _points;
    }
    virtual input_t eval(std::vector<input_t> a) override{
        //if(cache.find(&a) != cache.end()){
        //    return cache[&a];
        //}
        std::vector<input_t> neighbourlist;
        SquaredEuclideanMetric m;
        for(auto& p:points){
            neighbourlist.push_back(m.eval(std::vector<input_t>(&a[0],&a[p.size()]), p));
        }
        std::sort(neighbourlist.begin(), neighbourlist.end());
        input_t density = 0;
        for(size_t i=0; i<neighbourlist.size(); i++){
            if(neighbourlist[i] > 0){
                input_t d = (i+1.0)/std::pow(neighbourlist[i], points[0].size());
                if(d>density){
                    density = d;
                }
            }
        }
        //cache[&a] = density;
        return density;
    }
};


template <typename T> class Simplex{
public:
    std::vector<size_t> vertices;
    std::vector<T> grade;
    size_t signature = -1;
    
    Simplex(std::vector<size_t> _vertices, std::vector<T> _grade){
        vertices = _vertices;
        grade = _grade;
    }
    
    Simplex(std::vector<size_t> _vertices, std::vector<T> _grade, size_t _signature){
        vertices = _vertices;
        grade = _grade;
        signature = _signature;
    }
};

class VertexContainer : public std::vector<size_t>{
public:
    VertexContainer() : std::vector<size_t>(){}
    VertexContainer(std::vector<size_t>& v) : std::vector<size_t>(v){}
    
    bool operator==(const VertexContainer& other) const {
        for(size_t index=0; index<this->size(); index++){
            if(this->at(index) != other[index]){
                return false;
            }
        }
        return true;
    }
};

struct VectorHasher
{
    /*size_t operator()(std::vector<size_t> const& vec) const {
        size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }*/
    std::size_t operator()(const VertexContainer& k) const
    {
        return boost::hash_range(k.begin(), k.end());
    }
};

Metric* parse_metric(int metric_index);
Filter* parse_filter(int filter_index);

void build_VR_subcomplex(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim);

void build_VR_subcomplex_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& point_grades, std::vector<size_t>& vertices, grade_t& prev_grade, grade_t& max_grade, std::vector<Simplex<index_t>>& low_simplices, std::vector<Simplex<index_t>>& mid_simplices, std::vector<Simplex<index_t>>& high_simplices, int hom_dim);

void build_VR_subcomplex_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<std::vector<input_t>>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim);

template <typename T> grade_t transform_grade(const std::vector<T>& input_grade, std::vector<hash_map<T, index_t>>& grade_map){
    /* Discretises a real-valued grade vector to its integer-index form by looking up each
     component in the per-axis grade_map (built by compute_grade_map). The output grade_t holds
     int64 indices into the per-axis sorted-distinct-values list.

     input_grade {std::vector<T>} -- the real-valued grade, one component per axis.
     grade_map {std::vector<hash_map<T, index_t>>} -- per-axis lookup from real value to its
        index in the sorted-distinct-values list for that axis.

     Returns:
     grade_t -- the discretised grade, with int64 indices in place of the real values.
     */
    grade_t grade;
    for (size_t i = 0; i < input_grade.size(); i++) {
        grade.push_back(grade_map[i][input_grade[i]]);
    }
    return grade;
}

template <typename T> void add_boundary_columns(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& low_simplices, Matrix& matrix, std::vector<hash_map<T, index_t>>& grade_map){
    /* Helper function to construct and add the columns of the boundary matrix from lists of simplices. */
    
    std::unordered_map<VertexContainer, size_t, VectorHasher> map;
    for(size_t simplex_index=0; simplex_index<low_simplices.size(); simplex_index++){
        map[VertexContainer(low_simplices[simplex_index].vertices)] = simplex_index;
    }
    
    for(size_t simplex_index=0; simplex_index < high_simplices.size(); simplex_index++){
        grade_t grade(transform_grade<T>(high_simplices[simplex_index].grade, grade_map));
        SignatureColumn column(grade, simplex_index);
        sort(high_simplices[simplex_index].vertices.begin(), high_simplices[simplex_index].vertices.end());
        int coeff = 1;

        for(size_t vertex_index = 0; vertex_index < high_simplices[simplex_index].vertices.size() ; vertex_index++){

            VertexContainer vertices;
            for(size_t k=0; k < high_simplices[simplex_index].vertices.size(); k++){
                if(k!=vertex_index){
                    vertices.push_back(high_simplices[simplex_index].vertices[k]);
                }
            }
            column.push(column_entry_t(coeff, map[vertices]));
            //coeff *= -1; //TODO: handle more types of coefficients than F_2.
        }
        matrix.push_back(column);
    }
}

template <typename T> void add_boundary_columns_multicritical(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& low_simplices, Matrix& matrix, std::vector<hash_map<T, index_t>>& grade_map){
    /* Helper function to construct and add the columns of the boundary matrix from lists of simplices.
     
     Assumption: low_simplices is ordered reverse lexicographically w.r.t each signature index
     
     */
    
    // Maps a collection of vertices to the lex-minimal simplex that has those vertices
    std::unordered_map<VertexContainer, size_t, VectorHasher> map; // we gotta rename this mfer
    for(size_t simplex_index=0; simplex_index<low_simplices.size(); simplex_index++){
        map[VertexContainer(low_simplices[simplex_index].vertices)] = simplex_index;
    }
    
    // Map d
    // For each high simplex
    for(size_t simplex_index=0; simplex_index < high_simplices.size(); simplex_index++){

        // because the grades are stored by indices, we need to get the actual grade value for it.
        grade_t grade(transform_grade<T>(high_simplices[simplex_index].grade, grade_map));

        SignatureColumn column(grade, simplex_index); // create object to store.
        int coeff = 1; // §kate says this will be 1 forever and we can ignore it, because we're working over F2

        for(size_t vertex_index = 0; vertex_index < high_simplices[simplex_index].vertices.size() ; vertex_index++){
            VertexContainer vertices;
            std::vector<size_t> high_vertices;
            for(size_t k=0; k < high_simplices[simplex_index].vertices.size(); k++){
                if(k!=vertex_index){
                    vertices.push_back(high_simplices[simplex_index].vertices[k]);
                    high_vertices.push_back(high_simplices[simplex_index].vertices[k]);
                }
            }

            std::vector<size_t> low_vertices = low_simplices[map[vertices]].vertices;
            sort(high_vertices.begin(), high_vertices.end());
            sort(low_vertices.begin(), low_vertices.end());
            // if( low_vertices != high_vertices){
                // std::cout << "Finished computing boundary matrices. add_boundary_columns_multicritical" << std::endl;
            // }
            size_t diff = 0;
            while ( true ) {
                bool should_add = true;
                for (size_t k=0; k<low_simplices[map[vertices]-diff].grade.size(); k++) {
                    if ( low_simplices[map[vertices]-diff].grade[k]> high_simplices[simplex_index].grade[k]) {
                        should_add = false;
                        break;
                    }
                }
                if (should_add) {
                    break;
                } else {
                    diff++;
                }
            }
            size_t m = map[vertices];
            low_vertices = low_simplices[map[vertices]-diff].vertices;
            if( low_vertices != high_vertices){
                throw std::runtime_error(" Inconsistent basis ");
            }
            column.push(column_entry_t(coeff, map[vertices]-diff));
            //coeff *= -1; //TODO: handle more types of coefficients than F_2.
        }
        matrix.push_back(column);
    }
    
    // Map pi
    std::vector<std::vector<Simplex<T>>> group_by_signature;
    for (size_t simplex_index=0; simplex_index < low_simplices.size(); simplex_index++) {
        size_t signature = low_simplices[simplex_index].signature;
        while (group_by_signature.size() <= signature) {
            group_by_signature.push_back(std::vector<Simplex<T>>());
        }
        bool should_add = true;
        for (size_t i=0; i < group_by_signature[signature].size(); i++) {
            bool is_lower = true;
            for (size_t j=0; j < group_by_signature[signature][i].grade.size(); j++) {
                if ( group_by_signature[signature][i].grade[j] > low_simplices[simplex_index].grade[j] ) {
                    is_lower = false;
                    break;
                }
            }
            if (is_lower) {
                should_add = false;
                break;
            }
        }
        if (should_add) {
            low_simplices[simplex_index].signature = simplex_index;
            group_by_signature[signature].push_back(low_simplices[simplex_index]);
        }
    }
    
    for (size_t signature=0; signature < group_by_signature.size(); signature++) {
        if (group_by_signature[signature].size() > 1) {
            for (size_t i=0; i < group_by_signature[signature].size(); i++) {
                for (size_t j=0; j < i; j++) {
                    std::vector<T> max_grade;
                    for (size_t t=0; t < group_by_signature[signature][i].grade.size(); t++) {
                        max_grade.push_back( group_by_signature[signature][i].grade[t] < group_by_signature[signature][j].grade[t] ? group_by_signature[signature][j].grade[t] : group_by_signature[signature][i].grade[t] );
                    }
                    grade_t grade(transform_grade<T>(max_grade, grade_map));
                    SignatureColumn column(grade, matrix.size());
                    column.push(column_entry_t(1, group_by_signature[signature][i].signature));
                    column.push(column_entry_t(1, group_by_signature[signature][j].signature));
                    matrix.push_back(column);
                }
            }
        }
    }
    
}


template <typename T> std::pair<std::vector<hash_map<T, index_t>>,std::vector<std::vector<T>>> compute_grade_map(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& mid_simplices, std::vector<Simplex<T>>& low_simplices){
    /* Computes the map needed to transform a 'input_t' type grading of the simplices to the 'grade_t' type used for computations. */

    std::vector<std::set<T>> index_values;
    size_t num_parameters = 0;

    if(low_simplices.size()>0){
        num_parameters = low_simplices[0].grade.size();
    }
    if(mid_simplices.size() > 0 && mid_simplices[0].grade.size()>num_parameters){
        num_parameters = mid_simplices[0].grade.size();
    }

    for(size_t parameter=0; parameter<num_parameters; parameter++){
        index_values.push_back(std::set<T>());
    }

    for (size_t simplex_index = 0; simplex_index < high_simplices.size(); simplex_index++){
        for(size_t parameter=0; parameter< high_simplices[simplex_index].grade.size(); parameter++){
            index_values[parameter].insert(high_simplices[simplex_index].grade[parameter]);
        }
    }

    for (size_t simplex_index = 0; simplex_index < mid_simplices.size(); simplex_index++){
        for(size_t parameter=0; parameter< mid_simplices[simplex_index].grade.size(); parameter++){
            index_values[parameter].insert(mid_simplices[simplex_index].grade[parameter]);
        }
    }

    for (size_t simplex_index = 0; simplex_index < low_simplices.size(); simplex_index++){
        for(size_t parameter=0; parameter< low_simplices[simplex_index].grade.size(); parameter++){
            index_values[parameter].insert(low_simplices[simplex_index].grade[parameter]);
        }
    }
    
    // convert sets to lists and sort them
    std::vector<std::vector<T>> index_value_lists;
    for(size_t parameter=0; parameter<num_parameters; parameter++){
        index_value_lists.push_back(std::vector<T>(index_values[parameter].begin(), index_values[parameter].end()));
        sort(index_value_lists[parameter].begin(), index_value_lists[parameter].end());
    }
    
    std::vector<hash_map<T, index_t>> grade_map;

    // construct a reverse lookup table, so given a crit value can get the index of it.
    for(size_t parameter=0; parameter<index_value_lists.size(); parameter++){
        grade_map.push_back(hash_map<T, index_t>()); // empty hash map

        for(size_t crit_val_ind=0; crit_val_ind<index_value_lists[parameter].size(); crit_val_ind++){
            grade_map[parameter] [index_value_lists[parameter][crit_val_ind]] = crit_val_ind;
            //         ^^ an integer indexing the paramters
            //                      ^^ keys in the hashmap are critical values
            //                                                                   ^^ values are indices
        }
    }
    return std::pair<std::vector<hash_map<T, index_t>>, std::vector<std::vector<T>>>(grade_map, index_value_lists);
}

std::pair<Matrix, Matrix> compute_boundary_matrices(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<input_t>& max_metric_values, int hom_dim);

std::vector<std::vector<input_t>> compute_distance_matrix(const std::vector<std::vector<input_t>>& point_cloud, Metric* metric);

std::pair<Matrix, Matrix> compute_boundary_matrices_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& grades, grade_t& max_grade, int hom_dim);

std::pair<Matrix, Matrix> compute_boundary_matrices_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim);

void verify_kernel(Matrix& kernel, Matrix& map);


std::vector<std::vector<input_t>> read_point_cloud(std::istream& input_stream, int precision);

template <typename SignatureColumn, typename Column> void read_input_file(std::istream& input_stream, std::vector<CustomSignatureColumn<SignatureColumn, Column>>& high_matrix, std::vector<CustomSignatureColumn<SignatureColumn, Column>>& low_matrix, int precision=1){
    /**
            Read FIRep formatted file (same as rivet).
     */
    std::string line;
    for(int i=0; i<3; i++){
        std::getline(input_stream, line);
    }
    std::getline(input_stream, line);
    std::istringstream s(line);
    std::vector<int> dims;
    for (int i; s >> i;) {
        dims.push_back(i);
        if (s.peek() == ',' || s.peek() == ' '){
            s.ignore();
        }
    }
    std::vector<std::vector<double>> high_grades;
    std::vector<std::vector<double>> index_value_lists;
    std::vector<hash_map<double, index_t>> visited_map;
    //int precision = 10000000;
    for(int d=0; d<dims[0]; d++){
        std::getline(input_stream, line);
        s = std::istringstream(line);
        bool reading_grade = true;
        CustomSignatureColumn<SignatureColumn, Column> column(grade_t(), d, dims[1], dims[0]);
        high_grades.push_back(std::vector<double>());
        column.syzygy.push(column_entry_t(1, d));
        for (double i; s >> i;) {
            if(reading_grade){
                /*high_grades[d].push_back(i);
                while(index_value_lists.size() < high_grades[d].size()){
                    index_value_lists.push_back(std::vector<double>());
                    visited_map.push_back(hash_map<double, index_t>());
                }
                index_value_lists[high_grades[d].size()-1].push_back(i);*/
                column.grade.push_back((int)(precision*i));
            }else{
                column.push(column_entry_t(1, (int)i));
            }
            while (s.peek() == ' '){
                s.ignore();
            }
            if(s.peek() == ';'){
                s.ignore();
                reading_grade = false;
                while (s.peek() == ' '){
                    s.ignore();
                }
            }
        }
        high_matrix.push_back(column);
    }
    std::vector<std::vector<double>> low_grades;
    for(int d=0; d<dims[1]; d++){
        std::getline(input_stream, line);
        s = std::istringstream(line);
        CustomSignatureColumn<SignatureColumn, Column> column(grade_t(), d, dims[2], dims[1]);
        low_grades.push_back(std::vector<double>());
        column.syzygy.push(column_entry_t(1, d));
        bool reading_grade = true;
        for (double i; s >> i;) {
            if(reading_grade){
                /*low_grades[d].push_back(i);
                while(index_value_lists.size() < low_grades[d].size()){
                    index_value_lists.push_back(std::vector<double>());
                }
                index_value_lists[low_grades[d].size()-1].push_back(i);*/
                column.grade.push_back((int)(precision*i));
            }else{
                column.push(column_entry_t(1, (int)i));
            }
            while (s.peek() == ' '){
                s.ignore();
            }
            if(s.peek() == ';'){
                s.ignore();
                reading_grade = false;
                while (s.peek() == ' '){
                    s.ignore();
                }
            }
        }
        low_matrix.push_back(column);
    }
    
    /*for(size_t i=0; i<index_value_lists.size(); i++){
        sort(index_value_lists[i].begin(), index_value_lists[i].end());
    }
    
    std::vector<hash_map<double, index_t>> grade_map;
    for(size_t i=0; i<index_value_lists.size(); i++){
        grade_map.push_back(hash_map<double, index_t>());
        for(size_t j=0; j<index_value_lists[i].size(); j++){
            grade_map[i][index_value_lists[i][j]] = j;
        }
    }
    for(size_t i=0; i<high_grades.size(); i++){
        for(size_t j=0; j<high_grades[i].size(); j++){
            high_matrix[i].grade.push_back(grade_map[i][high_grades[i][j]]);
        }
    }
    for(size_t i=0; i<low_grades.size(); i++){
        for(size_t j=0; j<low_grades[i].size(); j++){
            low_matrix[i].grade.push_back(grade_map[i][low_grades[i][j]]);
        }
    }*/
}

void apply_diagonal_grade_transform(Matrix& columns);

#endif // MPH_RIPS_INCLUDED
