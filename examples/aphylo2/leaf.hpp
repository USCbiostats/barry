#include "../../include/barry/barry.hpp"

#ifndef LEAF_HPP
#define LEAF_HPP 1

// The same need to be locked
RULE_FUNCTION(rule_blocked) {
    if (Array->get_cell(i, j) == 9u)
        return false;
    return true;
};

using namespace phylocounters;

template<typename T1, typename T2>
std::vector< T1 > caster(const std::vector< T2 > & vec) {

    Vec< T1 > ans;
    ans.reserve(vec.size());

    for (auto &i : vec) {
        ans.push_back(*i);
    }

    return ans;

}

// Hasher
inline std::vector< double > tip_keygen(const PhyloArray & array) {
    
    // Baseline data: nrows and columns
    std::vector< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };
    
    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // Free cells
    for (auto i = 0u; i < array.nrow(); ++i)
        for (auto j = 0u; j < array.ncol(); ++j)
            dat.push_back((double) array.get_cell(i, j, false));
    
    return dat;
}

// Hasher
inline std::vector< double > tip_keygen_baseline(const PhyloArray & array) {
    
    // Baseline data: nrows and columns
    std::vector< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };
    
    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // // Free cells
    // for (auto i = 0u; i < array.nrow(); ++i)
    //     for (auto j = 0u; j < array.ncol(); ++j)
    //         dat.push_back((double) array.get_cell(i, j, false));
    
    return dat;
}

class Node {
public:
    unsigned int id;
    std::vector< unsigned int > dat = {};
    Node * parent                   = nullptr;
    std::vector< Node* > offspring  = {};

    Node(unsigned int id_) : id(id_) {};
    Node(unsigned int id_, std::vector< unsigned int > dat_) :
        id(id_), dat(dat_) {};
    ~Node() {};
};

class Leafs {
public:

    PhyloModel model_const;
    PhyloModel model_full;
    unsigned int nfuns;
    barry::Map< unsigned int, Node > dat;

    Leafs();

    Leafs(
        std::vector< std::vector<double> > & annotations,
        std::vector< unsigned int > & geneid,
        std::vector< unsigned int> & parent
        );

    ~Leafs();

    double operator()(std::vector< double > & par, unsigned int & i);

};

Leafs::Leafs() : model_const(), model_full(), dat() {
    return;
};

Leafs::Leafs(
    std::vector< std::vector<double> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< unsigned int> & parent
) : model_const(), model_full(), dat() {

    // Check the lengths
    if (annotations.size() == 0u)
        throw std::logic_error("Annotations is empty");

    nfuns = annotations.size();

    unsigned int n = annotations.at(0u).size();
    for (auto& iter : annotations) {
        if (iter.size() != n)
            throw std::length_error("Not all the annotations have the same length");
    }

    // Grouping up the data by parents
    for (unsigned int i = 0u; i < geneid.size(); ++i) {
        
        // Temp vector with the annotations
        std::vector< unsigned int > funs(nfuns);
        for (unsigned int j = 0u; j < nfuns; ++j)
            funs.at(j) = annotations.at(j).at(i);

        // Registered?
        auto iter = dat.find(parent.at(i));
        
        if (iter == dat.end()) {

            // Adding parent
            auto key_par = dat.insert({
                parent.at(i),
                Node({geneid.at(i)}) 
            });

            // Adding offspring
            auto key_off = dat.insert({
                geneid.at(i),
                Node({geneid.at(i), funs})
                });

            // Adding the offspring to the parent
            key_par.first->second.offspring.push_back(
                &key_off.first->second
            );

            // Adding the parent to the offspring
            key_off.first->second.parent = &key_par.first->second;

        } else { 
            // In this case, the parent exists, so we only need to assing the
            // offspring
            // Adding offspring
            auto key_off = dat.insert({
                geneid.at(i),
                Node({geneid.at(i), funs})
                });

            // Adding the offspring to the parent
            iter->second.offspring.push_back(
                &key_off.first->second
            );

            // Adding the parent to the offspring
            key_off.first->second.parent = &iter->second;
        }

    }




}

#endif