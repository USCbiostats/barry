#include "../../include/barry/barry.hpp"

#ifndef LEAF_HPP
#define LEAF_HPP 1

// The same need to be locked
RULE_FUNCTION(rule_blocked) {
    if (Array->get_cell(i, j) == 9u)
        return false;
    return true;
}

using namespace phylocounters;

template<typename T1, typename T2>
std::vector< T1 > caster(const std::vector< T2 > & vec) {

    std::vector< T1 > ans;
    ans.reserve(vec.size());

    for (auto &i : vec) {
        ans.push_back(*i);
    }

    return ans;

}

// Hasher
inline std::vector< double > keygen_full(const PhyloArray & array) {
    
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
inline std::vector< double > keygen_const(const PhyloArray & array) {
    
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

    int get_parent() const {
        if (parent == nullptr)
            return -1;
        else
            return (int) parent->id;
    };

    bool is_leaf() const {
        return offspring.size() == 0u;
    };

};

typedef std::vector< PhyloArray > InnerNode_type;
typedef std::vector< InnerNode_type > InnerNodes_type ;

class Leafs {
public:

    PhyloModel model_const;
    PhyloModel model_full;
    unsigned int nfuns;
    barry::Map< unsigned int, Node > dat;
    InnerNodes_type InnerNodes;
    InnerNode_type rawdata;
    std::vector< unsigned int > idx_const;
    std::vector< unsigned int > idx_full;

    Leafs();

    Leafs(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > & geneid,
        std::vector< unsigned int> & parent,
        PhyloCounters & counters
        );

    ~Leafs() {};



    double operator()(std::vector< double > & par, unsigned int & i);
    void print() const;

};

Leafs::Leafs() : model_const(), model_full(), dat() {
    return;
}

Leafs::Leafs(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< unsigned int> & parent,
    PhyloCounters & counters
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

    // Setting up the rules for the model
    model_const.set_keygen(keygen_const);
    model_full.set_keygen(keygen_full);

    model_const.add_rule(rule_blocked<PhyloArray,PhyloRuleData>);

    model_const.set_counters(&counters);
    model_full.set_counters(&counters);

    // Grouping up the data by parents -----------------------------------------
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

    // Generating the model data -----------------------------------------------

    // All combinations of the function
    PhyloPowerSet pset(nfuns, 1u);
    pset.calc();

    std::vector< std::vector< bool > > states;
    states.reserve(pset.data.size());
    unsigned int i = 0u;
    for (auto& iter : pset.data) {

        states.push_back(std::vector< bool >(nfuns, false));

        for (auto iter2 = iter.get_col(0u)->begin(); iter2 != iter.get_col(0u)->end(); ++iter2)
            states.at(i).at(iter2->first) = true;
        
        i++;
    }

    unsigned int narrays = 0u;

    // Iterating throught the nodes
    for (auto& iter : dat) {

        // Only parents get a node
        if (!iter.second.is_leaf()) {

            // Creating the phyloarray, nfuns x noffspring
            rawdata.push_back(PhyloArray(nfuns, iter.second.offspring.size()));

            // Adding the data, first through functions
            for (unsigned int k = 0u; k < nfuns; ++k) {

                // Then through the offspring
                unsigned int j = 0;
                for (auto& o : iter.second.offspring) {
                    
                    if (annotations.at(k).at(o->id) != 0) {
                        rawdata.at(narrays).insert_cell(
                            k, j, annotations.at(k).at(o->id), false, false
                            );
                    }
                    
                    j++;

                }

            }

            // We then need to set the powerset            
            InnerNodes.push_back(InnerNode_type(pset.data.size()));
            unsigned int i = 0u;
            std::vector< double > blen(iter.second.offspring.size(), 1.0);
            for (auto& iter : states) {
                
                InnerNodes.at(narrays).push_back(rawdata.at(narrays));
                InnerNodes.at(narrays).at(i).set_data(
                    new NodeData(blen, iter),
                    true
                );

                // Once the array is ready, we can add it to the model
                this->idx_const.push_back(
                    model_const.add_array(InnerNodes.at(narrays).at(i))
                    );
                this->idx_full.push_back(
                    model_full.add_array(InnerNodes.at(narrays).at(i))
                    );

            }

            narrays++;

        }
    }

}

void Leafs::print() const {
    for (auto& iter : this->dat) {
        std::cout << " Gene id: " << iter.first;
        std::cout << " Parent id: " << iter.second.get_parent();
        std::cout << " offsprings: {";
        for (auto& o : iter.second.offspring) {
            std::cout << o->id << ", ";
        }
        std::cout << "}" << std::endl;
    }
    return;
}

#endif