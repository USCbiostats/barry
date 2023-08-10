#ifndef GEESE_TYPES_HPP
#define GEESE_TYPES_HPP
/**
 * @name Convenient typedefs for Node objects.
 * */
/**
 * @brief Data definition for the `PhyloArray` class.
 * 
 * This holds basic information about a given node.
 * 
 */
class NodeData {
public:
  
    /**
     * Branch length.
     */
    std::vector< double > blengths = {};
    
    /**
     * State of the parent node.
     */
    std::vector< bool > states = {};
    
    /**
     * 
     */
    bool duplication = true;
    
    // NodeData() : blengths(0u), states(0u) {};
    
    NodeData(
        const std::vector< double > & blengths_,
        const std::vector< bool > & states_,
        bool duplication_ = true
    ) : blengths(blengths_), states(states_), duplication(duplication_) {};
    
    // ~NodeData() {};
  
};

class PhyloCounterData {
private:
    std::vector< size_t > data;
    std::vector< double > * counters;

public:
    PhyloCounterData(
        std::vector< size_t > data_,
        std::vector< double > * counters_ = nullptr
        ) : data(data_), counters(counters_) {};

    PhyloCounterData() : data(0u) {};

    size_t at(size_t d) {return data.at(d);};
    size_t operator()(size_t d) {return data.at(d);};
    size_t operator[](size_t d) {return data[d];};
    void reserve(size_t x) {return data.reserve(x);};
    void push_back(size_t x) {return data.push_back(x);};
    void shrink_to_fit()  {return data.shrink_to_fit();};
    size_t size() {return data.size();};

    std::vector< size_t >::iterator begin() {return data.begin();};
    std::vector< size_t >::iterator end() {return data.end();};

    bool empty() {return data.empty();};
    std::vector< double > * get_counters() {return counters;};

};

class PhyloRuleDynData {
public:
    const std::vector< double > * counts;
    size_t pos;
    size_t lb;
    size_t ub;
    size_t duplication;

    PhyloRuleDynData(
        const std::vector< double > * counts_,
        size_t pos_,
        size_t lb_,
        size_t ub_,
        size_t duplication_
        ) :
        counts(counts_), pos(pos_), lb(lb_), ub(ub_), duplication(duplication_) {};

    const double operator()() const
    {
        return (*counts)[pos];
    }
    
    ~PhyloRuleDynData() {};
    
};


typedef std::vector< std::pair< size_t, size_t > > PhyloRuleData;

///@{
typedef barry::BArrayDense<size_t, NodeData> PhyloArray;
typedef barry::Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef barry::Counters< PhyloArray, PhyloCounterData> PhyloCounters;

typedef barry::Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef barry::Rules<PhyloArray,PhyloRuleData> PhyloRules;

typedef barry::Rule<PhyloArray,PhyloRuleDynData> PhyloRuleDyn;
typedef barry::Rules<PhyloArray,PhyloRuleDynData> PhyloRulesDyn;

typedef barry::Support<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloSupport;
typedef barry::StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef barry::Model<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloModel;
typedef barry::PowerSet<PhyloArray, PhyloRuleData> PhyloPowerSet;
///@}

#endif