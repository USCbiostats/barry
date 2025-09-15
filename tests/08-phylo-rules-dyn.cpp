#include "../include/barry/barry.hpp"
#include "../include/barry/models/geese.hpp"

using namespace geese;

// Rule definition --Notice this has to be initialized with the
// counter itself.
class RuleDynD {
public:
    const std::vector< double > * counts;
    size_t pos;
    size_t max_changes;
    RuleDynD(const std::vector< double > * counts_, size_t pos_, size_t max_changes_) :
        counts(counts_), pos(pos_), max_changes(max_changes_) {};
    
    ~RuleDynD() {};
};

bool check_max_gains(const PhyloArray & A, size_t i, size_t j, RuleDynD * d)
{

    // The count must be included iff the number of changes
    // is less than the max allowed
    return d->max_changes > (d->counts->operator[](d->pos));

}

#include "tests.hpp"

BARRY_TEST_CASE("Phylo dynamic rules", "[phylo-dyn-rules]") {

    

    PhyloArray d(3, 2);
    d.set_data(new NodeData({1.0, 1.0, 1.0}, {true, false, false}), true);

    std::cout << "---- Support with reduced pset ---" << std::endl;
    // Generating the support
    barry::Support<PhyloArray,PhyloCounterData,PhyloRuleData,PhyloRuleDynData> S(d);
    counter_overall_changes(S.get_counters());
    
    // Creating a rule, we start with the data
    RuleDynD rd(S.get_current_stats(), 0, 2u);
    
    // barry::Rule<PhyloArray,RuleDynD> reduce_gains(check_max_gains, &rd, false); 
    // S.add_rule_dyn(reduce_gains);
    
    rule_dyn_limit_changes(&S, 0, 0, 2);
    S.calc();
    S.print();    

    std::cout << "---- Support with full pset ------" << std::endl;
    barry::Support<PhyloArray,PhyloCounterData,PhyloRuleData,RuleDynD> S2(d);
    counter_overall_changes(S2.get_counters());
    S2.calc();
    S2.print();

    // Computing differences
    size_t matches = 0u;

    const auto& s = S.get_data().get_index();
    auto D1 = S.get_data().get_data();
    auto D2 = S2.get_data().get_data();

    for (const auto& s2 : S2.get_data().get_index())
    {

        // Can we find it?
        auto iter = s.find(s2.first);

        if (iter == s.end()) 
            continue;

        // If we found it, is this matching?
        if (D1[(*iter).second] == D2[s2.second])
            ++matches;

    }

    #ifdef CATCH_CONFIG_MAIN
        REQUIRE(matches == S.get_data().size());
    #endif

}
