#include "exceptions.h"
#include "ordering.h"
#include "options.h"
#include <stdio.h>
#include <bn-to-cnf/bayesnet.h>
#include <bn-to-cnf/cnf.h>
#include <set>
#include <queue>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include "branchnbound.h"
#include "xary.h"
#include <cassert>
#include "bound.h"
#include "types.h"
#include "timer.h"

using namespace std;
using namespace bnc;

ordering_t::ordering_t(){
}

ordering_t& ordering_t::operator=(const ordering_t &o){
    std::vector<literal_t>::operator=(o);
    return *this;
}

ordering_t& ordering_t::operator=(ordering_t &o){
    std::vector<literal_t>::operator=(o);
    return *this;
}

void ordering_t::write(std::string filename, bayesnet *bn) const {
    FILE *file = fopen(filename.c_str(), "w");
    if(file){
        if(bn != NULL){
            string_ordering_t elim_ordering = get_variable_string_ordering(bn);
            for(auto it = elim_ordering.rbegin(); it != elim_ordering.rend(); it++)
                fprintf(file, "%s ", it->c_str());
        } else {
            for(auto it = begin(); it != end(); it++)
                fprintf(file, "%d ", *it);
        }
        fclose(file);
    } else throw compiler_ordering_exception("Could not write to ordering file '%s'", filename.c_str());
}


void ordering_t::read(std::string filename, bayesnet *bn){
    clear();

    FILE *file = fopen(filename.c_str(),"r");
    if(file){
        std::set<int> vars;
        int v;
        if(bn == NULL){
            while(fscanf(file,"%d",&v) != EOF){
                if(vars.find(v) != vars.end()){
                    fprintf(stderr, "Variable '%u' appears more than once in the ordering\n",v);
                    exit(1);
                }
                vars.insert(v);
                push_back(v);
            }
        } else {
            char word[100];
            while (fscanf(file,"%s",word) != EOF){
                v = bn->get_node_id(word);
                if(v >= 0){
                    if(vars.find(v) != vars.end()){
                        fprintf(stderr, "Variable '%s' appears more than once in the ordering\n",bn->get_node_name(v).c_str());
                        exit(1);
                    }
                    vars.insert(v);
                    push_back(v);
                } else {
                    fprintf(stderr, "Unknown variable in ordering: '%s'\n", word);
                    exit(1);
                }
            }
            std::reverse(begin(),end());
        }
        fclose(file);
    } else throw compiler_write_exception("Could not open ordering file '%s'", filename.c_str());
}

void all_topological_orders(bayesnet *bn, std::vector<unsigned int> &kVariables, std::vector<bool> &visited, std::vector<unsigned int> &indegree, const Bound &kBound, ordering_t &ordering, ordering_t &best_ordering, size_t &best_score, size_t &iterations){
    //if(iterations > 1000)
    //    return;

    bool flag = false;

    for (auto it = kVariables.begin(); it != kVariables.end(); it++){
        const unsigned int kVariable = *it;
        //  If indegree is 0 and not yet visited then
        //  only choose that vertex
        if (indegree[kVariable] == 0 && !visited[kVariable]){
            //  reducing indegree of adjacent vertices
            const auto *kChildren = bn->get_child(kVariable);
            const auto kNrChildren = bn->get_child_size(kVariable);
            for(unsigned int i = 0; i < kNrChildren; i++)
                indegree[kChildren[i]]--;

            // including in result
            ordering.push_back(kVariable);
            visited[kVariable] = true;
            all_topological_orders(bn, kVariables, visited, indegree, kBound, ordering, best_ordering, best_score,iterations);

            // resetting visited, res and indegree for
            // backtracking
            visited[kVariable] = false;
            ordering.pop_back();
            for(unsigned int i = 0; i < kNrChildren; i++)
                indegree[kChildren[i]]++;

            flag = true;
        }
    }

    //  We reach here if all vertices are visited.
    //  So we print the solution here
    if (!flag){
        iterations++;
        auto score = kBound.Compute(ordering);
        if(score < best_score){
            best_score = score;
            best_ordering = ordering;
        }
    }
}


void ordering_t::generate_variable_ordering_8(bayesnet *bn, const partition_t &kPartition){
    // all possible topological sorts
    const size_t kNrVariables = bn->get_nr_variables();
    std::vector<bool> visited(kNrVariables,false);

    Bound bound;
    bound.SetBayesnet(bn);
    bound.Init(kPartition);

    std::vector<unsigned int> variables;
    std::copy(kPartition.cutset.begin(),kPartition.cutset.end(), std::back_inserter(variables));
    std::copy(kPartition.set.begin(),kPartition.set.end(), std::back_inserter(variables));

    std::vector<unsigned int> indegree(bn->get_nr_variables(),0);
    for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++)
        indegree[*it] = bn->get_parent_size(*it);

    ordering_t ordering;

    generate_variable_ordering(bn, kPartition);
    size_t iterations = 0;
    size_t best_score = bound.Compute((*this));
    all_topological_orders(bn, variables, visited, indegree, bound, ordering, (*this),best_score, iterations);
}

void ordering_t::generate_variable_ordering(bayesnet *bn, const partition_t &kPartition){
    // algorithm providing topological sort of DAG, DFS and BFs, starting from all roots
    std::vector<unsigned int> roots;
    std::vector<unsigned int> indegree_st(bn->get_nr_variables(),0);
    for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++){
        indegree_st[*it] = bn->get_parent_size(*it);
        if(indegree_st[*it]  == 0)
            roots.push_back(*it);
    }
    for(auto it = kPartition.cutset.begin(); it != kPartition.cutset.end(); it++)
        roots.push_back(*it);

    Bound bound;
    bound.SetBayesnet(bn);
    bound.Init(kPartition);

    assert(!roots.empty());
    printf("Number of roots: %lu\n", roots.size());
    double best_score = std::numeric_limits<double>::max();

    Timer timer;
    timer.Start();
    // DFS
    unsigned int iteration = 0;
    auto dfs_roots = roots;
    do {
        ordering_t ordering;
        auto indegree = indegree_st;
        std::stack<unsigned int> s;
        for(auto it = dfs_roots.rbegin(); it != dfs_roots.rend(); it++)
            s.push(*it);

        while(!s.empty()){
            const unsigned int kNode = s.top();
            s.pop();

            const auto *kChildren = bn->get_child(kNode);
            for(int i = bn->get_child_size(kNode)-1; i >= 0; i--){
                const unsigned int kChild = kChildren[i];

                if(indegree[kChild] != 0){ // is in partition
                    indegree[kChild]--;
                    if(indegree[kChild] == 0)
                        s.push(kChild);
                }
            }
            ordering.push_back(kNode);
        }

        auto score = bound.Compute<Bound::kScore,double>(ordering);
        if(score < best_score){
            best_score = score;
            (*this) = ordering;
        }

        ++iteration;
        timer.Stop();
    } while ( std::next_permutation(dfs_roots.begin(),dfs_roots.end()) && timer.GetDuration<Timer::Seconds>() < 2);
    printf("Topological sorts tried: %lu\n", iteration);
}

void ordering_t::generate_variable_ordering_1(bayesnet *bn,const partition_t &kPartition){
    // breath first, from roots
    assert(kPartition.cutset.empty());
    queue<unsigned int> q;
    for(unsigned int v = 0; v < bn->get_nr_variables(); v++)
        if(bn->get_parent_size(v) == 0)
            q.push(v);

    set <unsigned int> done;
    while(!q.empty()){
        unsigned int v = q.front();
        if(done.find(v) == done.end()){
            for(unsigned int i = 0; i < bn->get_child_size(v); i++)
                q.push(bn->get_child(v)[i]);

            push_back(v);
            done.insert(v);
        }
        q.pop();
    }
}


void ordering_t::generate_variable_ordering_2(bayesnet *bn,const partition_t &kPartition){
    // reverse topological sort of DAG
    generate_variable_ordering(bn,kPartition);
    std::reverse(begin(), end());
}


void ordering_t::generate_variable_ordering_3(bayesnet *bn,const partition_t &kPartition){
    assert(kPartition.cutset.empty());
    // branch and bound
    branchnbound bb;
    bb.init(bn,kPartition);
    bb.bestguess();
    bb.bestfirst();

    (*this) = bb.get_ordering();
}

void ordering_t::generate_variable_ordering_4(bayesnet *bn,const partition_t &kPartition){
    // brute force all combinations
    assert(kPartition.cutset.empty());
    const size_t VARIABLES = bn->get_nr_variables();
    ordering_t var_order;

    for(unsigned int i = 0; i < VARIABLES; i++)
        var_order.push_back(i);

    std::vector<ordering_t> min_orders;
    size_t min_operators = std::numeric_limits<size_t>::max();

    bound_t<false> bound;
    bound.init(bn);

    bbstate state;
    bound.init(state);
    do {
        bound.reset(state);
        bool prune = false;
        size_t upper_bound;

        for(auto it = var_order.begin(); it != var_order.end(); it++){
            const unsigned int variable = *it;
            bound.condition(state, variable);
            // TODO: check for overflow
            upper_bound = bound.get_value(state).get_ui();
            if(upper_bound > min_operators){
                prune = true;
                break;
            }
        }

        if(!prune && upper_bound < min_operators){
            printf("%lu\n", upper_bound);
            min_operators = upper_bound;
            (*this) = var_order;
        } //else if(upper_bound == min_operators)
           // min_orders.push_back(var_order);
    } while (next_permutation(var_order.begin(), var_order.end()));
}

void ordering_t::generate_variable_ordering_5(bayesnet *bn,const partition_t &kPartition){
    // order all weakly connected component topologically
    assert(kPartition.cutset.empty());
    const unsigned int VARIABLES =  bn->get_nr_variables();
    size_t UPPER_BOUND = std::numeric_limits<size_t>::max();
    std::vector<ordering_t> best_orders;


    // find roots
    std::vector<unsigned int> roots;
    std::vector<unsigned int> leaves;
    for(unsigned int v = 0; v < VARIABLES; v++){
        if(bn->get_parent_size(v) == 0)
            roots.push_back(v);
        else if(bn->get_child_size(v) == 0)
            leaves.push_back(v);
    }

    // if we have one root and one leaf, topological sort is the always the best
    if(roots.size() == 1 && leaves.size() == 1){
        generate_variable_ordering(bn,kPartition);
        return;
    }

    // create all combinations of weakly connected components
    const unsigned int COMPONENTS = roots.size();
    bound_t<false> bound;
    bound.init(bn);
    bbstate state;
    bound.init(state);
    do {
        std::vector<unsigned int> component(VARIABLES,0);
        std::vector<unsigned int> encounter(VARIABLES,0);
        std::vector< std::set<unsigned int> > component_parent(COMPONENTS);
        std::vector< std::set<unsigned int> > component_child(COMPONENTS);

        // isolate all weakly components
        for(unsigned int c = 0; c < COMPONENTS; c++){
            queue<unsigned int> q;
            q.push(roots[c]);
            while(!q.empty()){
                unsigned int variable = q.front();
                q.pop();

                if(!component[variable]){
                    encounter[variable]++;
                    component[variable] = c+1;


                    for(unsigned int i = 0; i < bn->get_child_size(variable); i++)
                            q.push(bn->get_child(variable)[i]);

                } else if(component[variable] == c+1)
                    encounter[variable]++;
                else {
                    component_child[c].insert(component[variable]-1);
                    component_parent[component[variable]-1].insert(c);
                }
            }
        }

        // create topological sort of each component
        // NOTE: this creates only one possible topological sort ordering
        std::vector< ordering_t > orders(COMPONENTS);
        //std::vector< std::set<unsigned int> > component_p(COMPONENTS);
        for(unsigned int c = 0; c < COMPONENTS; c++){
            queue<unsigned int> q;
            q.push(roots[c]);
            while(!q.empty()){
                unsigned int variable = q.front();
                q.pop();

                for(unsigned int i = 0; i < bn->get_child_size(variable); i++){
                    unsigned int child = bn->get_child(variable)[i];
                    if(component[child] == c+1){
                        encounter[child]--;
                        if(encounter[child] == 0)
                            q.push(child);
                    }
                }
                orders[c].push_back(variable);
            }
        }


        { // generate topological sort of weakly connected components
            queue<unsigned int> q;
            std::vector<unsigned int> component_order;
            std::vector<unsigned int> component_encounter(COMPONENTS);
            for(unsigned int i = 0; i < COMPONENTS; i++){
                component_encounter[i] = component_parent[i].size();
                if(component_parent[i].size() == 0)
                    q.push(i);
            }

            while(!q.empty()){
                unsigned int c = q.front();
                q.pop();

                // assert(count[v] != 0)
                for(auto it = component_child[c].begin(); it != component_child[c].end(); it++){
                    unsigned int child = *it;

                    component_encounter[child]--;
                    if(component_encounter[child] == 0)
                        q.push(child);
                }
                component_order.push_back(c);
            }

            // determine upper bound for ordering based on weakly connected components
            ordering_t ordering;
            for(unsigned int i = 0; i < COMPONENTS; i++){
                unsigned int c = component_order[i];
                std::copy (orders[c].begin(), orders[c].end(), std::back_inserter(ordering));
            }

            bound.reset();
            for(auto it = ordering.begin(); it != ordering.end(); it++){
                const unsigned int variable = *it;
                bound.condition(variable);
            }

            if(bound.get_value() < UPPER_BOUND){
                // TODO: check for overflow
                UPPER_BOUND = bound.get_value().get_ui();
                best_orders.clear();
            }
            if(bound.get_value() <= UPPER_BOUND)
                best_orders.push_back(ordering);

        }

    } while (next_permutation(roots.begin(), roots.end()));

    (*this) = best_orders[0];
}

string_ordering_t ordering_t::get_variable_string_ordering(bayesnet *bn) const {
    ordering_t ordering = (*this);

    // change variable ordering to variable name ordering
    string_ordering_t string_ordering;
    for(auto it = ordering.begin(); it != ordering.end(); it++)
        string_ordering.push_back(bn->get_node_name(*it));

    return std::move(string_ordering);
}

void ordering_t::generate_variable_ordering_6(bayesnet *bn, const partition_t &kPartition){
    // branch and bound with lookahead
    assert(kPartition.cutset.empty());
    branchnbound bb;
    bb.init(bn,kPartition);
    bb.lookahead(OPT_LOOKAHEAD);

    (*this) = bb.get_ordering();
}

void ordering_t::generate_variable_ordering_7(bayesnet *bn, const partition_t &kPartition){
    // simulated annealing
    generate_variable_ordering_simulated_annealing(bn,kPartition);
}

void ordering_t::generate_variable_ordering_9(bayesnet *bn, const partition_t &kPartition){
    // simulated annealing

    // initial ordering
    if(OPT_SA_READ_ELIM_ORDERING){
        std::string filename;
        if(OPT_PARTITION)
            filename = files.get_filename(ELIM_ORDERING,stringf("%u",kPartition.id));
        else filename = files.get_filename(ELIM_ORDERING);
        this->read(filename, bn);
    } else
        *this = PseudoTree::MinFill(bn,kPartition);

    anneal(bn, kPartition, (*this),2);
}


void ordering_t::generate_variable_ordering_simulated_annealing(bayesnet *bn, const partition_t &kPartition){
    if(OPT_SA_READ_ELIM_ORDERING){
        std::string filename;
        if(OPT_PARTITION)
            filename = files.get_filename(ELIM_ORDERING,stringf("%u",kPartition.id));
        else filename = files.get_filename(ELIM_ORDERING);
        this->read(filename, bn);
    } else
        this->generate_variable_ordering(bn,kPartition);

    anneal(bn, kPartition, (*this),0);
}

void ordering_t::add_weights(bayesgraph &g, const partition_t &kPartition){
    assert(!OPT_DETERMINISM && "Unable to use determinism. Variable literals and weight literals would overlap");
    if(!OPT_ORDER_POST_WEIGHT){
        if(!g.defined())
            throw compiler_ordering_exception("Bayesgraph must be well defined in order to add weights to ordering");

        bayesnet *bn = g.get_bayesnet();
        const unsigned int kNrVariables = g.get_nr_variables();
        std::vector< std::vector<Variable> >variable_to_constraint(kNrVariables);
        std::vector< std::vector<Variable> >constraint_to_variable(kNrVariables);
        for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++){
            const auto kCpt = *it;
            const BayesNode *n = g.get_node(kCpt);
            const auto &kVariables = n->GetVariables();
            constraint_to_variable[kCpt] = n->GetVariables();

            for(auto vit = kVariables.begin(); vit != kVariables.end(); vit++)
                variable_to_constraint[*vit].push_back(kCpt);
        }

        // create ordering with weights
        ordering_t ordering;

        std::vector< std::vector<Variable> > variable_to_sat_value(kNrVariables);
        const auto &kLiteralToVariable = g.get_literal_to_variable();
        const auto &kVariableToLiteralBase = g.get_variable_to_literal();
        for(auto literal_it = begin(); literal_it != end(); literal_it++){
            const auto kLiteral = *literal_it;
            const auto kVariable = kLiteralToVariable[kLiteral];
            const unsigned int kValue = kLiteral-kVariableToLiteralBase[kVariable];
            ordering.push_back(kLiteral);
            variable_to_sat_value[kVariable].push_back(kValue);

            // check which cpts are involved at this level
            auto &cpt_ids = variable_to_constraint[kVariable];
            std::vector<unsigned int> cpts;
            for(auto it = cpt_ids.begin(); it != cpt_ids.end(); it++){
                auto &cpt = *it;
                bool sat = true;
                const auto &kCptVariables = constraint_to_variable[cpt];
                for(auto vit = kCptVariables.begin(); vit != kCptVariables.end(); vit++)
                    if(variable_to_sat_value[*vit].empty())
                        sat = false;

                if(sat)
                    cpts.push_back(cpt);
            }

            // add weights from those cpts
            std::vector<literal_t> weights;
            for(auto it = cpts.begin(); it != cpts.end(); it++){
                auto &cpt = *it;
                BayesNode *n = g.get_node(cpt);
                const auto &kWeights = n->GetWeights();

                XAry xary;
                XAry::Dimension dims;
                XAry::RestrictDimension restrict;
                const auto &kCptVariables = constraint_to_variable[cpt];
                for(auto vit = kCptVariables.begin(); vit != kCptVariables.end(); vit++){
                    dims.push_back(variable_to_sat_value[*vit].size());
                    if(kVariable == *vit){
                        restrict.push_back(true);
                        xary.push_back(kValue); // store dummy variable
                    } else {
                        restrict.push_back(false);
                        xary.push_back(0); // store index of satisified literal
                    }
                }
                xary.SetDimension(dims);

                XAry xary_cpt;
                xary_cpt.SetDimension(bn,kCptVariables);
                do {
                    // convert sat literals to their value and compute weight index
                    for(unsigned int i = 0; i < kCptVariables.size(); i++)
                        xary_cpt[i] = (kVariable == kCptVariables[i]?kValue:variable_to_sat_value[kCptVariables[i]][xary[i]]);

                    unsigned int kWeightIndex = xary_cpt.GetDecimal();
                    weights.push_back(kWeights[kWeightIndex]);
                } while(xary.RestrictedIncrement(restrict));
            }
            std::sort(weights.begin(), weights.end());
            ordering.insert(ordering.end(), weights.begin(), weights.end());

        }

        // remove double weights (local structure)
        std::vector<bool> done(g.get_nr_literals()+g.get_nr_weights()+1, false);
        for(int i = ordering.size() -1; i >= 0; i--){
            literal_t &l = ordering[i];
            if(done[l])
                l = -1;
            else done[l] = true;
        }

        clear();
        for(unsigned int i = 0; i < ordering.size(); i++)
            if(ordering[i] > 0)
                push_back(ordering[i]);
    } else {
        ordering_t &ordering = *this;
        const size_t SIZE = g.get_nr_literals() + g.get_nr_weights() + 1;
        for(unsigned int w = g.get_nr_literals()+1; w < SIZE; w++)
            ordering.push_back(w);

    }
}

ordering_t ordering_t::get_literal_to_variable_ordering(bayesgraph &g){
    ordering_t variable_ordering;
    if(!empty()){
        std::vector<bool> variable_processed(g.get_nr_variables(), false);
        const std::vector<unsigned int>& literal_to_variable =  g.get_literal_to_variable();
        for(auto it = begin(); it != end(); it++){
            unsigned int l = *it;

            // remove weights
            if(l == 0 || l > literal_to_variable.size())
                continue;

            // add variables
            int v = literal_to_variable[l];
            if(variable_ordering.empty() || variable_ordering.back() != v){
                if(!variable_processed[v]){
                    variable_ordering.push_back(v);
                    variable_processed[v] = true;
                } else throw compiler_ordering_exception("Cannot generate variable ordering from literals ordering if they are interleaved");
            }
        }
        if (std::find(variable_ordering.begin(), variable_ordering.end(), false) == variable_ordering.end())
            throw compiler_ordering_exception("Not all variables are represented in literal ordering");

    } else throw compiler_ordering_exception("First read variable variable ordering in order to map them to literals");

    return std::move(variable_ordering);
}

void ordering_t::literal_to_variable_ordering(bayesgraph &g){
    std::vector<literal_t>::operator=(get_literal_to_variable_ordering(g));
}

ordering_t ordering_t::get_variable_to_literal_ordering(bayesgraph &g, const partition_t& kPartition, const bool kIncludeWeightLiterals){
    ordering_t literal_ordering;
    //assert(g.get_nr_variables() == this->size() && "Mismatch in number of variables in ordering and BN");

    if(!empty()){
        for(auto it = begin(); it != end(); it++){
            unsigned int v = *it;
            for(unsigned int i = 0; i < g.get_dimension()[v]; i++){
                literal_t l = g.get_variable_to_literal()[v]+i;
                literal_ordering.push_back(l);
            }
        }

        if(kIncludeWeightLiterals)
            literal_ordering.add_weights(g, kPartition);

    } else throw compiler_ordering_exception("First read variable variable ordering in order to map them to literals");

    return std::move(literal_ordering);
}

void ordering_t::variable_to_literal_ordering(bayesgraph &g, const partition_t& kPartition, const bool kIncludeWeightLiterals){
    std::vector<literal_t>::operator=(get_variable_to_literal_ordering(g,kPartition, kIncludeWeightLiterals));
}

void ordering_t::set_ordering(ordering_t &sub1, ordering_t &sub2, ordering_t &global){
    unsigned int s1 = 0, s2 = 0, gc = 0;
    while(gc < global.size() && (s1 < sub1.size() || s2 < sub2.size())){
        if(global[gc] == sub1[s1]){
            push_back(global[gc]);
            if(global[gc] == sub2[s2])
                s2++;
            s1++;
        } else if(global[gc] == sub2[s2]){
            push_back(global[gc]);
            s2++;
        }
        gc++;
    }
}

void ordering_t::set_ordering(const ordering_t &ordering, const BayesNode *n, const bool kIncludeWeightLiterals){

    // get all variable, literals and dimensions involved
    unordered_set<literal_t> s;
    unsigned int variables = n->GetNrVariables();
    for(unsigned int i = 0; i < variables; i++){
        literal_t l = n->GetLiteral(i);
        unsigned int dim = n->GetDimension(i);
        for(unsigned d = 0; d < dim; d++)
            s.insert(l+d);
    }

    if(kIncludeWeightLiterals){
        const auto &kWeights = n->GetWeights();
        for(unsigned int i = 0; i < n->GetNrWeights(); i++){
            unsigned int w = kWeights[i];
            if(w != 0 && w != 1)
                s.insert(w);
        }
    }

    // make ordering
    clear();
    for(auto it = ordering.begin(); it != ordering.end(); it++){
        literal_t l = *it;
        if(s.find(l) != s.end())
            push_back(l);
    }
}

void biordering_t::set(const ordering_t &ordering){
    back = ordering;
    to.resize(back.size()+1);
    for(unsigned int i = 0; i < ordering.size(); i++){
        if((unsigned int) back[i] >= to.size())
            to.resize(back[i]+1);
        to[back[i]] = i;
    }
}

void ordering_t::stats(bayesnet *bn, const partition_t &kPartition){
    Bound bound;

    bound.SetBayesnet(bn);
    bound.Init(kPartition);
    printf("ordering stats\n");
    printf("==============\n");
    printf("    score     : %.3lf\n", bound.Compute<Bound::kScore,double>((*this)));
    printf("    nodes     : %lu\n", bound.Compute<Bound::kNodes>((*this)));
    printf("    edges     : %lu\n", bound.Compute<Bound::kEdges>((*this)));
    printf("    weights   : %lu\n", bound.Compute<Bound::kWeights>((*this)));
    printf("    ratio     : %.3lf\n", bound.Compute<Bound::kRatio,double>((*this)));
    printf("    operators : %lu\n", bound.Compute<Bound::kOperators>((*this)));
}

