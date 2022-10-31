#include "sat.h"
#include "misc.h"
#include <math.h>
#include <algorithm>
#include <map>
#include "options.h"
#include "xary.h"

#define negated signbit

using namespace bnc;
using namespace std;

sat::sat(){
    mapped = false;
    W = NULL;
    manager = NULL;
}


sat::~sat(){
    destroy_weights(manager, W);
}

void sat::set_manager(manager_t *manager){
    this->manager = manager;
}

node::weights* sat::get_weights(){
    node::weights *rW = W;
    W = NULL;
    return rW;
}


bool sat::weighted(){
    return W != NULL;
}

void sat::print(bool hist){
    hist = hist;
}

template <bool DETERMINISM>
void sat::init(BayesNode *n){

    // get all variable, literals and dimensions involved
    const auto kNrLiterals = n->GetNrVariables();
    literal_t literals[kNrLiterals];
    auto bnliterals = n->GetLiterals();
    for(unsigned int i = 0; i < kNrLiterals; i++){
        literal_t l = bnliterals[i];
        if(mapped)
            literals[i] = literal_map[l];
        else
            literals[i] = l;
    }

    // determine which literal is present in which clauses
    XAry xary;
    xary.SetDimension(n->GetDimensions());
    xary.SetZero();
    const auto &kWeights = n->GetWeights();
    literal_to_clause.resize(LITERALS+1);
    size_t clause_nr = 0;
    do {
        unsigned int w = kWeights[clause_nr++];

        clause_to_weight.push_back(w);
        unsat_literals.push_back(kNrLiterals);
        for(unsigned int i = 0; i < kNrLiterals; i++){
            literal_t l = literals[i]+xary[i];
            literal_to_clause[l].push_back(CLAUSES);
        }
        CLAUSES++;
    } while(xary.Increment());
}

template <bool DETERMINISM>
void sat::init(bayesgraph &g, int variable){
    LITERALS = 0;
    VARIABLES = 0;
    CLAUSES = 0;


    if(variable >= 0){
        // assume g encodes partial bn only, thus create mapping
        literal_to_variable.push_back(0); // do not use index 0 as literal
        BayesNode *n = g.get_node(variable);
        vector <literal_t> literals;
        vector <unsigned int> dims;
        for(unsigned int i = 0; i < n->GetNrVariables(); i++){
            literal_t l = n->GetLiteral(i);
            unsigned int dim = n->GetDimension(i);
            dims.push_back(dim);

            if(literal_map.find(l) == literal_map.end()){
                variable_to_literal.push_back(LITERALS+1);
                for(unsigned d = 0; d < dim; d++){
                    literal_map[l + d] = LITERALS+1;
                    literal_map[-1*(l + d)] = -1*(LITERALS+1);
                    literal_to_variable.push_back(VARIABLES);
                    LITERALS++;
                }
                VARIABLES++;
                dimension.push_back(dim);
            }
            literals.push_back(literal_map[l]);
        }
        mapped = true;

        init<DETERMINISM>(g.get_node(variable));
    } else {

        if(g.defined()){
            // no mapping neccesary
            variable_to_literal = g.get_variable_to_literal();
            literal_to_variable = g.get_literal_to_variable();
            dimension = g.get_dimension();
            LITERALS = literal_to_variable.size();
            VARIABLES = variable_to_literal.size();
        } else {
            // assume g encodes partial bn only, thus create mapping
            literal_to_variable.push_back(0); // do not use index 0 as literal
            for(unsigned int i = 0; i < g.size(); i++){
                BayesNode *n = g.get_node(i);
                vector <literal_t> literals;
                vector <unsigned int> dims;
                for(unsigned int i = 0; i < n->GetNrVariables(); i++){
                    literal_t l = n->GetLiteral(i);
                    unsigned int dim = n->GetDimension(i);
                    dims.push_back(dim);

                    if(literal_map.find(l) == literal_map.end()){
                        variable_to_literal.push_back(LITERALS+1);
                        for(unsigned d = 0; d < dim; d++){
                            literal_map[l + d] = LITERALS+1;
                            literal_map[-1*(l + d)] = -1*(LITERALS+1);
                            literal_to_variable.push_back(VARIABLES);
                            LITERALS++;
                        }
                        VARIABLES++;
                        dimension.push_back(dim);
                    }
                    literals.push_back(literal_map[l]);
                }

            }
            mapped = true;
        }

        for(unsigned int i = 0; i < g.size(); i++)
            init<DETERMINISM>(g.get_node(i));
    }

    // DEBUG: VARIALBES == variable_to_literal.size() == dimension.size();
    // DEBUG: CLAUSES == clauses_to_weight.size()
    // DEBUG: LITERALS == literal_to_variable.size() == literal_to_clause.size()

    total_unsat_clauses = CLAUSES;
    total_unsat_variables = VARIABLES;

    sat_variable.resize(VARIABLES);
    literal_constraint.resize(LITERALS+1);
    sat_literal_clause.resize(CLAUSES);

    literals_left = dimension;
    fill(sat_variable.begin(), sat_variable.end(), false);
    fill(literal_constraint.begin(), literal_constraint.end(), false);
    fill(sat_literal_clause.begin(), sat_literal_clause.end(), 0);
}

satisfy_t sat::impose_constraint(literal_t l){
    literal_t L = abs(l);
    variable_t v = literal_to_variable[L];
    if(sat_variable[v]){
        return redundant;
    } else {
        history.resize(history.size()+1);
        history.back().push_back(l);
        #ifdef DEBUG
        if(literal_constraint[abs(l)])
            fprintf(stderr, "already conditioned on literal %d\n", l);
        #endif
        literal_constraint[abs(l)] = true;
        if(negated(l)) {
            unsigned int &count = literals_left[literal_to_variable[L]];
            count--;

            if(count == 0)
                return unsatisfiable;
            else return unsatisfied;
        } else {
            for(unsigned int i = 0; i < dimension[v]; i++){
                literal_t l = variable_to_literal[v] + i;
                if(!literal_constraint[l])
                    history.back().push_back(-1*l);
            }
            sat_variable[v] = true;
            total_unsat_variables--;
            return satisfiable;
        }
    }
}

void sat::release_constraint(literal_t l){
    literal_constraint[abs(l)] = false;
    if(negated(l)) {
        unsigned int &count = literals_left[literal_to_variable[abs(l)]];
        count++;
    } else {
        sat_variable[literal_to_variable[l]] = false;
        total_unsat_variables++;
    }
}

template <bool DETERMINISM>
satisfy_t sat::condition(literal_t l){
    if(mapped){
        #ifdef DEBUG
        if(literal_map.find(l) == literal_map.end()){
            fprintf(stderr, "Literal %d not mapped\n", l);
            exit(1);
        }
        #endif
        l = literal_map[l];
    }

    satisfy_t c_l = impose_constraint(l);
    if(c_l != redundant){
        for(auto it = history.back().begin(); it != history.back().end(); it++){
            literal_t l = *it;
            if(negated(l)){
                // condition f on negated l
                std::vector<literal_t> &l_to_c = literal_to_clause[-1*l];
                for(auto it = l_to_c.begin(); it != l_to_c.end(); it++){
                    unsigned int c = *it;
                    literal_t &sat_l = sat_literal_clause[c];

                    if(!sat_l){
                        sat_l = l;
                        total_unsat_clauses--;
                    }
                }
            } else {
                // condition f on l
                std::vector<literal_t> &l_to_c = literal_to_clause[l];
                for(auto it = l_to_c.begin(); it != l_to_c.end(); it++){
                    unsigned int c = *it;

                    literal_t &sat_l = sat_literal_clause[c];
                    if(!sat_l){
                        uint8_t &unsat_l = unsat_literals[c];
                        unsat_l--;
                        if(unsat_l == 0){
                            total_unsat_clauses--;
                            sat_l = l;

                            weight_t w = clause_to_weight[c];
                            if(DETERMINISM){
                                if(w == 0)
                                    c_l |= unsatisfiable;
                                else if(w != 1){
                                    if(!W) W = bnc::create_weights(manager);
                                    W->insert(w);
                                }
                            } else {
                                if(!W) W = bnc::create_weights(manager);
                                W->insert(w);
                            }
                        }
                    }
                }
            }
        }
        #ifdef DEBUG
        if((total_unsat_clauses != 0 && total_unsat_variables == 0))
            fprintf(stderr, "inconsistency between theory and sat solver\n");
        #endif
        if(total_unsat_clauses == 0)
            return satisfiable | c_l;
        else return unsatisfied | c_l;
    } else return c_l; // redundant
}

void sat::undo(){
    if(!history.empty()){
        vector<literal_t> &L = history.back();
        if(W) delete W;

        release_constraint(L[0]);
        for(auto it = L.rbegin(); it != L.rend(); it++){
            literal_t l = *it;

            if(negated(l)){
                // sat for negated l
                std::vector<literal_t> &l_to_c = literal_to_clause[-1*l];
                for(auto it = l_to_c.begin(); it != l_to_c.end(); it++){
                    unsigned int c = *it;

                    literal_t &sat_l = sat_literal_clause[c];
                    if(sat_l == l){
                        sat_l = 0;
                        total_unsat_clauses++;
                    }
                }
            } else {
                // sat for l
                std::vector<literal_t> &l_to_c = literal_to_clause[l];
                for(auto it = l_to_c.begin(); it != l_to_c.end(); it++){
                    unsigned int c = *it;

                    literal_t &sat_l = sat_literal_clause[c];
                    uint8_t &unsat_l = unsat_literals[c];
                    #ifdef DEBUG
                    if(!sat_l && !unsat_l)
                        fprintf(stderr, "cannot have unsatisfied clause with no literals\n");
                    #endif
                    if (sat_l == l && unsat_l == 0){
                        sat_l = 0;
                        total_unsat_clauses++;
                        unsat_l++;
                    } else if(!sat_l)
                        unsat_l++;

                }
            }
        }
        history.pop_back();
    }
}

template satisfy_t sat::condition<false>(literal_t);
template satisfy_t sat::condition<true>(literal_t);

template void sat::init<false>(bayesgraph&, int);
template void sat::init<true>(bayesgraph&, int);
template void sat::init<false>(BayesNode*);
template void sat::init<true>(BayesNode*);

