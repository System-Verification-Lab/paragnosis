#include "domainclosure.h"
#include "exceptions.h"
#include <math.h>

#define negated signbit

domain_closure& domain_closure::operator=(domain_closure &closure){
    g = closure.g;
    ALO = closure.ALO;
    AMO = closure.AMO;
    total_unsat_variables = closure.total_unsat_variables;
    return *this;
}

void domain_closure::init(bayesgraph *g){
    if(!g || !g->defined())
        throw compiler_debug_exception("Domain closure requires bayesgraph to be defined");

    this->g = g;
    const unsigned int VARIABLES = g->get_nr_variables();
    ALO = g->get_dimension();
    AMO.resize(VARIABLES);
    fill(AMO.begin(), AMO.end(), false);
    #ifdef DEBUG
    total_unsat_variables = VARIABLES;
    #endif
}

bool domain_closure::is_redundant(literal_t l){
    const std::vector<unsigned int> &literal_to_variable = g->get_literal_to_variable();
    variable_t variable = literal_to_variable[abs(l)];
    if(AMO[variable])
        return true;
    else return false;

}

satisfy_t domain_closure::condition(literal_t l){
    #ifdef DEBUG
    if(total_unsat_variables == 0)
        throw compiler_debug_exception("All variable domains closed prior to conditioning");
    #endif
    const std::vector<unsigned int> &literal_to_variable = g->get_literal_to_variable();

    variable_t variable = literal_to_variable[abs(l)];
    if(AMO[variable]){
        return redundant;
    } else {
        if(negated(l)) {
            unsigned int &closure = ALO[variable];
            closure--;

            if(closure == 0)
                return unsatisfiable;
            else return satisfiable;
        } else {
            AMO[variable] = true;
            #ifdef DEBUG
            total_unsat_variables--;
            #endif
            return satisfiable;
        }
    }
}

void domain_closure::undo(literal_t l){
    const std::vector<unsigned int> &literal_to_variable = g->get_literal_to_variable();

    if(negated(l)) {
        unsigned int &closure = ALO[literal_to_variable[-1*l]];
        closure++;
    } else {
        AMO[literal_to_variable[l]] = false;
        #ifdef DEBUG
        total_unsat_variables++;
        #endif
    }
}

