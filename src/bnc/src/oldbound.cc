#include "oldbound.h"
#include <cassert>

using namespace bnc;

void bbstate::allocate(const unsigned int VARIABLES){
    constraint.resize(VARIABLES);
    lock.resize(VARIABLES);
}

bbstate& bbstate::operator=(bbstate &state){
    constraint = state.constraint;
    lock = state.lock;
    history = state.history;
    span_variables = state.span_variables;
    value = state.value;
    return (*this);
}

bool bbstate::operator!=(bbstate &state){
    if(constraint == state.constraint
        && lock == state.lock
        && history == state.history
        && span_variables == state.span_variables
        && value == state.value)
        return true;
    else return false;
}


template <bool EXT>
bound<EXT>::bound(){
    bn = NULL;
}

template <bool EXT>
void bound<EXT>::init(bbstate &state, bayesnet *bn){
    const unsigned int VARIABLES = bn->get_nr_variables();
    this->VARIABLES = VARIABLES;
    this->bn = bn;

    variable_to_constraint.resize(VARIABLES);
    constraint_to_variable.resize(VARIABLES);
    for(unsigned int i = 0; i < VARIABLES; i++){
        variable_to_constraint[i].push_back(i);
        constraint_to_variable[i].push_back(i);

        const uint32_t* parent = bn->get_parent(i);
        const unsigned int PARENTS = bn->get_parent_size(i);
        for(unsigned int j = 0; j < PARENTS; j++){
            variable_to_constraint[parent[j]].push_back(i);
            constraint_to_variable[i].push_back(parent[j]);
        }
    }

    init(state);
}

template <bool EXT>
void bound<EXT>::init(bayesnet *bn){
    init(state, bn);
}

template <bool EXT>
void bound<EXT>::init(bbstate &state){
    state.allocate(VARIABLES);
    reset(state);
}

template <bool EXT>
void bound<EXT>::reset(bbstate &state){
    while(EXT && !state.history.empty())
        state.history.pop();

    while(!state.value.empty())
        state.value.pop();

    state.value.push(0);

    state.span_variables.clear();
    std::fill(state.lock.begin(), state.lock.end(), 0);
    for(unsigned int i = 0; i < constraint_to_variable.size(); i++)
        state.constraint[i] = constraint_to_variable[i].size();
}

template <bool EXT>
void bound<EXT>::reset(){
    reset(state);
}

template <bool ext>
void bound<ext>::set_partition(const partition_t &partition, map_backup_t &backup){
    backup.resize(partition.cutset.size());
    unsigned int i = 0;
    for(auto it = partition.cutset.begin(); it != partition.cutset.end(); it++, i++){
        unsigned int variable = *it;

        // update variable_to_constraint to include only partitioned constraints
        std::vector<unsigned int> &constraints = variable_to_constraint[variable];
        backup[i] = constraints;
        constraints.clear();
        for(auto it = backup[i].begin(); it != backup[i].end(); it++){
            unsigned int constraint = *it;
            if(partition.set.find(constraint) != partition.set.end())
                constraints.push_back(constraint);
        }
    }
}

template <bool ext>
void bound<ext>::unset_partition(const partition_t &partition, map_backup_t &backup){
    unsigned int i = 0;
    for(auto it = partition.cutset.begin(); it != partition.cutset.end(); it++, i++){
        unsigned int variable = *it;

        std::vector<unsigned int> &constraints = variable_to_constraint[variable];
        constraints = backup[i];
    }
}

template <bool EXT>
void bound<EXT>::condition(bbstate &state, const unsigned int variable){
    // compute level width
    const uint32_t *dimension = bn->get_states();
    state.span_variables.insert(variable);

    mpz_class edges = 1;
    for(auto it = state.span_variables.begin(); it != state.span_variables.end(); it++){
        unsigned int variable = *it;
        edges *= dimension[variable];
    }

    // update span variables
    std::vector<unsigned int> &constraints = variable_to_constraint[variable];
    size_t weights = 0;
    for(auto it = constraints.begin(); it != constraints.end(); it++){
        const unsigned int i = *it;
        state.constraint[i]--;
        state.lock[variable]++;
        if(state.constraint[i] == 0){
            weights++;

            // remove span variables
            std::vector<unsigned int> &variables = constraint_to_variable[i];
            for(auto itc = variables.begin(); itc != variables.end(); itc++){
                const unsigned int v = *itc;
                state.lock[v]--;
                if(v != variable){
                    if(state.lock[v] == 0)
                        state.span_variables.erase(v);
                }
            }
        }
    }
    if(state.lock[variable] == 0)
        state.span_variables.erase(variable);

    // compute current size
    edges *= 2+weights;
    if(EXT){
        state.history.push(variable);
        state.value.push(state.value.top() + edges);
    } else state.value.top() = state.value.top() + edges;
}


template <bool EXT>
void bound<EXT>::condition(const unsigned int variable){
    condition(state, variable);
}

template <bool EXT>
void bound<EXT>::set_value(size_t value){
    return set_value(state, value);
}

template <bool EXT>
void bound<EXT>::set_value(bbstate& state, size_t value){
    state.value.top() = value;
}

template <bool EXT>
mpz_class bound<EXT>::get_value(){
    return get_value(state);
}

template <bool EXT>
mpz_class bound<EXT>::get_value(bbstate& state){
    return state.value.top();
}

template <bool EXT>
void bound<EXT>::undo(bbstate &state, unsigned int l){
    std::vector<unsigned int> &constraints = variable_to_constraint[l];
    for(auto it = constraints.begin(); it != constraints.end(); it++){
        const unsigned int constraint = *it;
        unsigned int &c = state.constraint[constraint];
        if(c == 0){
            std::vector<unsigned int> &variables = constraint_to_variable[constraint];
            for(auto vit = variables.begin(); vit != variables.end(); vit++){
                const unsigned int variable = *vit;
                if(variable != l){
                    if(state.lock[variable] == 0)
                        state.span_variables.insert(variable);

                    state.lock[variable]++;
                }
            }
        } else state.lock[l]--;
        c++;
    }
    state.span_variables.erase(l);
}

template <bool EXT>
void bound<EXT>::undo(unsigned int l){
    undo(state, l);
}

template <bool EXT>
void bound<EXT>::undo(bbstate &state){
    if(!state.history.empty()){
        undo(state, state.history.top());
        state.history.pop();
    }
}

template <bool EXT>
void bound<EXT>::undo(){
    undo(state);
}

template class bound<false>;
template class bound<true>;


