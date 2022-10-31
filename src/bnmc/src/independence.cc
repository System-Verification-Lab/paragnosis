#include "independence.h"
#include "options.h"
#include <algorithm>
#include <climits>

#define EVIDENCE UINT_MAX

namespace bnmc {

Independence::Independence(){
    nr_components_ = 0;
}

void Independence::DetermineAncestorNeighbors(std::queue<unsigned int> &queue, const unsigned int kComponent){
    bayesnet *bn = manager.bn;
    unsigned int current = queue.front();
    uint32_t *parents = bn->get_parent(current);
    unsigned int nr_parents = bn->get_parent_size(current);
    for(unsigned int i = 0; i < nr_parents; i++){
        unsigned int variable = parents[i];
        unsigned int component = variable_to_component_[variable];
        if(component == kComponent)
            queue.push(variable);
    }
}

void Independence::DetermineNeighbors(std::queue<unsigned int> &queue, const unsigned int kComponent){
    bayesnet *bn = manager.bn;
    unsigned int current = queue.front();
    uint32_t *children = bn->get_child(current);
    uint32_t *parents = bn->get_parent(current);
    unsigned int nr_children = bn->get_child_size(current);
    unsigned int nr_parents = bn->get_parent_size(current);

    for(unsigned int i = 0; i < nr_parents; i++){
        unsigned int variable = parents[i];
        unsigned int component = variable_to_component_[variable];
        if(component == 0)
            queue.push(variable);
        else if(component == EVIDENCE)
            component_evidence_[kComponent].insert(variable);
    }

    for(unsigned int i = 0; i < nr_children; i++){
        unsigned int variable = children[i];
        unsigned int component = variable_to_component_[variable];
        if(component == 0)
           queue.push(variable);
        else if(component == EVIDENCE){ // support for moralized parents
            component_evidence_[kComponent].insert(variable);

            uint32_t *parents = bn->get_parent(variable);
            unsigned int nr_parents = bn->get_parent_size(variable);
            for(unsigned int i = 0; i < nr_parents; i++){
                unsigned int variable = parents[i];
                unsigned int component = variable_to_component_[variable];
                if(component == 0)
                    queue.push(variable);
                else if(component == EVIDENCE)
                    component_evidence_[kComponent].insert(variable);
            }
        }
    }
}

void Independence::DetermineComponent(unsigned int root){
    std::queue<unsigned int> queue;

    if(variable_to_component_[root] == 0){
        queue.push(root);
        nr_components_++;
    }

    const unsigned int kComponent = nr_components_;
    component_evidence_.resize(kComponent+1);

    while(!queue.empty()){
        unsigned int i = queue.front();
        variable_to_component_[i] = kComponent;
        DetermineNeighbors(queue, kComponent);
        queue.pop();
    }
}

void Independence::DetermineConditionalIndependence(const Evidence &evidence){
    bayesnet *bn = manager.bn;

    // clear evidence variables that border components
    component_evidence_.clear();
    evidence_variables_.clear();
    nr_components_ = 0;

    // initialize variable component identification
    // set evidence to EVIDENCE (1), other to 0
    variable_to_component_.resize(bn->get_nr_variables());
    std::fill(variable_to_component_.begin(), variable_to_component_.end(), 0);
    const EvidenceVariableSet &evidence_set = evidence.GetEvidenceVariableSet();
    for(auto it = evidence_set.begin(); it != evidence_set.end(); it++){
        if(!evidence.IsQueryVariable(it)){
            variable_to_component_[it->variable] = EVIDENCE;
            evidence_variables_.insert(it->variable);
        }
    }

    // determine which variable belongs to which component
    for(unsigned int i = 0; i < variable_to_component_.size(); i++)
        if(variable_to_component_[i] == 0)
            DetermineComponent(i);

}

VariableSet Independence::GetDependentVariables(const Evidence &kEvidence){
    DetermineConditionalIndependence(kEvidence);
    return std::move(GetDependentVariables(kEvidence.GetQueryVariable()));
}

VariableSet Independence::GetConditionalDependentVariables(const Evidence &kEvidence){
    DetermineConditionalIndependence(kEvidence);
    return std::move(GetConditionalDependentVariables(kEvidence.GetQueryVariable()));
}

VariableSet Independence::GetConditionalDependentVariables(const variable_t kQueryVariable){
    const unsigned int kComponent = variable_to_component_[kQueryVariable];
    VariableSet dependent;
    for(unsigned int variable = 0; variable < variable_to_component_.size(); variable++)
        if(variable_to_component_[variable] == kComponent)
            dependent.insert(dependent.end(), variable);

    return std::move(dependent);
}

VariableSet Independence::GetDependentVariables(const variable_t kQueryVariable){
    // PRECONDITION: run DetermineConditionalIndependence() with evidence first
    const unsigned int kComponent = variable_to_component_[kQueryVariable];
    std::queue<unsigned int> queue;
    queue.push(kQueryVariable);

    VariableSet dependent;
    VariableSet &variables = component_evidence_[kComponent];
    for(auto it = variables.begin(); it != variables.end(); it++)
        queue.push(*it);

    while(!queue.empty()){
        unsigned int i = queue.front();
        if(dependent.find(i) == dependent.end()){
            dependent.insert(i);

            DetermineAncestorNeighbors(queue, kComponent);
        }
        queue.pop();
    }
    return std::move(dependent);
}

bool Independence::ConditionallyIndependent(const variable_t a, const variable_t b){
    // PRECONDITION: run DetermineConditionalIndependence() with evidence first
    if(variable_to_component_[a] == EVIDENCE || variable_to_component_[a] == EVIDENCE){
        if(variable_to_component_[a] == EVIDENCE && variable_to_component_[a] == EVIDENCE){
            // check if the border the same component
            for(auto it = component_evidence_.begin(); it != component_evidence_.end(); it++)
                if(it->find(a) != it->end() && it->find(b) != it->end())
                    return false;
            return true;
        } else {
            // check if either one borders the component of the other
            if(variable_to_component_[a] == EVIDENCE){
                const unsigned int kComponent = variable_to_component_[b];
                return component_evidence_[kComponent].find(a) == component_evidence_[kComponent].end();
            } else {
                const unsigned int kComponent = variable_to_component_[a];
                return component_evidence_[kComponent].find(b) == component_evidence_[kComponent].end();
            }
        }
    } else {
        // check if the are in the same component
        return variable_to_component_[a] != variable_to_component_[b];
    }
}

VariableSet Independence::GetConditionalIndependentEvidenceVariables(const Evidence &kEvidence){
    DetermineConditionalIndependence(kEvidence);
    return std::move(GetConditionalIndependentEvidenceVariables(kEvidence.GetQueryVariable()));
}

VariableSet Independence::GetConditionalIndependentEvidenceVariables(const variable_t kQueryVariable){
    const unsigned int kComponent = variable_to_component_[kQueryVariable];
    VariableSet &component_evidence_variables = component_evidence_[kComponent];
    VariableSet independent;

    std::set_difference(evidence_variables_.begin(), evidence_variables_.end(),
        component_evidence_variables.begin(), component_evidence_variables.end(),
        std::inserter(independent, independent.end()));

    return std::move(independent);
}

}

