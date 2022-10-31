#ifndef BNMC_INCLUDE_INDEPENDENCE_H_
#define BNMC_INCLUDE_INDEPENDENCE_H_

#include <queue>
#include <set>
#include <bn-to-cnf/bayesnet.h>
#include "evidence.h"
#include "types.h"

namespace bnmc {

class Independence {
    public:
        Independence();
        void DetermineConditionalIndependence(const Evidence &);

        VariableSet GetDependentVariables(const Evidence &);
        VariableSet GetDependentVariables(const variable_t);
        VariableSet GetConditionalDependentVariables(const Evidence &);
        VariableSet GetConditionalDependentVariables(const variable_t);

        bool ConditionallyIndependent(const variable_t, const variable_t);
        VariableSet GetConditionalIndependentEvidenceVariables(const Evidence &);
        VariableSet GetConditionalIndependentEvidenceVariables(const variable_t);
    private:
        void DetermineComponent(unsigned int);
        void DetermineNeighbors(std::queue<unsigned int>&, const unsigned int);
        void DetermineAncestorNeighbors(std::queue<unsigned int>&, const unsigned int);

        std::vector<unsigned int> variable_to_component_;
        std::vector< VariableSet > component_evidence_;
        VariableSet evidence_variables_;
        unsigned int nr_components_;
};

}

#endif
