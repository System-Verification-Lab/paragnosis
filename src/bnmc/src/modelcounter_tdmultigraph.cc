#include <unistd.h>
#include <algorithm>
#include <stack>
#include <bn-to-cnf/config.h>
#include <bn-to-cnf/exceptions.h>
#include <bnc/bayesgraph.h>
#include <bnc/exceptions.h>
#include "modelcounter.h"
#include "io.h"
#include "options.h"
#include "exceptions.h"
#include <climits>
#include "mstack.h"
#include "debug.h"
#include "multigraph.h"

namespace bnmc {

using namespace std;
using namespace mstack;

template <>
template <>
probability_t ModelCounter<ModelType::TDMULTIGRAPH>::Traverse<1>(Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList){
    const unsigned int kTier = 0;
    const unsigned int kPartition = 0;

    MultiGraph::Circuit &kCircuit = partition_[kPartition].tdmgc_;
    const size_t kNrNodes = kCircuit.GetSize();
    DynamicArray<probability_t> probabilities(kNrNodes);
    //std::fill(probabilities.begin(), probabilities.begin()+kNrNodes, kNotTraversed);
    bool traversed[kNrNodes] = {false};
    for(unsigned int i = 0; i < kCircuit.GetNrTerminals(); i++){
        traversed[i] = true;
        probabilities[i] = 1;
    }

    std::stack<const MultiGraph::Node*> s;
    s.push(kCircuit.GetRoot());

    while(!s.empty()){
        // depth first traversal
        const MultiGraph::Node *kNode = s.top();

        if(traversed[kCircuit.GetIndex(kNode)]){
            // compute probability based on children
            s.pop();

            probability_t &probability = probabilities[kCircuit.GetIndex(kNode)];
            if(kNode->IsAnd()){  // AND node
                probability = 1;
                const MultiGraph::Edge* kEnd = kNode->EdgeEnd();
                const MultiGraph::Edge* kEdge = kNode->EdgeBegin();
                while(kEdge != kEnd){
                    probability *= probabilities[kCircuit.GetIndex(kEdge->to)];
                    ++kEdge;
                }
            } else { // OR node
                probability = 0;
                if(kConditionTierList[kNode->variable] <= kTier){ // conditioned OR node
                    const MultiGraph::Edge *kEdge = &(kNode->edges[kEvidenceList[kNode->variable]]);
                    probability += kEdge->probability * probabilities[kCircuit.GetIndex(kEdge->to)];
                } else { // unconditioned OR node
                    const MultiGraph::Edge* kEnd = kNode->EdgeEnd();
                    const MultiGraph::Edge* kEdge = kNode->EdgeBegin();
                    while(kEdge != kEnd){
                        probability += kEdge->probability * probabilities[kCircuit.GetIndex(kEdge->to)];
                        ++kEdge;
                    }
                }
            }
        } else {
            // push children on stack
            traversed[kCircuit.GetIndex(kNode)] = true;
            if(kNode->IsAnd() || !(kConditionTierList[kNode->variable] <= kTier)){
                const MultiGraph::Edge* kEnd = kNode->rEdgeEnd();
                const MultiGraph::Edge* kEdge = kNode->rEdgeBegin();
                while(kEdge != kEnd){
                    if(!traversed[kCircuit.GetIndex(kEdge->to)])
                        s.push(kEdge->to);
                    --kEdge;
                }
            } else {
                const MultiGraph::Edge *kEdge = &(kNode->edges[kEvidenceList[kNode->variable]]);
                if(!traversed[kCircuit.GetIndex(kEdge->to)])
                    s.push(kEdge->to);
            }
        }
    }
    return probabilities[kCircuit.GetRootIndex()];
}

template <>
template <>
probability_t ModelCounter<ModelType::TDMULTIGRAPH>::Posterior<1>(){
    probability_t p, pq;

    pq = Traverse<1>(cache_,evidence_list_,condition_tier_);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["MULTINODE"];
    probs.p = 1;
    probs.pq = pq;
    #endif

    if(has_query_variable_){
        p = Traverse<1>(cache2_,evidence_list2_,condition_tier2_);
        #ifdef DEBUG
        probs.p = p;
        #endif

        if(p == 0)
            return -1;
        else return pq/p;
    } else return pq;
}

}

