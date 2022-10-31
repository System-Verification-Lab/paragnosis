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

struct StackItem {
    double &probability;
    const double kWeight;
    const MultiGraph::Node *kNode;

    template<class A, class B, class C>
    StackItem(A& probability, B kWeight, C kNode) : probability(probability), kWeight(kWeight), kNode(kNode) {};
};

template <>
template <>
probability_t ModelCounter<ModelType::MULTIGRAPH>::Traverse<2>(Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList){
    const unsigned int kTier = 0;
    const unsigned int kPartition = 0;

    MultiGraph::Circuit &kCircuit = partition_[kPartition].mgc_;
    const auto kNrNodes = kCircuit.GetSize();
    const auto kNrVariables = kConditionTierList.size();
    ProbabilityList probabilities(kNrNodes,0);

    //std::fill(probabilities.begin(), probabilities.begin()+kNrNodes, kNotTraversed);
    bool traversed[kNrNodes] = {false};
    for(unsigned int i = 0; i < kCircuit.GetNrTerminals(); i++){
        traversed[i] = true;
        probabilities[i] = 1;
    }

    probability_t result = 0;
    const probability_t kOne = 1;

    //std::stack<const MultiNode*> s;
    //s.push(kCircuit.GetRoot());
    std::vector < StackItem > s;
    s.reserve(kNrVariables*2);
    s.emplace_back(result, kOne, kCircuit.GetRoot());

    while(!s.empty()){
        // depth first traversal
        auto &kItem = s.back();
        const MultiGraph::Node * const kNode = kItem.kNode;
        const auto kIndex = kCircuit.GetIndex(kNode);

        if(traversed[kIndex]){
            // compute probability based on children
            s.pop_back();
            kItem.probability += kItem.kWeight * probabilities[kIndex];
        } else {
            traversed[kIndex] = true;
            probability_t &probability = probabilities[kIndex];
            probability = 0;

            // push children on stack
            const bool kIsConditioned = kConditionTierList[kNode->variable] <= kTier;
            if(kIsConditioned){
                const MultiGraph::Edge *kEdge = &(kNode->edges[kEvidenceList[kNode->variable]]);
                s.emplace_back(probability,kEdge->probability,kEdge->to);
            } else {
                const MultiGraph::Edge* kEnd = kNode->rEdgeEnd();
                const MultiGraph::Edge* kEdge = kNode->rEdgeBegin();
                while(kEdge != kEnd){
                    s.emplace_back(probability,kEdge->probability,kEdge->to);
                    --kEdge;
                }
            }
        }
    }
    return result;
}

template <>
template <>
probability_t ModelCounter<ModelType::MULTIGRAPH>::Posterior<2>(){
    probability_t p, pq;

    pq = Traverse<2>(cache_,evidence_list_,condition_tier_);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["MULTIGRAPH2"];
    probs.p = 1;
    probs.pq = pq;
    #endif

    if(has_query_variable_){
        p = Traverse<2>(cache2_,evidence_list2_,condition_tier2_);
        #ifdef DEBUG
        probs.p = p;
        #endif

        if(p == 0)
            return -1;
        else return pq/p;
    } else return pq;
}



}

