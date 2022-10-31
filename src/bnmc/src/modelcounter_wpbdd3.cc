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
#include <cassert>

namespace bnmc {

using namespace std;
using namespace mstack;


#define N 3

template <>
template <>
probability_t ModelCounter<ModelType::WPBDD>::Traverse<N>(Cache &cache, const EvidenceList &kEvidenceList,const ConditionTierList &kConditionTierList){
    const unsigned int kTier = 0;
    const unsigned int kPartitionId = 0;
    const Partition::ArithmeticCircuit &kCircuit = partition_[kPartitionId].ac_;
    //std::vector<probability_t> probabilities(kCircuit.size());
    //std::fill(probabilities.begin(), probabilities.end(), kNotTraversed);
    probability_t* probabilities = (probability_t*) alloca(sizeof(probability_t)*kCircuit.size());
    assert(probabilities);
    std::fill_n(probabilities,kCircuit.size(), kNotTraversed);

    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    #ifdef DEBUG
    bool kWriteDot = false;
    int kMaxDepth = -1;
    //std::vector<unsigned int> identity(probabilities.size(),0);
    std::string filename = "debug.wpbdd.dot";
    #endif

    StackTop i;
    StackTopInit(i);
    //Stack s(cache.GetStackSize(kPartition));
    NodeIndex* s = (NodeIndex*) alloca(sizeof(NodeIndex)*cache.GetStackSize(kPartitionId));
    assert(s);
    Push(s,i,kRootIndex,kCircuit,kEvidenceList,kConditionTierList);

    traverse: {
        const unsigned int &kParentIndex = s[i-2];
        const unsigned int &kIndex = s[i-1];
        probability_t &probability = probabilities[kIndex];

        #ifdef DEBUG
        //if(kWriteDot)
        //    WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex, kMaxDepth);
        #endif

        if(probability == kNotTraversed){
            //#ifdef DEBUG
            //if(i >= s.size()-1)
            //    throw ModelCounterException("Custom stack is going out of bounds");
            //#endif
            probability = 0;
            Push(s,i,kIndex,kCircuit,kEvidenceList,kConditionTierList);

            goto traverse;
        }

        const Partition::Node &kParentNode = kCircuit[kParentIndex];
        const bool kIsPositiveCofactor = kParentNode.t == kIndex;

        if(kIsPositiveCofactor){
            // this child is a positive cofactor
            probabilities[kParentIndex] += kParentNode.w * probabilities[kIndex];
        } else {
            // this child is a negative cofactor
            probabilities[kParentIndex] += probabilities[kIndex];
        }

        #ifdef DEBUG
        //if(kWriteDot)
        //    WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex, kMaxDepth);
        #endif

        Pop(i);
        if(!Empty(i))
            goto traverse;
    }
    #ifdef DEBUG
    //if(kWriteDot)
    //    WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kTrueTerminalIndex, kMaxDepth);
    #endif

    return probabilities[kRootIndex];
}

template <>
template <>
probability_t ModelCounter<ModelType::WPBDD>::Posterior<N>(){
    probability_t p, pq;

    pq = Traverse<N>(cache_,evidence_list_,condition_tier_);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["WPBDD2"];
    probs.p = 1;
    probs.pq = pq;
    #endif

    if(has_query_variable_){
        p = Traverse<N>(cache2_,evidence_list2_,condition_tier2_);
        #ifdef DEBUG
        probs.p = p;
        #endif

        if(p == 0)
            return -1;
        else return pq/p;
    } else return pq;
}



}

