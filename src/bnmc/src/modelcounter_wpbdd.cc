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

namespace bnmc {

using namespace std;
using namespace mstack;

template <>
template <>
probability_t ModelCounter<ModelType::WPBDD>::Traverse<1>(Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList){
    const unsigned int kTier = 0;
    const unsigned int kPartition = 0;
    const Partition::ArithmeticCircuit &kCircuit = partition_[kPartition].ac_;

    ProbabilityList &probabilities = cache.GetProbabilityList(0);
    std::fill(probabilities.begin(), probabilities.begin()+kCircuit.size(), kNotTraversed);
    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    #ifdef DEBUG
    bool kWriteDot = false;
    int kMaxDepth = -1;
    std::vector<unsigned int> identity(probabilities.size(),0);
    std::string filename = "debug.wpbdd.dot";
    #endif

    StackTop i;
    StackTopInit(i);
    Stack &s = cache.GetStack(0);
    Push(s,i,kRootIndex,kCircuit,kEvidenceList,kConditionTierList);

    traverse: {
        const unsigned int &kParentIndex = s[i-2];
        const unsigned int &kIndex = s[i-1];
        probability_t &probability = probabilities[kIndex];

        #ifdef DEBUG
        if(kWriteDot)
            WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex, kMaxDepth);
        #endif

        if(probability == kNotTraversed){
            #ifdef DEBUG
            if(i >= s.size()-1)
                throw ModelCounterException("Custom stack is going out of bounds");
            #endif
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
        if(kWriteDot)
            WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex, kMaxDepth);
        #endif

        Pop(i);
        if(!Empty(i))
            goto traverse;
    }
    #ifdef DEBUG
    if(kWriteDot)
        WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kTrueTerminalIndex, kMaxDepth);
    #endif

    return probabilities[kRootIndex];
}

template <>
template <>
probability_t ModelCounter<ModelType::WPBDD>::Posterior<1>(){
    probability_t p, pq;

    pq = Traverse<1>(cache_,evidence_list_,condition_tier_);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["WPBDD"];
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

