#include <unistd.h>
#include <algorithm>
#include <stack>
#include <thread>
#include <bn-to-cnf/config.h>
#include <bn-to-cnf/exceptions.h>
#include <bnc/bayesgraph.h>
#include <bnc/exceptions.h>
#include <bnc/xary.h>
#include <climits>
#include "modelcounter.h"
#include "io.h"
#include "options.h"
#include "exceptions.h"
#include "mstack.h"
#include "debug.h"

namespace bnmc {


using namespace std;
using namespace mstack;
using namespace bnc;

#define N 6


template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::TraverseArchitecture<N>(const Composition::Node * const kCompositionNode, Cache &cache, EvidenceList &evidence_list, const ConditionTierList &kConditionTierList);

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::TraversePartition<N>(const Composition::Node * const kCompositionNode, Cache &cache, EvidenceList &evidence_list, const ConditionTierList &kConditionTierList, const unsigned int kNodeId){
    const unsigned int kTier = kCompositionNode->id;
    const unsigned int kPartition = kCompositionNode->partition_id;
    const Partition::ArithmeticCircuit &kCircuit = partition_[kPartition].ac_;

    #ifdef DEBUG
    bool kWriteDot = false;
    bool kWriteDot2 = false;
    unsigned int kMaxDepth = -1;
    unsigned int depth = 1;
    std::string filename = stringf("debug.pwpbdd.%u.dot", kTier);
    #endif

    ProbabilityList &probabilities = cache.GetProbabilityList(kTier);
    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    IdentityList &identity = cache.GetIdentityList(kTier); // identity[node] = [0...kMaxTier], implementation of IsTraversed(node)
    std::fill(identity.begin(),identity.begin()+kCircuit.size(),0);

    unsigned int unique_id = 1;
    identity[kFalseTerminalIndex] = unique_id;

    StackTop i;
    StackTopInit(i);
    Stack &s = cache.GetStack(kTier);
    Push(s,i,kFalseTerminalIndex,kRootIndex);

    bool recompute = false;
    traverse:
    const unsigned int &kParentIndex = s[i-2];
    const unsigned int &kIndex = s[i-1];
    {
        const unsigned int &kParentId = identity[kParentIndex];
        unsigned int &id = identity[kIndex];

        #ifdef DEBUG
        if(kWriteDot)
            WriteDot(filename, kCircuit, evidence_list, probabilities, identity, kTier, kIndex,kMaxDepth,i,&s);
        #endif

        if(kIndex == kTrueTerminalIndex){
            // this is a positive cofactor (implicitly by language definition)
            // true terminal is linked to root of next partition

            const probability_t kNextTierProbability = TraverseArchitecture<N>(kCompositionNode, cache, evidence_list, kConditionTierList);
            const Partition::Node &kParentNode = kCircuit[kParentIndex];
            probabilities[kParentIndex] += kParentNode.w * kNextTierProbability;
        } else if(kIndex != kFalseTerminalIndex){

            if(recompute || id < kParentId){
                id = (recompute?unique_id++:kParentId);

                const Partition::Node &kNode = kCircuit[kIndex];
                const Variable &kVariable = kNode.v;
                const TierId &kVariableTier = kConditionTierList[kVariable];
                const bool kIsConditioned = kVariableTier < kTier;
                VariableValue &value = evidence_list[kVariable];
                const bool kIsTrue = value == kNode.i;
                recompute = !kIsConditioned && is_persistence_[kVariable];

                if(!kIsConditioned){
                    // add both cofactors
                    Push(s,i,kIndex,kNode.e);
                    Push(s,i,kIndex,kNode.t);

                    // record decision i.o.t. to the same in the following partition(s)
                    // this is only important for persistence variables
                    value = kNode.i;
                } else if(kIsTrue) {
                    // add only positive cofactor
                    Push(s,i,kIndex,kNode.t);
                } else {
                    // add only negative cofactor
                    Push(s,i,kIndex,kNode.e);
                }

                // initialize probability
                probabilities[kIndex] = 0;
                goto traverse;

            } else {
                // this subgraph has been traversed
                // a child node computes the probability of the parent node
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
                    WriteDot(filename, kCircuit, evidence_list, probabilities, identity, kTier, kIndex,kMaxDepth,i,&s);
                #endif
            }
        }
    }

    recompute = false;
    Pop(i);
    // stop when encountering root node
    if(i > kStackItemSize) // root element should remain on stack
        goto traverse;

    return probabilities[kRootIndex];
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::Traverse<N>(const Composition::Node * const kCompositionNode, Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList, const unsigned int kNodeId){
    const unsigned int kTier = kCompositionNode->id;
    const unsigned int kPartition = kCompositionNode->partition_id;
    const Partition::ArithmeticCircuit &kCircuit = partition_[kPartition].ac_;

    ProbabilityList &probabilities = cache.GetProbabilityList(kTier);
    std::fill(probabilities.begin(), probabilities.begin()+kCircuit.size(), kNotTraversed);
    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    #ifdef DEBUG
    bool kWriteDot = false;
    unsigned int kMaxDept = -1;
    std::vector<unsigned int> identity(probabilities.size(),0);
    std::string filename = stringf("debug.pwpbdd.%u.dot", kTier);
    #endif

    StackTop i;
    StackTopInit(i);
    Stack &s = cache.GetStack(kTier);
    Push(s,i,kRootIndex,kCircuit,kEvidenceList,kConditionTierList,kTier);

    traverse: {
        const unsigned int &kParentIndex = s[i-2];
        const unsigned int &kIndex = s[i-1];
        probability_t &probability = probabilities[kIndex];

        #ifdef DEBUG
        if(kWriteDot)
            WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex,kMaxDept,i,&s);
        #endif

        if(probability == kNotTraversed){
            #ifdef DEBUG
            if(i >= s.size()-1){
                WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex,kMaxDept,i,&s);
                throw ModelCounterException("Custom stack is going out of bounds");
            }
            #endif
            probability = 0;
            Push(s,i,kIndex,kCircuit,kEvidenceList,kConditionTierList,kTier);

            goto traverse;
        }

        const Partition::Node &kParentNode = kCircuit[kParentIndex];
        const Partition::Node &kNode = kCircuit[kIndex];
        if(kParentNode.v != kNode.v){
            // this is a t child
            probabilities[kParentIndex] += kParentNode.w * probability;
        } else {
            // this is an e child
            probabilities[kParentIndex] += probability;
        }

        #ifdef DEBUG
        if(kWriteDot)
            WriteDot(filename, kCircuit, kEvidenceList, probabilities, identity, kTier, kIndex,kMaxDept,i,&s);
        #endif

        Pop(i);
        if(!Empty(i))
            goto traverse;
    }

    return probabilities[kRootIndex];
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::TraverseArchitecture<N>(
        const Composition::Node * const kCompositionNode,
        Cache              &cache,
        EvidenceList       &evidence_list,
        const ConditionTierList &kConditionTierList){

    probability_t mar_probability = 1;
    for(auto it = kCompositionNode->children.begin(); it != kCompositionNode->children.end(); it++){
        const Composition::Node * const kChildNode = *it;
        const unsigned int kNodeId = XAry::GetDecimal(manager.bn, kChildNode->context, evidence_list);
        const unsigned int kTier = kChildNode->id;

        Probability &probability = cache.GetCacheEntry(kNodeId, kTier);
        if(probability == Cache::kNotCached){
            if(!Composition::IsLeaf(kChildNode))
                probability = TraversePartition<N>(kChildNode, cache, evidence_list, kConditionTierList);
            else probability = Traverse<N>(kChildNode, cache, evidence_list, kConditionTierList);
        }
        mar_probability *= probability;
    }
    return mar_probability;
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::Posterior<N>(){
    probability_t p, pq;

    const Composition::Node * const kCompositionRoot = architecture_.GetComposition().GetCompositionOrdering();


    pq = TraverseArchitecture<N>(kCompositionRoot,cache_, evidence_list_, condition_tier_);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["CPWPBDD"];
    probs.pq = pq;
    probs.p = 1;
    #endif
    if(has_query_variable_){
        //cache_.Reset();
        p = TraverseArchitecture<N>(kCompositionRoot,cache2_, evidence_list2_, condition_tier2_);
        #ifdef DEBUG
        probs.p = p;
        #endif

        if(p == 0)
            return -1;
        else return pq/p;
    } else return pq;
}

}

