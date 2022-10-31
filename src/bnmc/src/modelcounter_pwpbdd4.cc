#include <unistd.h>
#include <algorithm>
#include <stack>
#include <bn-to-cnf/config.h>
#include <bn-to-cnf/exceptions.h>
#include <bnc/bayesgraph.h>
#include <bnc/exceptions.h>
#include <climits>
#include "modelcounter.h"
#include "io.h"
#include "options.h"
#include "exceptions.h"
#include "mstack.h"
#include "debug.h"
#include <bnc/xary.h>
#include <pthread.h>

namespace bnmc {

using namespace std;
using namespace mstack;
using namespace bnc;

#define N 4

template <>
template <>
std::string ModelCounter<ModelType::PWPBDD>::Description<N>(){
    return "Bottom up traverals of partitions, using pre-initialized evidence cache and local runtime arrays. Partitions are traversed per tier, starting at the bottom tier.";
}

template<>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::GetPartitionResult<N>(
        const Architecture &kArchitecture,
        Cache              &cache,
        const EvidenceList &kEvidenceList,
        const unsigned int kTier){

    const size_t kMaxTier = kArchitecture.Size()-1;
    if(kTier <= kMaxTier){
        const Spanning &kSpanning = kArchitecture.GetSpanning();
        const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kTier);
        const unsigned int kNodeId = XAry::GetDecimal(manager.bn, kSpanningSet, kEvidenceList);

        Probability &probability = cache.GetCacheEntry(kNodeId, kTier);
        #ifdef DEBUG
        if(probability < 0)
            throw ModelCounterException("Child component was not computed");
        #endif
        return probability;
    }
    return 1;
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::TraversePartition<N>(const Architecture &kArchitecture, Cache &cache, EvidenceList &evidence_list, const ConditionTierList &kConditionTierList, const unsigned int kTier, const unsigned int kNodeId){
    #ifdef DEBUG
    bool kWriteDot = false;
    bool kWriteDot2 = false;
    unsigned int kMaxDepth = -1;
    unsigned int depth = 1;
    std::string filename = stringf("debug.pwpbdd.%u.dot", kTier);
    #endif

    const ordering_t &kPartitionOrdering = kArchitecture.GetPartitionOrdering();
    const unsigned int kPartitionId = kPartitionOrdering[kTier];
    const Partition &kPartition = partition_[kPartitionId];
    const Partition::ArithmeticCircuit &kCircuit = kPartition.ac_;
    std::vector<probability_t> probabilities(kCircuit.size());

    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    std::vector<unsigned int> identity(kCircuit.size(),0);

    unsigned int unique_id = 1;
    identity[kFalseTerminalIndex] = unique_id;

    StackTop i;
    StackTopInit(i);
    Stack s(cache.GetStackSize(kPartitionId));
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
            const probability_t kNextTierProbability = GetPartitionResult<N>(kArchitecture, cache, evidence_list,  kTier+1);
            const Partition::Node &kParentNode = kCircuit[kParentIndex];
            probabilities[kParentIndex] += kParentNode.w * kNextTierProbability;
        } else if(kIndex != kFalseTerminalIndex){

            if(recompute || id < kParentId){
                id = (recompute?unique_id++:kParentId);

                const Partition::Node &kNode = kCircuit[kIndex];
                const Variable &kVariable = kNode.v;
                const TierId &kVariableTier = kConditionTierList[kVariable];
                const bool kIsConditioned = kVariableTier <= kTier;
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
probability_t ModelCounter<ModelType::PWPBDD>::Traverse<N>(const Architecture &kArchitecture, Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList, const unsigned int kTier, const unsigned int kNodeId){
    const ordering_t &kPartitionOrdering = kArchitecture.GetPartitionOrdering();
    const unsigned int kPartitionId = kPartitionOrdering[kTier];
    const Partition &kPartition = partition_[kPartitionId];
    const Partition::ArithmeticCircuit &kCircuit = kPartition.ac_;
    std::vector<probability_t> probabilities(kCircuit.size(),kNotTraversed);
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
    Stack s(cache.GetStackSize(kPartitionId));
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
void ModelCounter<ModelType::PWPBDD>::TraverseTiers<N>(
        const Architecture &kArchitecture,
        Cache              &cache,
        const EvidenceList &kEvidenceList,
        const ConditionTierList &kConditionTierList){


    const size_t kTiers = kArchitecture.Size();
    // do traversals bottom up
    for(int tier = kTiers-1; tier >= 0; tier--){
        const NodeIdList &kNodeIdList = cache.GetNodeIds(tier);
        for(auto it = kNodeIdList.begin(); it != kNodeIdList.end(); it++){
            const unsigned int &kNodeId = *it;
            EvidenceList &evidence_list = cache.GetEvidenceList(kNodeId, tier);

            Probability &probability = cache.GetCacheEntry(kNodeId, tier);
            if(tier < kTiers-1)
                probability = TraversePartition<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier);
            else probability = Traverse<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier);
        }
    }
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::Posterior<N>(){
    probability_t probability_without_query_variable, probability_with_query_variable;

    // NOTE 'cache' is not used.
    TraverseTiers<N>(architecture_,cache_,evidence_list_,condition_tier_);
    probability_with_query_variable = cache_.GetCacheEntry(0,0);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["PWPBDD_4_"];
    probs.p = 1;
    probs.pq = probability_with_query_variable;
    #endif
    if(has_query_variable_){
        TraverseTiers<N>(architecture_,cache2_,evidence_list2_,condition_tier2_);
        probability_without_query_variable = cache2_.GetCacheEntry(0,0);
    } else return probability_with_query_variable;

    #ifdef DEBUG
    probs.p = probability_without_query_variable;
    #endif

    if(probability_without_query_variable == 0)
        return -1;
    else return probability_with_query_variable/probability_without_query_variable;
}

}

