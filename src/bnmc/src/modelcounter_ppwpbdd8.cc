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
#include <thread>
#include <iostream>
#include <cmath>
namespace bnmc {

using namespace std;
using namespace mstack;
using namespace bnc;

#define N 8

namespace bnmc_8_ {
    ModelCounter<ModelType::PWPBDD> *gMc;
    const Architecture *gArchitecture;
    const EvidenceList *gEvidenceList;
    const EvidenceList *gEvidenceList2;
    Cache *gCache;
    Cache *gCache2;
    ConditionTierList *gConditionTierList;
    ConditionTierList *gConditionTierList2;
    bool gHasQueryVariable;
    pthread_barrier_t start_barrier;
    pthread_barrier_t barrier;
    unsigned int gNrThreads;

    __thread unsigned int gThreadId;

}

using namespace bnmc_8_;

template <>
template <>
std::string ModelCounter<ModelType::PWPBDD>::ParallelDescription<N>(){
    return "Based on PWPBDD4 with pthreads. Bottom up traverals of partitions, using pre-initialized evidence cache. Partitions are traversed per tier, starting at the bottom tier.";
}

template<>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::ParallelGetPartitionResult<N>(
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
            throw ModelCounterException("Child component was not computed (Node (%u,%u))",kTier,kNodeId);
        #endif
        return probability;
    }
    return 1;
}

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::ParallelTraversePartition<N>(const Architecture &kArchitecture, Cache &cache, EvidenceList &evidence_list, const ConditionTierList &kConditionTierList, const unsigned int kTier, const unsigned int kNodeId){
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
    //ProbabilityList &probabilities = cache.GetProbabilityList(gThreadId);
    std::vector<probability_t> probabilities(kCircuit.size());
    probabilities[kFalseTerminalIndex] = 0;
    probabilities[kTrueTerminalIndex]  = 1;
    probabilities[kRootIndex] = 0;

    //IdentityList &identity = cache.GetIdentityList(gThreadId);
    //std::fill(identity.begin(),identity.begin()+kCircuit.size(),0);
    std::vector<unsigned int> identity(kCircuit.size(),0);

    unsigned int unique_id = 1;
    identity[kFalseTerminalIndex] = unique_id;

    StackTop i;
    StackTopInit(i);
    //Stack &s = cache.GetStack(gThreadId);
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
            const probability_t kNextTierProbability = ParallelGetPartitionResult<N>(kArchitecture, cache, evidence_list,  kTier+1);
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
probability_t ModelCounter<ModelType::PWPBDD>::ParallelTraverse<N>(const Architecture &kArchitecture, Cache &cache, const EvidenceList &kEvidenceList, const ConditionTierList &kConditionTierList, const unsigned int kTier, const unsigned int kNodeId){
    const ordering_t &kPartitionOrdering = kArchitecture.GetPartitionOrdering();
    const unsigned int kPartitionId = kPartitionOrdering[kTier];
    const Partition &kPartition = partition_[kPartitionId];
    const Partition::ArithmeticCircuit &kCircuit = kPartition.ac_;

    //ProbabilityList &probabilities = cache.GetProbabilityList(gThreadId);
    //std::fill(probabilities.begin(), probabilities.begin()+kCircuit.size(), kNotTraversed);
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
    //Stack &s = cache.GetStack(gThreadId);
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

namespace bnmc_8_ {

void pthread_TraverseTiers(
        const Architecture &kArchitecture,
        Cache              &cache,
        Cache              &cache2,
        const ConditionTierList &kConditionTierList,
        const ConditionTierList &kConditionTierList2,
        const unsigned int kThreadId,
        const unsigned int kNrThreads){

    const size_t kLastTier = kArchitecture.Size()-1;
    for(int tier = kLastTier; tier >= 0; tier--){
        const NodeIdList &kNodeIdList = cache.GetNodeIds(tier);

        for(unsigned int i = kThreadId; i < kNodeIdList.size(); i += kNrThreads){
            const unsigned int &kNodeId = kNodeIdList[i];
            EvidenceList &evidence_list = cache.GetEvidenceList(kNodeId, tier);

            //printf("[%u] (%u,%u)\n",kThreadId,tier,kNodeId);
            Probability &probability = cache.GetCacheEntry(kNodeId, tier);
            if(tier < kLastTier)
                probability = gMc->ParallelTraversePartition<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier, kNodeId);
            else probability = gMc->ParallelTraverse<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier, kNodeId);
        }

        const unsigned int kRemainder = kNodeIdList.size()%kNrThreads;
        const unsigned int kThreadId2 = (kThreadId<kRemainder?kNrThreads-kRemainder:kThreadId-kRemainder);
        const NodeIdList &kNodeIdList2 = cache2.GetNodeIds(tier);
        for(unsigned int i = kThreadId2; i < kNodeIdList2.size(); i += kNrThreads){
            const unsigned int &kNodeId = kNodeIdList2[i];
            EvidenceList &evidence_list2 = cache2.GetEvidenceList(kNodeId, tier);

            //printf("[%u] (%u,%u) (2)\n",kThreadId,tier,kNodeId);
            Probability &probability = cache2.GetCacheEntry(kNodeId, tier);
            if(tier < kLastTier)
                probability = gMc->ParallelTraversePartition<N>(kArchitecture, cache2, evidence_list2, kConditionTierList2, tier, kNodeId);
            else probability = gMc->ParallelTraverse<N>(kArchitecture, cache2, evidence_list2, kConditionTierList2, tier, kNodeId);
        }

        if(tier != 0)
            pthread_barrier_wait(&barrier);
    }
}

void pthread_TraverseTiers(
        const Architecture &kArchitecture,
        Cache              &cache,
        const ConditionTierList &kConditionTierList,
        const unsigned int kThreadId,
        const unsigned int kNrThreads){

    const size_t kLastTier = kArchitecture.Size()-1;
    for(int tier = kLastTier; tier >= 0; tier--){
        const NodeIdList &kNodeIdList = cache.GetNodeIds(tier);

        for(unsigned int i = kThreadId; i < kNodeIdList.size(); i += kNrThreads){
            const unsigned int &kNodeId = kNodeIdList[i];
            EvidenceList &evidence_list = cache.GetEvidenceList(kNodeId, tier);

            //printf("[%u] (%u,%u)\n",kThreadId,tier,kNodeId);
            Probability &probability = cache.GetCacheEntry(kNodeId, tier);
            if(tier < kLastTier)
                probability = gMc->ParallelTraversePartition<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier, kNodeId);
            else probability = gMc->ParallelTraverse<N>(kArchitecture, cache, evidence_list, kConditionTierList, tier, kNodeId);
        }

        if(tier != 0)
            pthread_barrier_wait(&barrier);
    }
}

void* pthread_workers(void *in){
    const unsigned int kThreadId = (size_t)in;
    gThreadId = kThreadId;

    pthread_barrier_wait(&start_barrier);

    if(gHasQueryVariable)
        pthread_TraverseTiers(*gArchitecture, *gCache, *gCache2, *gConditionTierList, *gConditionTierList2, kThreadId, gNrThreads);
    else pthread_TraverseTiers(*gArchitecture, *gCache, *gConditionTierList, kThreadId, gNrThreads);
}

} // namespace bnmc_8_

template <>
template <>
probability_t ModelCounter<ModelType::PWPBDD>::ParallelPosterior<N>(const unsigned int kWorkers, Timer* t){
    // globals
    gArchitecture = &architecture_;
    gMc = this;
    gEvidenceList = &evidence_list_;
    gEvidenceList = &evidence_list2_;
    gCache = &cache_;
    gCache2 = &cache2_;
    gConditionTierList = &condition_tier_;
    gConditionTierList2 = &condition_tier2_;
    gHasQueryVariable = has_query_variable_;

    //static int run = 0;
    //printf("\nrun: %u\n",run++);

    // start threads
    gNrThreads = kWorkers;
    const unsigned int kMaxThreads = std::thread::hardware_concurrency()-1;
    if(gNrThreads == 0 || gNrThreads > kMaxThreads)
        gNrThreads = kMaxThreads;

    pthread_barrier_init(&start_barrier, NULL, gNrThreads+1);
    pthread_barrier_init(&barrier, NULL, gNrThreads);
    std::vector<pthread_t> threads(gNrThreads);

    for(unsigned int i = 0; i < threads.size(); i++){
        pthread_t &thread = threads[i];
        if(pthread_create(&thread, NULL, pthread_workers, (void*)i)){
            fprintf(stderr, "Error creating thread\n");
            exit(1);
        }
    }

    pthread_barrier_wait(&start_barrier);
    if(t) t->Start();

    for(auto it = threads.begin(); it != threads.end(); it++){
        if(pthread_join(*it, NULL)){
            fprintf(stderr, "Error joining thread PPWPBDD6\n");
            exit(1);
        }
    }
    if(t) {
        t->Stop();
        t->Add();
    }

    pthread_barrier_destroy(&start_barrier);
    pthread_barrier_destroy(&barrier);

    probability_t probability_with_query_variable = cache_.GetCacheEntry(0,0);
    #ifdef DEBUG
    QueryProbabilities &probs = manager.probabilities["PPWPBDD_8_"];
    probs.pq = probability_with_query_variable;
    probs.p = 1;
    #endif

    if(has_query_variable_) {

        probability_t probability_without_query_variable = cache2_.GetCacheEntry(0,0);
        #ifdef DEBUG
        probs.p = probability_without_query_variable;
        #endif

        if(probability_without_query_variable == 0)
            return -1;
        else return probability_with_query_variable/probability_without_query_variable;
    } else return cache_.GetCacheEntry(0,0);
}

} // namespace bnmc

