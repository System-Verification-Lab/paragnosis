#include "cache.h"
#include "exceptions.h"
#include "options.h"
#include <algorithm>
#include <cmath>
#include <bnc/xary.h>
#include <cassert>

namespace bnmc {

using namespace bnc;
using namespace mstack;

//const CacheValue Cache::kNotCached = -1;
const probability_t Cache::kNotCached = -1;
const probability_t Cache::kTraversed = -2;
const unsigned int Cache::kMaxBufferSize = 128;
Cache::Cache(){
    buffer_size_ = 0;
    max_circuit_size_ = 0;
    max_stack_size_ = 0;
    max_alloc_size_ = 0;
}

void Cache::SetTierSizes(const Architecture &kArchitecture){
    if(kArchitecture.IsPartitioned()){
        const auto &kComposition = kArchitecture.GetComposition();
        assert(kComposition.HasOrdering() && "Composition does not contain an ordering");
        const auto &kNodes = kComposition.GetNodes();
        tier_sizes_.resize(kNodes.size());
        for(auto it = kNodes.begin(); it != kNodes.end(); it++){
            const Composition::Node * const kNode = *it;
            tier_sizes_[kNode->id] = XAry::Max(manager.bn, kNode->context);
        }
    }
}

void Cache::SetCircuitSizes(const std::vector< Partition > &kPartitions){
    max_circuit_size_ = 0;
    circuit_sizes_.clear();
    for(auto it = kPartitions.begin(); it != kPartitions.end(); it++){
        size_t size = it->GetCircuitSize();
        if(size > max_circuit_size_)
            max_circuit_size_ = size;
        circuit_sizes_.push_back(size);
    }
}

void Cache::SetStackSizes(const bn_partitions_t &kBnPartitions){
    stack_sizes_.clear();
    if(kBnPartitions.size() == 0){
        const unsigned int kVariables = manager.bn->get_nr_variables();
        max_stack_size_ = 0;
        for(unsigned int i = 0; i < kVariables; i++)
            max_stack_size_ += manager.mapping.get_dimension()[i];

        max_stack_size_ *= 4;
        stack_sizes_.push_back(max_stack_size_);
    } else {
        max_stack_size_ = 0;
        for(unsigned int i = 0; i < kBnPartitions.size(); i++){
            size_t nr_literals = 0;
            const partition_t &kPartition = kBnPartitions[i].partition;
            for(auto variable_it = kPartition.set.begin(); variable_it != kPartition.set.end(); variable_it++)
                nr_literals += manager.mapping.get_dimension()[*variable_it];
            for(auto variable_it = kPartition.cutset.begin(); variable_it != kPartition.cutset.end(); variable_it++)
                nr_literals += manager.mapping.get_dimension()[*variable_it];

            size_t size = nr_literals*4;
            stack_sizes_.push_back(size);
            if(size > max_stack_size_)
               max_stack_size_ = size;
        }
    }
}

const size_t Cache::GetMaxAllocSize() const{
    return max_alloc_size_;
}

void Cache::SetAllocSizes(){
    assert(circuit_sizes_.size() == stack_sizes_.size());
    max_alloc_size_ = 0;
    alloc_sizes_.clear();
    for(unsigned int partition_id = 0; partition_id < circuit_sizes_.size(); partition_id++){
        const size_t &kCircuitSize = circuit_sizes_[partition_id];
        const size_t kProbabilitySize = sizeof(probability_t)*kCircuitSize;
        const size_t kIdentitySize = sizeof(NodeIndex)*kCircuitSize;
        const size_t kStackSize = sizeof(NodeIndex)*stack_sizes_[partition_id];
        const size_t kBytes = kProbabilitySize + kIdentitySize + kStackSize;

        if(kBytes > max_alloc_size_)
            max_alloc_size_ = kBytes;

        alloc_sizes_.push_back(kBytes);
    }
}

void Cache::SetSizes(const bn_partitions_t &kBnPartitions,const std::vector< Partition > &kPartitions, const Architecture &kArchitecture){
    SetTierSizes(kArchitecture);
    SetCircuitSizes(kPartitions);
    SetStackSizes(kBnPartitions);
    SetAllocSizes();
}

void Cache::Prepare(){
    PrepareTierProbabilities();
}

void Cache::InitTierProbabilitiers(){
    for(auto tier_it = cache_.begin(); tier_it != cache_.end(); tier_it++)
        std::fill(tier_it->begin(),tier_it->end(),kNotCached);
}

void Cache::Init(){
    InitTierProbabilitiers();
}

void Cache::PrepareTierProbabilities(){
    assert(tier_sizes_.size() != 0);
    cache_.resize(tier_sizes_.size());
    for(unsigned int tier = 0; tier < tier_sizes_.size(); tier++)
        cache_[tier] = TierCache(tier_sizes_[tier]);
}

void Cache::PrepareEvidenceLists(){
    assert(tier_sizes_.size() != 0);
    const unsigned int kTiers = tier_sizes_.size();
    const unsigned int kVariables = manager.bn->get_nr_variables();

    evidence_cache_.resize(kTiers);
    for(unsigned int tier = 0; tier < kTiers; tier++){
        evidence_cache_[tier].resize(tier_sizes_[tier]);
        for(auto it = evidence_cache_[tier].begin(); it != evidence_cache_[tier].end(); it++)
            it->resize(kVariables);
    }
}

void Cache::Resize(const unsigned int kBufferSize){
    assert(kBufferSize <= kMaxBufferSize);
    assert(max_circuit_size_ != 0 && max_stack_size_ != 0);
    buffer_size_ = kBufferSize;

    // allocate probability lists
    probability_cache_.resize(kBufferSize);
    for(auto it = probability_cache_.begin(); it != probability_cache_.end(); it++)
        it->resize(max_circuit_size_);

    // allocate identity lists
    identity_cache_.resize(kBufferSize);
    for(auto it = identity_cache_.begin(); it != identity_cache_.end(); it++)
        it->resize(max_circuit_size_);

    // allocate stacks
    stack_cache_.resize(kBufferSize);
    for(auto it = stack_cache_.begin(); it != stack_cache_.end(); it++)
        it->resize(max_stack_size_);

    // fill the stacks
    stack_probabilities_.clear();
    stack_identity_.clear();
    stack_stack_.clear();
    for(int i = kBufferSize -1; i >= 0; i--){
        stack_probabilities_.push_back(i);
        stack_identity_.push_back(i);
        stack_stack_.push_back(i);
    }
}

IdentityList& Cache::GetIdentityList(const unsigned int kId) {
    return identity_cache_[kId];
}

ProbabilityList& Cache::GetProbabilityList(const unsigned int kId)  {
    return probability_cache_[kId];
}

Stack& Cache::GetStack(const unsigned int kId) {
    return stack_cache_[kId];
}

const unsigned int Cache::ObtainProbabilityListId(){
    #ifdef DEBUG
    assert(!stack_probabilities_.empty());
    #endif
    unsigned int id = stack_probabilities_.back();
    stack_probabilities_.pop_back();
    return id;
}

void Cache::ReleaseProbabilityList(const unsigned int kId){
    stack_probabilities_.push_back(kId);
}

const unsigned int Cache::ObtainStackId(){
    #ifdef DEBUG
    assert(!stack_stack_.empty());
    #endif
    unsigned int id = stack_stack_.back();
    stack_stack_.pop_back();
    return id;
}

void Cache::ReleaseStack(const unsigned int kId){
    stack_stack_.push_back(kId);
}

const unsigned int Cache::ObtainIdentityListId(){
    #ifdef DEBUG
    assert(!stack_identity_.empty());
    #endif
    unsigned int id = stack_identity_.back();
    stack_identity_.pop_back();
    return id;
}

void Cache::ReleaseIdentityList(const unsigned int kId){
    stack_identity_.push_back(kId);
}

const NodeIdList& Cache::GetNodeIds(const unsigned int kTier) const{
    return id_cache_[kTier];
};

void Cache::SetNodeIds(const Architecture &kArchitecture, const Evidence &kEvidence, bool with_query_variable){
    const Persistence &kPersistence = kArchitecture.GetPersistence();
    const Spanning &kSpanning = kArchitecture.GetSpanning();
    const ConditionTierList kConditionTierList = kPersistence.GetConditionTierList(kEvidence);
    const EvidenceList kEvidenceList = kEvidence.GetEvidenceList();
    const unsigned int kTiers = kArchitecture.Size();

    id_cache_.resize(kTiers);
    for(unsigned int tier = 0; tier < kTiers; tier++){
        const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(tier);

        // create traversal counter
        XAry xary;
        xary.SetDimension(manager.bn,kSpanningSet);

        // restrict counter based on evidence
        variable_t query_variable = 9999999;
        if(!with_query_variable && kEvidence.HaveQueryVariable())
            query_variable = kEvidence.GetQueryVariable();

        xary.SetZero();
        XAry::RestrictDimension restrict_dimension;
        for(unsigned int i = 0; i < kSpanningSet.size(); i++){
            const variable_t &kVariable = kSpanningSet[i];
            bool is_conditioned = kConditionTierList[kVariable] == 0;
            if(!with_query_variable && kEvidence.HaveQueryVariable() && kEvidence.GetQueryVariable() == kVariable)
                is_conditioned = false;

            restrict_dimension.push_back(is_conditioned);
            xary[i] = (is_conditioned?kEvidenceList[kVariable]:0);
        }

        // create input evidence
        id_cache_[tier].resize(0);
        do {
            // store ids to be traversed, based on evidence
            const unsigned int kNodeId = xary.GetDecimal();
            id_cache_[tier].push_back(kNodeId);
        } while(xary.RestrictedIncrement(restrict_dimension));
    }
}

void Cache::SetEvidence(const Architecture &kArchitecture, const Evidence &kEvidence, bool with_query_variable){
    const Persistence &kPersistence = kArchitecture.GetPersistence();
    const Spanning &kSpanning = kArchitecture.GetSpanning();
    const ConditionTierList kConditionTierList = kPersistence.GetConditionTierList(kEvidence);
    const EvidenceList kEvidenceList = kEvidence.GetEvidenceList();
    const unsigned int kTiers = kArchitecture.Size();

    evidence_cache_.resize(kTiers);
    id_cache_.resize(kTiers);
    for(unsigned int tier = 0; tier < kTiers; tier++){
        const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(tier);

        // create traversal counter
        XAry xary;
        xary.SetDimension(manager.bn,kSpanningSet);
        const unsigned int kTierWidth = xary.Max();
        evidence_cache_[tier].resize(kTierWidth);

        // restrict counter based on evidence
        variable_t query_variable = 9999999;
        if(!with_query_variable && kEvidence.HaveQueryVariable())
            query_variable = kEvidence.GetQueryVariable();

        xary.SetZero();
        XAry::RestrictDimension restrict_dimension;
        for(unsigned int i = 0; i < kSpanningSet.size(); i++){
            const variable_t &kVariable = kSpanningSet[i];
            bool is_conditioned = kConditionTierList[kVariable] == 0;
            if(!with_query_variable && kEvidence.HaveQueryVariable() && kEvidence.GetQueryVariable() == kVariable)
                is_conditioned = false;

            restrict_dimension.push_back(is_conditioned);
            xary[i] = (is_conditioned?kEvidenceList[kVariable]:0);
        }

        // create input evidence
        id_cache_[tier].resize(0);
        do {
            // store ids to be traversed, based on evidence
            const unsigned int kNodeId = xary.GetDecimal();
            id_cache_[tier].push_back(kNodeId);

            // prepare evidence
            EvidenceList &evidence_list = evidence_cache_[tier][kNodeId];
            evidence_list.resize(kEvidenceList.size());
            std::copy(kEvidenceList.begin(),kEvidenceList.end(),evidence_list.begin());
            for(unsigned int i = 0; i < kSpanningSet.size(); i++)
                evidence_list[kSpanningSet[i]] = xary[i];
        } while(xary.RestrictedIncrement(restrict_dimension));
    }
}

const TaskGroup& Cache::GetTaskGroup(const unsigned int kTaskId) const {
    return task_cache_[kTaskId];
}

const unsigned int Cache::GetTaskGroupId(const TierId kTier, const NodeIndex kNodeId) const {
    return task_id_cache_[kTier][kNodeId];
}

TaskCounter& Cache::GetTaskCounter(const unsigned int kTaskId){
    return task_dependency_cache_[kTaskId];
}


size_t Cache::GetTaskCount() const {
    return task_cache_.size();
}

void Cache::SetTasks(const Architecture &kArchitecture, const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList){
    const Spanning &kSpanning = kArchitecture.GetSpanning();
    const unsigned int kTiers = kArchitecture.Size();

    task_cache_.resize(0);
    task_reverse_cache_.resize(0);
    task_id_cache_.resize(kTiers);
    task_dependency.resize(0);
    for(unsigned int tier = kTiers-1; tier > 0; tier--){
        const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(tier);
        const unsigned int kMaxTasks = XAry::Max(manager.bn,kSpanningSet);
        std::vector<bool> done(kMaxTasks,false);


        task_id_cache_[tier].resize(kMaxTasks);
        #ifdef DEBUG
        std::fill(task_id_cache_[tier].begin(),task_id_cache_[tier].end(),-1);
        #endif

        TaskGroup task_group;
        task_group.tier = tier;
        task_group.nodes = kArchitecture.CreateNodeList(tier, kConditionTierList, kEvidenceList);
        if(tier == kTiers-1){
            task_cache_.push_back(task_group);
            task_reverse_cache_.push_back({tier,});
            task_dependency.push_back(0);
        }

        #ifdef DEBUG
        unsigned int kTaskGroupIdStart = task_cache_.size();
        #endif
        for(auto it = task_group.nodes.begin(); it != task_group.nodes.end(); it++){
            const NodeIndex &kNodeId = *it;
            if(!done[kNodeId]){
                unsigned int kTaskGroupId = task_cache_.size();
                task_cache_.resize(kTaskGroupId+1);
                task_reverse_cache_.resize(kTaskGroupId+1);
                TaskGroup &tasks = task_cache_[kTaskGroupId];
                TaskGroup &rtasks = task_reverse_cache_[kTaskGroupId];
                tasks.tier = tier - 1;
                tasks.nodes = kArchitecture.CreateParentList(tier,kNodeId,kConditionTierList,kEvidenceList);

                rtasks.tier = tier;
                rtasks.nodes = kArchitecture.CreateSiblingList(tier,kNodeId,kConditionTierList, kEvidenceList);
                NodeIdList &siblings = rtasks.nodes;
                task_dependency.push_back(siblings.size());

                #ifdef DEBUG
                if(tasks.nodes.size() > 0){
                    NodeIdList children = kArchitecture.CreateChildList(tier-1,tasks.nodes[0],kConditionTierList, kEvidenceList);
                    assert(children == siblings);
                }
                for(unsigned int i = kTaskGroupIdStart; i < task_cache_.size()-1; i++){
                    auto &other_tasks = task_cache_[i];
                    assert(!(other_tasks.nodes == tasks.nodes && other_tasks.tier == tasks.tier));
                    NodeIdList intersection;
                    std::set_intersection(other_tasks.nodes.begin(),other_tasks.nodes.end(), tasks.nodes.begin(), tasks.nodes.end(),std::back_inserter(intersection));
                    assert(intersection.empty());
                }
                #endif

                // store parent group id
                for(unsigned int i = 0; i < siblings.size(); i++){
                    #ifdef DEBUG
                    assert(!done[siblings[i]]);
                    #endif
                    done[siblings[i]] = true;
                    task_id_cache_[tier][siblings[i]] = kTaskGroupId;
                }
                #ifdef DEBUG
                assert(done[kNodeId]);
                #endif
            }
        }
    }
    // add counter for node (0,0)
    task_id_cache_[0].resize(1);
    task_id_cache_[0][0] = 0;

    // allocate dependency counter
    if(task_dependency_cache_.size() < task_cache_.size())
        task_dependency_cache_ = TaskDependencyCache(task_cache_.size());

    if(task_issue_cache_.size() != task_cache_.size())
        task_issue_cache_ = TaskDependencyCache(task_cache_.size());

    for(auto it = task_issue_cache_.begin(); it != task_issue_cache_.end(); it++)
        *it = 0;

    // set dependency counters
    auto it_to = task_dependency_cache_.begin();
    auto it_from = task_dependency.begin();
    while(it_from != task_dependency.end()){
        *it_to = *it_from;
        it_to++;
        it_from++;
    }

    #ifdef DEBUG
    std::map<NodeIndex, std::vector<NodeIndex> > dependents;
    for(int tier = kTiers-1; tier >= 0; tier--){
        // create list of nodes that must be covered
        NodeIdList nodes = kArchitecture.CreateNodeList(tier, kConditionTierList, kEvidenceList);
        unsigned int cover_size = nodes.size();
        std::map<NodeIndex,bool> cover;

        for(auto it = nodes.begin(); it != nodes.end(); it++)
            cover[*it] = false;


        auto begin = task_cache_.begin();
        while(begin != task_cache_.end() && begin->tier != tier) begin++;
        auto end = begin;
        while(end != task_cache_.end() && end->tier == tier) end++;

        // traverse tasks on this tier and check cover
        for(auto it = begin; it != end; it++){

            for(auto vit = it->nodes.begin(); vit != it->nodes.end(); vit++){
                assert(!cover[*vit] && "node was already covered");
                const unsigned int kGroupId = task_id_cache_[tier][*vit];
                if(kGroupId != 0)
                    dependents[kGroupId].push_back(*vit);
                cover[*vit] = true;
                cover_size--;
            }
        }
        assert(cover_size == 0 && "not all nodes are covered");
        NodeIdList nodes2;
        for(auto it = cover.begin(); it != cover.end(); it++)
            nodes2.push_back(it->first);
        assert(nodes == nodes2 && "not all nodes are covered");

        // traverse tasks on this tier and dependencies
        for(auto it = nodes.begin(); it != nodes.end(); it++){
            unsigned int group_id = task_id_cache_[tier][*it];
            assert(task_dependency_cache_[group_id] == dependents[group_id].size());
        }
    }
    #endif
}



void Cache::ShallowCopy(const Cache &kCache) {
    for(auto it = kCache.cache_.begin(); it != kCache.cache_.end(); it++)
        cache_.push_back(TierCache(it->size()));

    evidence_cache_.resize(kCache.evidence_cache_.size());
    for(unsigned int i = 0; i < kCache.evidence_cache_.size(); i++){
        evidence_cache_[i].resize(kCache.evidence_cache_[i].size());
        for(unsigned int j = 0; j < kCache.evidence_cache_[i].size(); j++)
            evidence_cache_[i][j].resize(kCache.evidence_cache_[i][j].size());
    }

    stack_sizes_      = kCache.stack_sizes_;
    circuit_sizes_    = kCache.circuit_sizes_;
    tier_sizes_       = kCache.tier_sizes_;
    max_circuit_size_ = kCache.max_circuit_size_;
    max_stack_size_   = kCache.max_stack_size_;
}

EvidenceList& Cache::GetEvidenceList(const unsigned int kNodeId, const unsigned int kTier){
    return evidence_cache_[kTier][kNodeId];
}

void Cache::SetEvidenceList(EvidenceList &dest, const Architecture &kArchitecture, const unsigned int &kTier, const unsigned int &kNodeId) const {
    const Spanning &kSpanning = kArchitecture.GetSpanning();
    const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kTier);

    XAry xary;
    xary.SetDimension(manager.bn,kSpanningSet);
    xary.SetDecimal(kNodeId);
    for(unsigned int i = 0; i < kSpanningSet.size(); i++)
        dest[kSpanningSet[i]] = xary[i];
}

void Cache::SetEvidenceList(EvidenceList &dest, const EvidenceList &source, const Architecture &kArchitecture, const unsigned int &kTier, const unsigned int &kNodeId) const {
    dest = source;
    SetEvidenceList(dest, kArchitecture,kTier,kNodeId);
}

Probability& Cache::GetCacheEntry(const unsigned int kNodeId, const unsigned int kTier){
    return cache_[kTier][kNodeId];
}

unsigned int Cache::GetNodeId(const Spanning &kSpanning, const EvidenceList &kEvidenceList, const unsigned int kTier) const{
    const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kTier);
    return XAry::GetDecimal(manager.bn, kSpanningSet, kEvidenceList);
}

Probability& Cache::GetCacheEntry(const Spanning &kSpanning, const EvidenceList &kEvidenceList, const unsigned int kTier){
    return GetCacheEntry(GetNodeId(kSpanning,kEvidenceList,kTier), kTier);
}

const size_t Cache::GetStackSize(const unsigned int kPartitionId) const {
    return stack_sizes_[kPartitionId];
}

const size_t Cache::GetTierSize(const unsigned int kTier) const {
    return tier_sizes_[kTier];
}

const size_t Cache::GetCircuitSize(const unsigned int kPartitionId) const {
    return circuit_sizes_[kPartitionId];
}

const size_t Cache::GetArchitectureSize() const {
    size_t size = 0;
    for(auto it = tier_sizes_.begin(); it != tier_sizes_.end(); it++)
        size += *it;

    return size;
}


}

