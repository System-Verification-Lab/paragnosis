#include "architecture.h"
#include "options.h"
#include "exceptions.h"
#include <bnc/architecturebound.h>
#include <algorithm>
#include <bnc/xary.h>
#include <cassert>
namespace bnmc {

using namespace bnc;

Architecture::Architecture(){
}

Architecture::~Architecture(){
}

bool Architecture::IsPartitioned() const {
    return comp_.HasOrdering();
}

const Composition& Architecture::GetComposition() const {
    return comp_;
}

const ordering_t& Architecture::GetOrdering() const{
    return GetPartitionOrdering();
}

const ordering_t& Architecture::GetPartitionOrdering() const{
    return comp_.GetOrdering();
}

size_t Architecture::Size() const {
    return comp_.Size()-1;
}


const Spanning& Architecture::GetSpanning() const {
    return spanning_;
}

const Persistence& Architecture::GetPersistence() const {
    return persistence_;
}

void Architecture::Init(const bn_partitions_t &kPartitions, const Evidence &kEvidence){
    assert(false && "not implemented");
    //ComputePartitionOrdering(kPartitions, kEvidence);

    //persistence_.Init(kPartitions, comp_.GetOrdering(),kEvidence);
    //spanning_.Init(persistence_);
}

void Architecture::Init(const bn_partitions_t &kPartitions){
    DeterminePartitionOrdering(kPartitions);

    persistence_.Init(kPartitions, comp_.GetOrdering());
    spanning_.Init(persistence_);
}

void Architecture::InitConditionTierList(ConditionTierList &condition_tier_list) const {
    const unsigned int kVariables = manager.bn->get_nr_variables();

    // create list of in which tier a variable is conditioned
    condition_tier_list.resize(kVariables);
    std::fill(condition_tier_list.begin(), condition_tier_list.end(),TIER_INIT_VALUE);

    const auto &kPartitionNodes = comp_.GetPartitionNodes();
    const auto &kShared = comp_.GetShared();
    for(auto it = kPartitionNodes.begin(); it != kPartitionNodes.end(); it++){
        const Composition::Node * const kNode = (*it);
        for(auto its = kShared[kNode->partition_id].begin(); its != kShared[kNode->partition_id].end(); its++)
            if(condition_tier_list[*its] > kNode->id)
                condition_tier_list[*its] = kNode->id;
    }
}

void Architecture::ApplyEvidenceToConditionTierList(ConditionTierList &condition_tier_list, const Evidence &kEvidence, const bool kNoQuery) const {
    const EvidenceVariableSet& kEvidenceVariableSet = kEvidence.GetEvidenceVariableSet();
    for(auto variable_it = kEvidenceVariableSet.begin(); variable_it != kEvidenceVariableSet.end(); variable_it++)
        if(!kNoQuery || !kEvidence.IsQueryVariable(variable_it))
            condition_tier_list[variable_it->variable] = 0;
}

ConditionTierList Architecture::GetConditionTierList(const Evidence &kEvidence, const bool kNoQuery) const {
    ConditionTierList condition_tier_list;

    InitConditionTierList(condition_tier_list);
    ApplyEvidenceToConditionTierList(condition_tier_list,kEvidence,kNoQuery);

    return std::move(condition_tier_list);
}

ConditionTierList Architecture::GetConditionTierList() const {
    ConditionTierList  condition_tier_list;

    InitConditionTierList(condition_tier_list);
    return std::move(condition_tier_list);
}

void Architecture::SetPartitionOrdering(const ordering_t &kPartitionOrdering){
    comp_.BuildOrdering(kPartitionOrdering);
}

std::vector<size_t> Architecture::GetTierWidths() const {
    const auto &kNodes = comp_.GetNodes();

    const unsigned int kNrNodes = kNodes.size();
    std::vector<size_t> tier_widths(kNrNodes,0);
    for(unsigned int i = 0; i < kNrNodes; i++){
        assert(i == kNodes[i]->id && "Node ids not in pre-order");
        tier_widths[kNodes[i]->id] = bnc::XAry::Max(manager.bn, kNodes[i]->context);
    }

    return std::move(tier_widths);
}

template<class InputIt1, class InputIt2>
bool overlap(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2){
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
            continue;
        }
        if (*first2 < *first1) {
            ++first2;
            continue;
        }
        return true;
    }
    return false;
}

template<class Input1, class Input2>
bool overlap(Input1 &input1, Input2 &input2){
    return overlap(input1.begin(), input1.end(), input2.begin(), input2.end());
}

void Architecture::DeterminePartitionOrdering(const bn_partitions_t &kPartitions){
    comp_.SetBayesnet(manager.bn);
    comp_.SetPartitions(kPartitions);
    ordering_t ordering;

    if(manager.files.is_set(file_t::COMPOSITION_ORDERING)){
        // read composition ordering
        std::string filename = manager.files.get_filename(file_t::COMPOSITION_ORDERING);
        ordering.read(filename);

        if(ordering.size() != kPartitions.size())
            throw ArchitectureException("Composition ordering in file '%s' has %u elements but there are %u partitions",
                filename.c_str(), ordering.size(), kPartitions.size());
    } else ordering = comp_.FindOrdering();

    SetPartitionOrdering(ordering);
}

const NodeIdList& Architecture::GetParentList(const TierId kTier, const NodeIndex kGroupId) const {
    assert(false && "depricated");
    return parent_lists_[kTier][kGroupId];
}

const unsigned int Architecture::GetParentGroup(const TierId kTier, const NodeIndex kNodeId) const {
    assert(false && "depricated");
    return parent_groups_[kTier][kNodeId];
}

const bool Architecture::HasGraph() const{
    assert(false && "depricated");
    return parent_groups_.size() != 0;
}

NodeIdList Architecture::CreateNodeList(const TierId kTier, const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList) const {
 assert(false && "depricated");
    const Spanning &kSpanning = spanning_;
    const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kTier);

    // create traversal counter
    XAry xary;
    xary.SetDimension(manager.bn,kSpanningSet);
    XAry::RestrictDimension restrict_dimension;
    for(unsigned int i = 0; i < kSpanningSet.size(); i++){
        const variable_t &kVariable = kSpanningSet[i];
        const bool kIsEvidence = kConditionTierList[kVariable] == 0;
        restrict_dimension.push_back(kIsEvidence);
        xary[i] = (kIsEvidence?kEvidenceList[kVariable]:0);
    }

    // store ids to be traversed, based on evidence
    NodeIdList nodes;
    do {
        const unsigned int kNodeId = xary.GetDecimal();
        nodes.push_back(kNodeId);
    } while(xary.RestrictedIncrement(restrict_dimension));

    return std::move(nodes);
}

NodeIdList Architecture::CreateParentList(const TierId kTier, const NodeIndex kNodeId,const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList) const {
    NodeIdList nodes;

    const TierId kEvidenceTier = 0;
    const TierId kParentTier = kTier-1;
    const TierId kChildTier = kTier;
    const Spanning &kSpanning = spanning_;
    const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kParentTier);
    const SpanningSet &kSpanningSetChild = kSpanning.GetSpanningSet(kChildTier);

    XAry xary;
    xary.SetDimension(manager.bn,kSpanningSet);

    XAry xary_child;
    xary_child.SetDimension(manager.bn,kSpanningSetChild);
    xary_child.SetDecimal(kNodeId);

    XAry::RestrictDimension restrict_dimension;
    restrict_dimension.resize(kSpanningSet.size());
    std::fill(restrict_dimension.begin(), restrict_dimension.end(), false);
    for(unsigned int i = 0; i < kSpanningSet.size(); i++){
        const variable_t &kVariable = kSpanningSet[i];
        if(kConditionTierList[kVariable] == kEvidenceTier){
            xary[i] = kEvidenceList[kVariable];
            restrict_dimension[i] = true;
        } else if(kConditionTierList[kVariable] <= kParentTier){
            unsigned int pos = std::distance(kSpanningSetChild.begin(),std::find(kSpanningSetChild.begin(), kSpanningSetChild.end(), kVariable ));
            if(pos < kSpanningSetChild.size()){
                xary[i] = xary_child[pos];
                restrict_dimension[i] = true;
            } else xary[i] = 0;
        } else xary[i] = 0;
    }

    // create parent list of child (kTierId,kNodeId)
    do {
        const unsigned int kParentId = xary.GetDecimal();
        nodes.push_back(kParentId);
    } while(xary.RestrictedIncrement(restrict_dimension));

    return std::move(nodes);
}


NodeIdList Architecture::CreateSiblingList(const TierId kTier, const NodeIndex kNodeId, const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList) const {
    NodeIdList nodes;
    //if(kTier >= Size()-1 || kTier == 0)
    //    return std::move(nodes);
    const TierId kEvidenceTier = 0;
    const TierId kParentTier = kTier-1;
    const TierId kChildTier = kTier;
    const Spanning &kSpanning = spanning_;
    const SpanningSet &kSpanningSetChild = kSpanning.GetSpanningSet(kChildTier);

    XAry xary;
    xary.SetDimension(manager.bn,kSpanningSetChild);
    xary.SetDecimal(kNodeId);

    XAry::RestrictDimension restrict_dimension;
    restrict_dimension.resize(kSpanningSetChild.size());
    std::fill(restrict_dimension.begin(), restrict_dimension.end(), false);
    for(unsigned int i = 0; i < kSpanningSetChild.size(); i++){
        const variable_t &kVariable = kSpanningSetChild[i];
        if(kConditionTierList[kVariable] == kEvidenceTier){
            xary[i] = kEvidenceList[kVariable];
            restrict_dimension[i] = true;
        } else if(kConditionTierList[kVariable] <= kParentTier){
            restrict_dimension[i] = true;
        }
    }

    // create sibling list node (kTierId,kNodeId)
    do {
        const unsigned int kSiblingId = xary.GetDecimal();
        nodes.push_back(kSiblingId);
    } while(xary.RestrictedIncrement(restrict_dimension));

    return std::move(nodes);
}

NodeIdList Architecture::CreateChildList(const TierId kTier, const NodeIndex kNodeId, const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList) const {
    NodeIdList nodes;
    //if(kTier >= Size()-1)
    //    return std::move(nodes);

    const TierId kEvidenceTier = 0;
    const TierId kParentTier = kTier;
    const TierId kChildTier = kTier + 1;
    const Spanning &kSpanning = spanning_;
    const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kParentTier);
    const SpanningSet &kSpanningSetChild = kSpanning.GetSpanningSet(kChildTier);

    XAry xary;
    xary.SetDimension(manager.bn,kSpanningSet);
    xary.SetDecimal(kNodeId);

    XAry xary_child;
    xary_child.SetDimension(manager.bn,kSpanningSetChild);

    XAry::RestrictDimension restrict_dimension;
    restrict_dimension.resize(kSpanningSetChild.size());
    std::fill(restrict_dimension.begin(), restrict_dimension.end(), false);
    for(unsigned int i = 0; i < kSpanningSetChild.size(); i++){
        const variable_t &kVariable = kSpanningSetChild[i];
        if(kConditionTierList[kVariable] == kEvidenceTier){
            xary_child[i] = kEvidenceList[kVariable];
            restrict_dimension[i] = true;
        } else if(kConditionTierList[kVariable] <= kParentTier){
            unsigned int pos = std::distance(kSpanningSet.begin(),std::find(kSpanningSet.begin(), kSpanningSet.end(), kVariable ));
            assert(pos < kSpanningSet.size());
            xary_child[i] = xary[pos];
            restrict_dimension[i] = true;
        } else xary_child[i] = 0;
    }

    // create child list of parent (kTierId,kNodeId)
    do {
        const unsigned int kChildId = xary_child.GetDecimal();
        nodes.push_back(kChildId);
    } while(xary_child.RestrictedIncrement(restrict_dimension));

    return std::move(nodes);
}

void Architecture::CreateGraph(){
    const Persistence &kPersistence = persistence_;
    const Spanning &kSpanning = spanning_;
    const ConditionTierList kConditionTierList = kPersistence.GetConditionTierList();
    EvidenceList evidence_list(kConditionTierList.size(),0);
    const unsigned int kTiers = Size();

    parent_lists_.clear();
    parent_groups_.clear();
    parent_lists_.resize(kTiers);
    parent_groups_.resize(kTiers);

    for(unsigned int tier = kTiers-1; tier > 0; tier--){
        const TierId kParentTier = tier-1;
        const SpanningSet &kSpanningSet = kSpanning.GetSpanningSet(kParentTier);
        const SpanningSet &kSpanningSetChild = kSpanning.GetSpanningSet(tier);

        // create traversal counter
        XAry xary;
        xary.SetDimension(manager.bn,kSpanningSetChild);
        const unsigned int kNrChildren = xary.Max();

        std::vector<bool> done(kNrChildren,false);
        parent_groups_[tier].resize(kNrChildren,-1);
        for(unsigned int node_id = 0; node_id < done.size(); node_id++){
            if(!done[node_id]){
                NodeIdList parents = CreateParentList(tier,node_id,kConditionTierList,evidence_list);
                NodeIdList siblings = CreateSiblingList(tier,node_id,kConditionTierList, evidence_list);
                #ifdef DEBUG
                if(parents.size() > 0){
                    NodeIdList children = CreateChildList(tier-1,parents[0],kConditionTierList, evidence_list);
                    assert(children == siblings);
                }
                for(unsigned int i = 0; i < parent_lists_[tier].size(); i++){
                    auto &other_parents = parent_lists_[tier][i];
                    assert(other_parents != parents);
                    NodeIdList intersection;
                    std::set_intersection(other_parents.begin(),other_parents.end(), parents.begin(), parents.end(),std::back_inserter(intersection));
                    assert(intersection.empty());
                }
                #endif
                // store parent group
                unsigned int kParentGroupId = parent_lists_[tier].size();
                parent_lists_[tier].emplace_back(parents);

                // store parent group id
                for(unsigned int i = 0; i < siblings.size(); i++){
                    #ifdef DEBUG
                    assert(!done[siblings[i]]);
                    assert(parent_groups_[tier][siblings[i]] >= kParentGroupId);
                    #endif
                    done[siblings[i]] = true;
                    parent_groups_[tier][siblings[i]] = kParentGroupId;
                }
                #ifdef DEBUG
                assert(done[node_id]);
                #endif
            }
        }
    }
}

}


