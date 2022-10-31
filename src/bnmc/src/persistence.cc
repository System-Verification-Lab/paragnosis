#include "persistence.h"
#include "options.h"
#include "exceptions.h"
#include <algorithm>

namespace bnmc {

void Persistence::Init(const bn_partitions_t &partitions, const ordering_t &kPartitionOrdering){
    Evidence evidence;
    DeterminePersistenceVariables(partitions,kPartitionOrdering,evidence);
}

void Persistence::Init(const bn_partitions_t &partitions, const ordering_t &kPartitionOrdering,const Evidence &kEvidence){
    DeterminePersistenceVariables(partitions,kPartitionOrdering,kEvidence);
}

const PersistenceSet& Persistence::GetPersistenceVariables() const {
    return persistence_variables_;
}

const PersistenceSet& Persistence::GetPersistenceSet(const unsigned int i) const {
    return persistence_sets_[i];
}

PersistenceList Persistence::GetPersistenceList() const {
    const unsigned int kVariables = manager.bn->get_nr_variables();
    PersistenceList persistence_list(kVariables,false);
    for(auto it = persistence_variables_.begin(); it != persistence_variables_.end(); it++)
        persistence_list[*it] = true;

    return std::move(persistence_list);
}

const PersistenceSets& Persistence::GetPersistenceSets() const {
    return persistence_sets_;
}

const PersistenceVariableCount& Persistence::GetPersistenceVariableCount() const {
    return persistence_variable_count_;
}

void Persistence::DeterminePersistenceVariables(const bn_partitions_t &kPartitions, const ordering_t &kPartitionOrdering,const Evidence &kEvidence){

     // store variable count from each partition (count > 1 == persistence variable)
    const unsigned int kVariables = manager.bn->get_nr_variables();
    persistence_variable_count_.resize(kVariables);
    std::fill(persistence_variable_count_.begin(), persistence_variable_count_.end(),0);
    for(auto it = kPartitionOrdering.begin(); it != kPartitionOrdering.end(); it++){
        const unsigned int kPartitionId = *it;
        const partition_t &partition = kPartitions[kPartitionId].partition;
        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++)
            persistence_variable_count_[*itv]++;
        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++)
            persistence_variable_count_[*itv]++;
    }

    // determine persistence variables for each partition
    persistence_variables_.clear();
    persistence_sets_.resize(kPartitionOrdering.size());
    for(unsigned int i = 0; i < kPartitionOrdering.size(); i++){
        const unsigned int kPartitionId = kPartitionOrdering[i];
        const partition_t &partition = kPartitions[kPartitionId].partition;
        PersistenceSet& persistence = persistence_sets_[i];
        persistence.clear();

        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++){
            const variable_t &kVariable = *itv;
            if(persistence_variable_count_[kVariable] > 1
                    && (not kEvidence.HasEvidence(kVariable)
                    || (kEvidence.HaveQueryVariable() && kVariable == kEvidence.GetQueryVariable()))){
                persistence.push_back(kVariable);
                persistence_variables_.push_back(kVariable);
            }
        }

        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++){
            const variable_t &kVariable = *itv;
            if(persistence_variable_count_[kVariable] > 1
                    && (not kEvidence.HasEvidence(kVariable)
                    || (kEvidence.HaveQueryVariable() && kVariable == kEvidence.GetQueryVariable()))){
                persistence.push_back(kVariable);
                persistence_variables_.push_back(kVariable);
            }
        }
    }

    // make list of persistent variables (sorted and unique)
    std::sort(persistence_variables_.begin(),persistence_variables_.end());
    auto end = std::unique(persistence_variables_.begin(),persistence_variables_.end());
    persistence_variables_.resize(end-persistence_variables_.begin());
}

size_t Persistence::Size() const {
    return persistence_sets_.size();
}

void Persistence::FillConditionTierList(ConditionTierList &condition_tier_list) const {
    const unsigned int kVariables = manager.bn->get_nr_variables();

    // create list of in which tier a variable is conditioned
    condition_tier_list.resize(kVariables);
    std::fill(condition_tier_list.begin(), condition_tier_list.end(),TIER_INIT_VALUE);
    for(unsigned int tier = 0; tier < persistence_sets_.size(); tier++){
        const PersistenceSet &kPersistenceSet =  persistence_sets_[tier];
        for(auto variable_it = kPersistenceSet.begin(); variable_it != kPersistenceSet.end(); variable_it++){
            TierId &condition_tier = condition_tier_list[*variable_it];
            if(condition_tier == TIER_INIT_VALUE)
                condition_tier = tier+1;
        }
    }
}

void Persistence::ApplyEvidenceToConditionTierList(ConditionTierList &condition_tier_list, const Evidence &kEvidence) const {
    const EvidenceVariableSet& kEvidenceVariableSet = kEvidence.GetEvidenceVariableSet();
    for(auto variable_it = kEvidenceVariableSet.begin(); variable_it != kEvidenceVariableSet.end(); variable_it++)
        condition_tier_list[variable_it->variable] = 0;
}

void Persistence::ApplyEvidenceToConditionTierListNoQuery(ConditionTierList &condition_tier_list, const Evidence &kEvidence) const {
    const EvidenceVariableSet& kEvidenceVariableSet = kEvidence.GetEvidenceVariableSet();
    for(auto variable_it = kEvidenceVariableSet.begin(); variable_it != kEvidenceVariableSet.end(); variable_it++)
        if(!kEvidence.IsQueryVariable(variable_it))
            condition_tier_list[variable_it->variable] = 0;
}


ConditionTierList Persistence::GetConditionTierList(const Evidence &kEvidence) const {
    ConditionTierList condition_tier_list;

    FillConditionTierList(condition_tier_list);
    ApplyEvidenceToConditionTierList(condition_tier_list,kEvidence);

    return std::move(condition_tier_list);
}

ConditionTierList Persistence::GetConditionTierListNoQuery(const Evidence &kEvidence) const {
    ConditionTierList condition_tier_list;

    FillConditionTierList(condition_tier_list);
    ApplyEvidenceToConditionTierListNoQuery(condition_tier_list,kEvidence);

    return std::move(condition_tier_list);
}


ConditionTierList Persistence::GetConditionTierList() const {
    ConditionTierList  condition_tier_list;

    FillConditionTierList(condition_tier_list);
    return std::move(condition_tier_list);
}


}
