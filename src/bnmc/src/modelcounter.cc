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
#include <cassert>
using namespace std;

namespace bnmc {

#define ENCODE(x) (-(1+x))
#define DECODE(x) ENCODE(x)

template <ModelType kModelType>
ModelCounter<kModelType>::ModelCounter() :
    kNotTraversed(-4),
    kTraversed(-3),
    kRootIndex(2),
    kTrueTerminalIndex(1),
    kFalseTerminalIndex(0) {
}


template <ModelType kModelType>
void ModelCounter<kModelType>::Prepare(){

    size_t circuit_size = partition_[0].GetCircuitSize();
    assert(circuit_size != 0);

    cache_.SetSizes(bn_partitions_,partition_,architecture_);
    cache_.Resize(1);

    cache2_.SetSizes(bn_partitions_,partition_,architecture_);
    cache2_.Resize(1);

    // set condition tier list with no evidence
    const Persistence &kPersistence = architecture_.GetPersistence();
    kPersistence.FillConditionTierList(condition_tier_no_evidence_);
}

template <ModelType kModelType>
void ModelCounter<kModelType>::SetEvidence(const Evidence &kEvidence){
    assert(condition_tier_no_evidence_.size() > 0);

    query_variable_ = kEvidence.GetQueryVariable();
    //has_query_variable_ = kEvidence.HaveQueryVariable() && kEvidence.Size() != 1;
    has_query_variable_ = kEvidence.HaveQueryVariable();

    // set evidence
    auto evidence_list = kEvidence.GetEvidenceList();
    evidence_list_.resize(evidence_list.size());
    std::copy(evidence_list.begin(), evidence_list.end(),evidence_list_.begin());

    // set condition tiers
    const Persistence &kPersistence = architecture_.GetPersistence();
    condition_tier_.resize(condition_tier_no_evidence_.size());
    std::copy(condition_tier_no_evidence_.begin(),condition_tier_no_evidence_.end(), condition_tier_.begin());
    kPersistence.ApplyEvidenceToConditionTierList(condition_tier_,kEvidence);

    // set evidence and condition tiers without query variable
    if(has_query_variable_){
        evidence_list2_.resize(evidence_list.size());
        std::copy(evidence_list.begin(), evidence_list.end(),evidence_list2_.begin());

        condition_tier2_.resize(condition_tier_.size());
        std::copy(condition_tier_.begin(),condition_tier_.end(), condition_tier2_.begin());
        condition_tier2_[query_variable_] = condition_tier_no_evidence_[query_variable_];
    }
}

template <ModelType kModelType>
const bn_partitions_t& ModelCounter<kModelType>::GetBnPartitions() const {
    return bn_partitions_;
}

template <ModelType kModelType>
void ModelCounter<kModelType>::Read(){
    partition_.resize(1);
    partition_[0].Read(kModelType);

    Prepare();
}


//template <ModelType kModelType>
//void ModelCounter<kModelType>::ExcludeQueryVariable(){
//    query_variable_condition_tier_ = condition_tier_[query_variable_];
//    condition_tier_[query_variable_] = condition_tier_no_evidence_[query_variable_];
//    query_variable_value_ = evidence_list_[query_variable_];
//}
//
//template <ModelType kModelType>
//void ModelCounter<kModelType>::RestoreQueryVariable(){
//    condition_tier_[query_variable_] = query_variable_condition_tier_;
//    evidence_list_[query_variable_] = query_variable_value_;
//}

template <ModelType kModelType>
const Architecture& ModelCounter<kModelType>::GetArchitecture() const {
    return architecture_;
}



template <ModelType kModelType>
template <std::size_t N>
probability_t ModelCounter<kModelType>::Posterior(){
    probability_t p, pq;

    pq = Traverse<N>(cache_,evidence_list_);
    if(has_query_variable_){
        p = Traverse<N>(cache_,evidence_list2_);
        if(p == 0)
            return -1;
        else return pq/p;
    } else return pq;
}

// template instantiations
template class ModelCounter<ModelType::WPBDD>;
template class ModelCounter<ModelType::PWPBDD>;
template class ModelCounter<ModelType::MULTIGRAPH>;
template class ModelCounter<ModelType::TDMULTIGRAPH>;

}

