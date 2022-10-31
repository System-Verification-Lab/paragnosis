#include "architecturebound.h"
#include <algorithm>
#include <unordered_set>


ArchitectureBound::ArchitectureBound(){
    bn_ = NULL;
}

void ArchitectureBound::SetBayesnet(bayesnet *bn){
    this->bn_ = bn;
}

void ArchitectureBound::SetPartitions(const bn_partitions_t &kPartitions){
    std::set<unsigned int> partition_ids;
    for(unsigned int i = 0; i < kPartitions.size(); i++)
        partition_ids.insert(partition_ids.end(),i);

    SetPartitions(kPartitions, partition_ids);
}

template <class T>
void intersection(T &a, const T &b){
    auto first1 = a.begin();
    auto last1 = a.end();
    auto first2 = b.begin();
    auto last2 = b.end();

    while (first1!=last1 && first2!=last2){
        if(*first1<*first2)
            first1 = a.erase(first1);
        else if(*first2<*first1)
            ++first2;
        else {
            ++first1;
            ++first2;
        }
    }
    a.erase(first1,last1);
}

void ArchitectureBound::SetPartitions(const bn_partitions_t &kPartitions, const std::set<unsigned int> &kPartitionIds){

    #ifdef DEBUG
    if(bn_ == NULL)
        ArchBoundException("Bayesian network not set in ArchBoundException");
    #endif

    // store variable count from each partition (count > 1 == persistence variable)
    const unsigned int kVariables = bn_->get_nr_variables();
    total_variable_occurrences_.resize(kVariables);
    std::fill(total_variable_occurrences_.begin(), total_variable_occurrences_.end(),0);
    for(auto it = kPartitionIds.begin(); it != kPartitionIds.end(); it++){
        const unsigned int kPartitionId = *it;
        const partition_t &partition = kPartitions[kPartitionId].partition;
        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++)
            total_variable_occurrences_[*itv]++;
        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++)
            total_variable_occurrences_[*itv]++;
    }

    // determine persistence variables for each partition
    component_variables_.resize(kPartitions.size());
    for(auto it = kPartitionIds.begin(); it != kPartitionIds.end(); it++){
        const unsigned int kPartitionId = *it;
        const partition_t &partition = kPartitions[kPartitionId].partition;
        auto &persistence = component_variables_[kPartitionId];

        persistence.clear();
        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++){
            const variable_t &kVariable = *itv;
            if(total_variable_occurrences_[kVariable] > 1)
                persistence.push_back(kVariable);
        }
        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++){
            const variable_t &kVariable = *itv;
            if(total_variable_occurrences_[kVariable] > 1)
                persistence.push_back(kVariable);
        }
        //std::sort(persistence.begin(),persistence.end());
    }

    copy_.resize(total_variable_occurrences_.size());
}

size_t ArchitectureBound::ComputeBound(const ordering_t &ordering) {
    #ifdef DEBUG
    if(bn_ == NULL)
        ArchBoundException("Bayesian network not set in ArchBoundException");
    if(ordering.size() != component_variables_.size())
        ArchBoundException("Ordering size does not match with number of registered components in ArchBoundException");
    #endif
    const uint32_t *kDimension = bn_->get_states();
    std::copy(total_variable_occurrences_.begin(), total_variable_occurrences_.end(),copy_.begin());
    std::vector<unsigned int> &variable_occurrences = copy_;

    size_t size = 0;
    std::set<variable_t> context;
    for(auto component_it = ordering.begin(); component_it != ordering.end(); component_it++){
        const unsigned int &component = *component_it;
        const std::vector<variable_t> &variables = component_variables_[component];

        // compute context of component
        for(auto variable_it = variables.begin(); variable_it != variables.end(); variable_it++){
            const variable_t &variable = *variable_it;

            unsigned int &count = variable_occurrences[variable];
            #ifdef DEBUG
            if(count == 0)
                ArchBoundException("Occurrence of variable about to go to -1 in bound computation");
            #endif

            count--;
            if(count != 0)
                context.insert(variable);
            else context.erase(variable);
        }

        // compute level width
        size_t level_size = 1;
        for(auto it = context.begin(); it != context.end(); it++)
            level_size *= kDimension[*it];

        // compute size
        size += level_size;
    }
    return size;
}

