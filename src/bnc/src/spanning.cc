#include "spanning.h"
#include <cstdlib>
#include <time.h>       /* time */
#include <queue>
#include <list>
#include <cassert>
#include <set>
#include <map>

namespace bnc {

Spanning::Spanning(){
    srand(time(NULL));
}

void Spanning::SetBayesnet(bayesnet *bn){
    bn_ = bn;
}

size_t Spanning::ComputeNrComponents(){
    const size_t kNrVariables = bn_->get_nr_variables();
    unsigned int partition_map[kNrVariables];
    std::fill(partition_map, partition_map + kNrVariables, 0);

    unsigned int id = 0;
    unsigned int nr_partitions = 0;
    unsigned int nr_variables = 0;
    for(auto it = edges_.begin(); it != edges_.end(); it++){
        const unsigned int from = it->first;
        const unsigned int to = it->second;

        if(partition_map[from] != 0){
            if(partition_map[to] != 0)
                --nr_partitions;
            else {
                partition_map[to] = partition_map[from];
                ++nr_variables;
            }
        } else if (partition_map[to] != 0) {
            partition_map[from] = partition_map[to];
            ++nr_variables;
        } else {
            ++id;
            ++nr_partitions;
            nr_variables += 2;
            partition_map[from] = id;
            partition_map[to] = id;
        }
    }

    return nr_partitions + (kNrVariables - nr_variables);
}

int Spanning::ComputePartitionMap(const bool* const kDeletedEdges, unsigned int * const partition_map, unsigned int * const partition_size, const unsigned int kNrPartitions) const {
    const size_t kNrVariables = bn_->get_nr_variables();
    std::fill(partition_map, partition_map + kNrVariables, 0);

    unsigned int id = 0;
    std::vector<unsigned int> group_map(1);
    std::vector< std::vector<unsigned int> > id_group;
    for(unsigned int i = 0; i < edges_.size(); i++){
        if(!kDeletedEdges[i]){
            const unsigned int from = edges_[i].first;
            const unsigned int to = edges_[i].second;

            if(partition_map[from] != 0){
                if(partition_map[to] != 0){
                    assert(partition_map[from] != partition_map[to] && group_map[partition_map[from]] != group_map[partition_map[to]]);
                    auto &from_group = id_group[group_map[partition_map[from]]];
                    auto &to_group = id_group[group_map[partition_map[to]]];
                    group_map[partition_map[to]] = group_map[partition_map[from]];
                    from_group.insert(from_group.end(), to_group.begin(), to_group.end());
                    to_group.clear();
                } else partition_map[to] = partition_map[from];
            } else if (partition_map[to] != 0){
                partition_map[from] = partition_map[to];
            } else {
                ++id;
                group_map.push_back(id);
                id_group.resize(id+1);
                id_group[id].push_back(id);
                partition_map[from] = id;
                partition_map[to] = id;
            }
        }
    }

    // build reindexing
    unsigned int reindex[id+1];
    id = 0;
    for(auto it = id_group.begin(); it != id_group.end(); it++){
        for(auto vit = it->begin(); vit != it->end(); vit++)
            reindex[*vit] = id;
        if(it->size() != 0)
            ++id;
    }

    // reindex partitions and compute sizes
    assert(id <= kNrPartitions);
    std::fill(partition_size, partition_size +kNrPartitions, 0);
    for(unsigned int i = 0; i < kNrVariables; i++){
        if(partition_map[i] == 0){
            partition_map[i] = id;
            ++partition_size[id];
            ++id;
        } else {
            partition_map[i] = reindex[partition_map[i]];
            ++partition_size[partition_map[i]];
        }
    }

    // return number of partitions
    return id;
}

unsigned int Spanning::GetSize() const {
    return edges_.size();
}

void Spanning::CreateSpanningTree(){
    const size_t kNrVariables = bn_->get_nr_variables();
    edges_.clear();

    bool done[kNrVariables] = {0};
    std::queue< unsigned int > q;
    unsigned int root_id = rand() % kNrVariables;
    done[root_id] = true;
    q.push(root_id);
    while(!q.empty()){
        const unsigned int kVariable = q.front();
        q.pop();

        // check neighbourhood
        const auto *kPtr = bn_->get_parent(kVariable);
        const auto *kEnd = kPtr + bn_->get_parent_size(kVariable);
        while(kPtr != kEnd){
            if(!done[*kPtr]){
                edges_.emplace_back(kVariable,*kPtr);
                done[*kPtr] = true;
                q.push(*kPtr);
            }
            ++kPtr;
        }

        kPtr = bn_->get_child(kVariable);
        kEnd = kPtr + bn_->get_child_size(kVariable);
        while(kPtr != kEnd){
            if(!done[*kPtr]){
                edges_.emplace_back(kVariable,*kPtr);
                done[*kPtr] = true;
                q.push(*kPtr);
            }
            ++kPtr;
        }
    }
}

const Spanning::EdgeList& Spanning::GetEdges() const {
    return edges_;
}

}

