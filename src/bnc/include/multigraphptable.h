#ifndef MULTIGRAPHP_TABLE_H
#define MULTIGRAPHP_TABLE_H

#include "multigraphpdef.h"
#include "hash.h"
#include <unordered_map>
#include <unordered_set>
namespace bnc {

inline Hash MultiNodeProbabilityHash(const MultiGraphProbabilityDef::Node *n){
    Hasher hasher;
    hasher.Seed(n->variable);
    assert(n->size > 1);

    const auto kEdgeEnd = &(n->edges[n->size]);
    auto edge = &(n->edges[0]);
    while(edge != kEdgeEnd){
        hasher.AddHash((size_t)edge->to);
        hasher.AddHash((size_t)edge->probability);
        ++edge;
    }

    return hasher.GetHash();
}

struct MultiGraphProbabilityHashPtr {
    Hash operator()(const MultiGraphProbabilityDef::Node *n) const {
        return MultiNodeProbabilityHash(n);
    }
};

struct MultiGraphProbabilityEqPtr {
    bool operator()(const MultiGraphProbabilityDef::Node *m, const MultiGraphProbabilityDef::Node *n) const {
        return *m == *n;
    }
};

struct MultiGraphProbabilityComputedTable : public std::unordered_set<MultiGraphProbabilityDef::Node*, MultiGraphProbabilityHashPtr, MultiGraphProbabilityEqPtr> {
    inline MultiGraphProbabilityDef::Node* find_or_insert(MultiGraphProbabilityDef::Node *n){
        auto succes = insert(n);
        if(succes.second)
            return NULL;
        return *succes.first;
    }

    inline size_t count_collisions() const {
        size_t collisions = 0;
        for (size_t bucket = 0; bucket != bucket_count(); bucket++){
            size_t size = bucket_size(bucket);
            if(size > 1){
                collisions += size-1;
            }
        }
        //for (size_t bucket = 0; bucket != bucket_count(); ++bucket)
        //    if (bucket_size(bucket) > 1)
        //        ++collisions;

        return collisions;
    }

    inline void stats() const {
        size_t oversized_buckets = 0;
        size_t active_buckets = 0;
        size_t max_bucket_entries = 0;
        size_t largest_bucket = 0;
        for (size_t bucket = 0; bucket != bucket_count(); bucket++){
            size_t size = bucket_size(bucket);
            if(size > 0){
                active_buckets++;
                if(size > 1)
                    oversized_buckets++;
                if(size > max_bucket_entries){
                    max_bucket_entries = size;
                    largest_bucket = bucket;
                }
            }
        }

        #ifdef DEBUG
        printf("\nHashmap distribution:\n");
        printf("=============\n");
        for (size_t bucket = 0; bucket != bucket_count(); bucket++){
            size_t size = bucket_size(bucket);
            printf("%3u] %4lu | ", bucket,size);
            for(auto it = begin(bucket); it != end(bucket); it++){
                MultiGraphProbabilityDef::Node *node = *it;
                printf("[%lu",node->variable);
                for(auto edge_it = node->EdgeBegin(); edge_it != node->EdgeEnd(); edge_it++)
                    printf(" -> %lu (%.3lf)", (size_t) edge_it->to, edge_it->probability);
                printf("] ");
            }
            printf("\n");
        }
        printf("=============\n");
        #endif

        printf("Entries         : %lu\n", size());
        printf("Buckets         : %lu\n", bucket_count());
        printf("Active buckets  : %lu\n", active_buckets);
        printf("Buckets > 1     : %lu\n", oversized_buckets);
        printf("Most entries/B  : %lu in bucket %lu\n", max_bucket_entries,largest_bucket);
        printf("Collisions      : %lu\n", count_collisions());
        printf("Load factor     : %lf\n", load_factor());
        printf("Active load     : %lf\n", (double)size()/(double)active_buckets);
        printf("Max load factor : %lf\n", max_load_factor());

        printf("\n");
    }
};

}

#endif

