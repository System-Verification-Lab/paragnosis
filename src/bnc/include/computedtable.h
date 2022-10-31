#ifndef COMPUTEDTABLE_H
#define COMPUTEDTABLE_H

// included in node.h

#include "defines.h"

#define USE_UNORDERED_MAP
#ifdef USE_UNORDERED_MAP
typedef std::pair<weighted_node*, weighted_node*> computed_table_key_t;
struct computed_table_hash {
    inline size_t compute_hash(weighted_node* const & n, weighted_node* const & m) const {
        //const uintptr_t N = ((uintptr_t) n)>>8; // divide by 8
        //const uintptr_t M = ((uintptr_t) m)>>8; // divide by 8

        //const unsigned int wordsize = sizeof(void*)*8;
        //const unsigned int halfwordsize = sizeof(void*)*4;
        //return ((((size_t) n)/wordsize) << halfwordsize) | ((((size_t) m)/wordsize) & (~(size_t)0 >> halfwordsize));

        //return std::hash<size_t>()((uintptr_t) n + (uintptr_t) m);
        return std::hash<size_t>()((uintptr_t) n + (((uintptr_t) m)>>(sizeof(void*)/4)));
    }
    size_t operator()(const computed_table_key_t key) const {
            #ifdef DEBUG
            if(((uintptr_t) key.first)%4 != 0)
                fprintf(stderr, "first node adress %p (%lu) not dividable by 4\n", key.first, ((uintptr_t) key.first));
            if(((uintptr_t) key.second)%4 != 0)
                fprintf(stderr, "second node adress %p (%lu) not dividable by 4\n", key.second, ((uintptr_t) key.second));
            if(((uintptr_t) key.first)%8 != 0)
                fprintf(stderr, "first node adress %p (%lu) not dividable by 8\n", key.first, ((uintptr_t) key.first));
            if(((uintptr_t) key.second)%8 != 0)
                fprintf(stderr, "second node adress %p (%lu) not dividable by 8\n", key.second, ((uintptr_t) key.second));

            #endif
            #if __SIZEOF_POINTER__ == __SIZEOF_SIZE_T__
            #define UINTPTR_T_IS_SIZE_T
            #endif
            return compute_hash(key.first,key.second);
    }
};

struct computed_table_equal {
    bool operator()(const computed_table_key_t key1, const computed_table_key_t key2) const {
        return key1 == key2;
    }
};

typedef std::unordered_map<computed_table_key_t, weighted_node*, computed_table_hash, computed_table_equal> computed_table_base_t;
#else
typedef std::pair<weighted_node*, weighted_node*> computed_table_key_t;
typedef std::map< computed_table_key_t, weighted_node**> computed_table_base_t;
#endif

struct weighted_node_computed_table : public computed_table_base_t {
    inline size_t count_collisions() const {
        size_t collisions = 0;
        for (size_t bucket = 0; bucket != bucket_count(); bucket++){
            size_t size = bucket_size(bucket);
            if(size > 1){
                collisions += size-1;
            }
        }
        for (size_t bucket = 0; bucket != bucket_count(); ++bucket)
            if (bucket_size(bucket) > 1)
                ++collisions;

        return collisions;
    }

    inline void stats() const {
        size_t collisions = 0;
        size_t oversized_buckets = 0;
        size_t active_buckets = 0;
        for (size_t bucket = 0; bucket != bucket_count(); bucket++){
            size_t size = bucket_size(bucket);
            if(size > 0){
                active_buckets++;
                collisions += size-1;
                if(size > 1)
                    oversized_buckets++;
            }
        }

        //#ifdef DEBUG
        //printf("\nHashmap distribution:\n");
        //printf("=============\n");
        //computed_table_hash hasher;
        //for (size_t bucket = 0; bucket != bucket_count(); bucket++){
        //    size_t size = bucket_size(bucket);
        //    printf("%3u] %4lu | ", bucket,size);
        //    for(auto it = begin(bucket); it != end(bucket); it++){
        //        auto &key = (*it).first;
        //        size_t hash = hasher.compute_hash(key.first,key.second);
        //        printf("[%lu] ",hash);
        //    }
        //    printf("\n");
        //}
        //printf("=============\n");
        //#endif

        printf("Entries         : %lu\n", size());
        printf("Buckets         : %lu\n", bucket_count());
        printf("Active buckets  : %lu\n", active_buckets);
        printf("Buckets > 1     : %lu\n", oversized_buckets);
        printf("Collisions      : %lu\n", collisions);
        printf("Collisions      : %lu\n", count_collisions());
        printf("Load factor     : %lf\n", load_factor());
        printf("Active load     : %lf\n", (double)size()/(double)active_buckets);
        printf("Max load factor : %lf\n", max_load_factor());

        printf("\n");
    }

    #ifdef USE_UNORDERED_MAP
    inline computed_table_key_t get_key(weighted_node *a,weighted_node *b) noexcept {
        return std::move(std::make_pair(a,b));
    }
    inline bool contains_key(weighted_node *a,weighted_node *b){
        if(computed_table_base_t::find(get_key(a,b)) != end())
            return true;
        else return false;
    }
    inline bool contains(weighted_node *n){
        for(auto it = begin(); it != end(); ++it )
            if (it->second == n)
                return true;
        return false;
    }
    inline weighted_node*& find_or_reference(weighted_node *a, weighted_node *b){
        computed_table_key_t key = get_key(a,b);
        return operator[](key);
    }
    inline weighted_node* find(weighted_node *a, weighted_node *b){
        computed_table_key_t key = get_key(a,b);
        auto hit = computed_table_base_t::find(key);
        if(hit != computed_table_base_t::end())
            return hit->second;
        else
            return NULL;
    }
    inline void erase_value(weighted_node *n){
        for(auto it = begin(); it != end(); ++it){
            if (it->second == n){
                erase(it);
                break;
            }
        }
    }
    #else
    inline computed_table_key_t get_key(weighted_node *a,weighted_node *b) noexcept {
        return std::move(std::make_pair(a,b));
    }

    inline weighted_node** find_or_insert(weighted_node *a, weighted_node *b, weighted_node **n){
        computed_table_key_t key = get_key(a,b);

        weighted_node **& m = operator[](key);
        if(!m){
            m = n;
            return NULL;
        } else return m;
    }

    inline weighted_node** find(weighted_node *a, weighted_node *b){
        computed_table_key_t key = get_key(a,b);
        auto hit = computed_table_base_t::find(key);
        if(hit != computed_table_base_t::end())
            return hit->second;
        else
            return NULL;
    }

    #endif
};

#endif
