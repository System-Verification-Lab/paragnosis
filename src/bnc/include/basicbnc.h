#ifndef BASICBNC_H
#define BASICBNC_H

#include "support.h"
#include "manager.h"
#include "node.h"
#include "domainclosure.h"

namespace bnc {
    typedef domain_closure domain_closure_t;
    typedef manager manager_t;
    typedef weighted_node node;

    inline node* create(manager_t *manager, const bool initialize = true){
        return manager->create(initialize);
    }

    inline node* create_terminal(manager_t *manager, const bool satisfiable){
        return manager->create_terminal(satisfiable);
    }

    inline node::weights* create_weights(manager_t *manager){
        return manager->create_weights();
    }

    inline void force_destroy(manager_t *manager, node *&n){
        manager->force_destroy(n);
    }

    inline void destroy(manager_t *manager, node *&n){
        manager->destroy(n);
    }

    inline void destroy_weights(manager_t *manager, node::weights *&w){
        manager->destroy_weights(w);
    }

    inline void recursive_destroy(manager_t *manager, node *&n){
        manager->recursive_destroy(n);
    }

    inline void destroy_support(manager_t *manager, node *n){
        manager->destroy_support(n);
    }

    inline void deallocate(manager_t *manager, node *n){
        manager->deallocate(n);
    }

    inline void add_support(manager_t *manager, node *n, support_t &support){
        manager->add_support(n, support);
    }

    inline support_t create_variable_support(manager_t *manager, unsigned int variable){
        return std::move(manager->create_variable_support(variable));
    }

    inline void add_support(manager_t *manager, node *n, unsigned int variable){
        manager->add_support(n, variable);
    }

    inline const support_t get_support(manager_t *manager, node *n){
        return std::move(manager->get_support(n));
    }

    inline const support_t create_support(manager_t *manager, node *n, node *m){
        return std::move(manager->create_support(n,m));
    }

    inline const ordering_t create_ordering(manager_t *manager, const unsigned int partition_id, const support_t &support){
        return std::move(manager->create_ordering(partition_id, support));
    }

    inline void computed_table_rehash(manager_t *manager,const size_t SIZE){
        // set number of buckets
        manager->ctable.rehash(SIZE);
    }

    inline void computed_table_reserve(manager_t *manager,const size_t SIZE){
        // set number of buckets
        manager->ctable.reserve(SIZE);
    }

    inline void computed_table_clear(manager_t *manager){
        manager->ctable.clear();
    }

    inline void computed_table_set_load_factor(manager_t *manager,float factor){
        manager->ctable.max_load_factor(factor);
    }

    inline float computed_table_get_load_factor(manager_t *manager, float factor){
        return manager->ctable.max_load_factor();
    }
    // inline domain_closure_t& domain_closure(manager_t *manager){
    //     return manager->get_domain_closure();
    // }
}

#endif
