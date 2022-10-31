#include "ite.h"
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "options.h"
#include "exceptions.h"
#include "defines.h"

#define negated signbit

using namespace bnc;

ite::ite(bnc::node::computed_table& ctable): table(ctable){
    init();
}

ite::ite(manager_t *m): table(m->ctable){
    init();
    manager = m;
}

void ite::init(){
    #ifdef DEBUG
    cache_hits = 0;
    hashmap_resizes = 0;
    #endif
    W = NULL;
    cached_node = NULL;
    BDDS = 0;

    stack.clear();
    stack.push_back({NULL, NULL});
}

ite::~ite(){

}


void ite::set_manager(manager_t *manager){
    this->manager = manager;
}

void ite::clear_weight(){
    bnc::destroy_weights(manager, W);
    W = NULL;
}

void ite::stats(){
    #ifdef DEBUG
    printf("cache hits      : %u/%u\n",cache_hits, table.size());
    printf("Hashmap resizes : %u\n",hashmap_resizes);
    #endif
    table.stats();
}

void ite::add_bdd(node* n){
    if(BDDS >= MAX_BDDS)
        throw compiler_debug_exception("ITE can only take %u wpbdds", MAX_BDDS);

    node_block_t &bdd = stack.back();
    bdd[BDDS] = n;
    BDDS++;
}

node::weights* ite::get_weights(){
    node::weights *tmp = W;
    W = NULL;
    return tmp;
}

const bnc::node::computed_table& ite::get_computed_table(){
    return table;
}

void ite::remove_cache(node *n){
    table.erase_value(n);
}

void ite::add_cache(node *n){
    #ifndef WITHOUT_CACHING
    const node_block_t & top = stack.back();

    #ifdef DEBUG
    size_t current_hashmap_size = table.bucket_count();
    #endif

    node *&hit = table.find_or_reference(top[0], top[1]);

    #ifdef DEBUG
    if(current_hashmap_size != table.bucket_count())
        hashmap_resizes++;
    #endif

    if(!hit)
        hit = n;
    #ifdef DEBUG
    else {
        if(n != hit)
            fprintf(stderr, "item overwritten in cache\n");
        else
            fprintf(stderr, "item already in cache\n");
    }
    #endif
    #endif
}

bool ite::is_cached(){
    #ifndef WITHOUT_CACHING
    const node_block_t & top = stack.back();
    cached_node = table.find(top[0], top[1]);
    #ifdef DEBUG
    if(cached_node)
        cache_hits++;
    #endif
    return cached_node != NULL;
    #else
    return false;
    #endif
}

node* ite::get_current(unsigned int i, unsigned int context){
    if(stack.size() <= context)
        return NULL;

    if(context > 0 && stack.size() > context)
        return stack.rbegin()[context][i];
    else return stack.rbegin()[0][i];
}

node* ite::get_cached(){
    #ifndef WITHOUT_CACHING
    return cached_node;
    #else
    return NULL;
    #endif
}

inline void ite::add_weight(const bnc::node * const &n){
    if(n->W){
        if(!W) W = bnc::create_weights(manager);
        W->insert(n->W->begin(), n->W->end());
    }
}

template <int i>
node * seek_parent(const ite_stack_t &stack){
    const node_block_t & top = stack.back();
    node* child = top[i];
    for(auto it = stack.rbegin()+1; it != stack.rend(); it++){
        node* parent = (*it)[i];
        if(parent != child)
            return parent;
    }
    return child;
}

template <>
satisfy_t ite::condition<true,false,false>(const literal_t l){
    // Positive = true, collapse = false, determinism = false
    const node_block_t & top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();
    if(bdd[0]->l == l){
        add_weight(bdd[0]);
        bdd[0] = bdd[0]->t;
    }
    if(bdd[1]->l == l){
        add_weight(bdd[1]);
        bdd[1] = bdd[1]->t;
    }

    #ifdef DEBUG
    if(node::is_unsatisfiable(bdd[0]) || node::is_unsatisfiable(bdd[1]))
        throw compiler_debug_exception("Negative cofactor at true edge");
    #endif

    if(is_cached())
        return cached;

    // Assumption: we never encounter a false terminal at a 'then' edge
    if(node::is_terminal(bdd[0]) || node::is_terminal(bdd[1])){
        if(node::is_satisfiable(bdd[0]) && node::is_satisfiable(bdd[1]))
            return satisfiable;
        else return trivial;
    } else return unsatisfied;
}

template <>
satisfy_t ite::condition<false,false,false>(const literal_t l){
    // Positive = false, collapse = false, determinism = false
    const node_block_t &top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();
    if(bdd[0]->l == l)
        bdd[0] = bdd[0]->e;
    if(bdd[1]->l == l)
        bdd[1] = bdd[1]->e;

    // Assumption: we never encounter a true terminal at an 'else' edge
    if(node::is_unsatisfiable(bdd[0]) || node::is_unsatisfiable(bdd[1]))
        return unsatisfiable;
    else return unsatisfied;
}

template <>
satisfy_t ite::condition<true,false,true>(const literal_t l){
    // Positive = true, collapse = false, determinism = true
    const node_block_t &top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();
    if(bdd[0]->l == l){
        add_weight(bdd[0]);
        bdd[0] = bdd[0]->t;
    }
    if(bdd[1]->l == l){
        add_weight(bdd[1]);
        bdd[1] = bdd[1]->t;
    }

    if(is_cached())
        return cached;

    if(node::is_terminal(bdd[0]) || node::is_terminal(bdd[1])){
        if(node::is_unsatisfiable(bdd[0]) || node::is_unsatisfiable(bdd[1]))
            return unsatisfiable;
        else if(node::is_satisfiable(bdd[0]) && node::is_satisfiable(bdd[1]))
            return satisfiable;
        else // one of both is a true terminal
            return trivial;
    } else return unsatisfied;
}

template <>
satisfy_t ite::condition<false,false,true>(const literal_t l){
    // Positive = false, collapse = false, determinism = true
    return condition<false,false,false>(l);
}


template <>
satisfy_t ite::condition<true,true,true>(const literal_t l){
    // Positive = true, collapse = true, determinism = true
    const node_block_t &top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();

    if(bdd[0]->l == l){
        add_weight(bdd[0]);
        bdd[0] = bdd[0]->t;
    } else {
        node *parent = seek_parent<0>(stack);
        if(node::is_negative_cofactor(parent,bdd[0])){
            add_weight(parent);
            bdd[0] = parent->t;
        }
    }

    if(bdd[1]->l == l){
        add_weight(bdd[1]);
        bdd[1] = bdd[1]->t;
    } else {
        node *parent = seek_parent<1>(stack);
        if(node::is_negative_cofactor(parent,bdd[1])){
            add_weight(parent);
            bdd[1] = parent->t;
        }
    }

    if(is_cached())
        return cached;

    if(node::is_terminal(bdd[0]) || node::is_terminal(bdd[1])){
        if(node::is_unsatisfiable(bdd[0]) || node::is_unsatisfiable(bdd[1]))
            return unsatisfiable;
        else if(node::is_satisfiable(bdd[0]) && node::is_satisfiable(bdd[1]))
            return satisfiable;
        else // one of both is a true terminal
            return trivial;
    } else return unsatisfied;
}


template <>
satisfy_t ite::condition<true,true,false>(const literal_t l){
    // Positive = true, collapse = true, determinism = false
    const node_block_t &top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();

    if(bdd[0]->l == l){
        add_weight(bdd[0]);
        bdd[0] = bdd[0]->t;
    } else {
        const node * const & parent = seek_parent<0>(stack);
        if(node::is_negative_cofactor(parent,bdd[0])){
            add_weight(parent);
            bdd[0] = parent->t;
        }
    }

    if(bdd[1]->l == l){
        add_weight(bdd[1]);
        bdd[1] = bdd[1]->t;
    } else {
        const node * const & parent = seek_parent<1>(stack);
        if(node::is_negative_cofactor(parent,bdd[1])){
            add_weight(parent);
            bdd[1] = parent->t;
        }
    }

    #ifdef DEBUG
    if(node::is_unsatisfiable(bdd[0]) || node::is_unsatisfiable(bdd[1]))
        throw compiler_debug_exception("Negative cofactor at true edge");
    #endif

    if(is_cached())
        return cached;

    // Assumption: we never encounter a false terminal at a 'then' edge
    if(node::is_satisfiable(bdd[0]) || node::is_satisfiable(bdd[1])){
        if(node::is_satisfiable(bdd[0]) && node::is_satisfiable(bdd[1]))
            return satisfiable;
        else return trivial;
    } else return unsatisfied;
}



template <>
satisfy_t ite::condition<false,true,true>(const literal_t l){
    // Positive = false, collapse = true, determinism = true
    const node_block_t &top = stack.back();
    if(top[0]->l != l && top[1]->l != l)
        return redundant;

    stack.push_back(top);
    node_block_t &bdd = stack.back();

    if(bdd[0]->l == l){
        bdd[0] = bdd[0]->e;
    }
    if(bdd[1]->l == l){
        bdd[1] = bdd[1]->e;
    }

    // Assumption: we never encounter a true terminal at an 'else' edge
    if((node::is_unsatisfiable(bdd[0]) && !node::is_negative_cofactor(seek_parent<1>(stack),bdd[1]))
            || (node::is_unsatisfiable(bdd[1]) && !node::is_negative_cofactor(seek_parent<0>(stack),bdd[0]))
            || (node::is_unsatisfiable(bdd[0]) && node::is_unsatisfiable(bdd[1])))
        return unsatisfiable;
    else return unsatisfied;
}

template <>
satisfy_t ite::condition<false,true,false>(const literal_t l){
    // Positive = false, collapse = true, determinism = false
    return condition<false,true,true>(l);
}

node* ite::get_trivial(){
    node_block_t &bdd = stack.back();
    if(node::is_satisfiable(bdd[0]))
        return bdd[1];
    else return bdd[0];
}

void ite::undo(){
    stack.pop_back();
}

template satisfy_t ite::condition<true,true,true>(const literal_t);
template satisfy_t ite::condition<true,true,false>(const literal_t);
template satisfy_t ite::condition<true,false,true>(const literal_t);
template satisfy_t ite::condition<true,false,false>(const literal_t);
template satisfy_t ite::condition<false,true,true>(const literal_t);
template satisfy_t ite::condition<false,true,false>(const literal_t);
template satisfy_t ite::condition<false,false,true>(const literal_t);
template satisfy_t ite::condition<false,false,false>(const literal_t);

