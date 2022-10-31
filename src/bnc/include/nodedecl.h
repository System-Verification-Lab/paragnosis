// ==== Hash functions =========================================================

typedef unsigned int hash_t;

inline static bool identical(const weighted_node *m, const weighted_node *n){
    if(m->l != n->l
        || m->t != n->t
        || m->e != n->e)
        return false;

    if(m->W != NULL){
        if(n->W != NULL){
            return *(m->W) == *(n->W);
        } else return false;
    } else if (n->W != NULL)
        return false;

    return true;
}


inline static bool equal(const weighted_node *m, const weighted_node *n){
    if(m->l != n->l)
        return false;

    if(m->W != NULL){
        if(n->W != NULL){
            return *(m->W) == *(n->W);
        } else return false;
    } else if (n->W != NULL)
        return false;

    return true;
}

struct ptr_hash {
    hash_t operator()(const weighted_node *n) const {
        return std::hash<hash_t>()((uintptr_t) n);
    }
};

struct ptr_equal {
    bool operator()(const weighted_node *m, weighted_node *n) const {
        return m == n;
    }
};

struct obj_hash {
    hash_t operator()(const weighted_node *n) const {
        return std::hash<hash_t>()(
                n->l
                + (uintptr_t) n->t
                + (uintptr_t) n->e
                // + n->w.get_sum()
                );
    }
};

struct obj_equal {
    bool operator()(const weighted_node *m, weighted_node *n) const {
        return identical(m,n);
    }
};

// ==== Data structures ========================================================

#include "computedtable.h"

struct weighted_node_map : std::unordered_map<weighted_node*,weighted_node*,ptr_hash, ptr_equal> {
    inline weighted_node*& find_or_reference(weighted_node *n){
        return operator[](n);
    }
};


struct weighted_node_set : public std::unordered_set<weighted_node*,ptr_hash, ptr_equal> {
    inline bool contains(weighted_node *n) const {
        return find(n) != end();
    }
};

struct weighted_node_table : public std::unordered_set<weighted_node*,obj_hash, obj_equal> {
    inline weighted_node* find_or_insert(weighted_node *n){
        auto it = this->find(n);
        if(it == this->end()){
            this->insert(n);
            return NULL;
        } else {
            if(n != *it)
                return *it;
            else return NULL;
        }
    }
};

// ==== Typedefs ===============================================================

typedef std::stack<weighted_node**> weighted_node_reference_stack;
typedef std::stack<weighted_node*> weighted_node_stack;
typedef std::queue<weighted_node*> weighted_node_queue;
//typedef std::unordered_map<weighted_node*, support_set,ptr_hash, ptr_equal> weighted_node_support_map;
typedef std::map<weighted_node*, support_set> weighted_node_support_map;

