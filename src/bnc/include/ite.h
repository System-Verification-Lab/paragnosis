#ifndef ITE_H
#define ITE_H

#include "satisfy.h"
#include <bn-to-cnf/cnf.h>
#include <stack>
#include <vector>
#include <array>
#include "basicbnc.h"

#define MAX_BDDS  2
struct node_block_t {
    bnc::node* nodes[MAX_BDDS];
    bnc::node* const & operator[](unsigned int i) const { return nodes[i]; };
    bnc::node*& operator[](unsigned int i) { return nodes[i]; };
};

typedef std::vector< node_block_t > ite_stack_t;

class ite {
    public:
        ite(bnc::node::computed_table&);
        ite(bnc::manager_t*);
        ~ite();
        void init();
        void add_bdd(bnc::node*);
        void set_manager(bnc::manager_t*);
        template <bool POSITIVE, bool COLLAPSE = true, bool DETERMINISM = false> satisfy_t condition_(const literal_t);
        template <bool POSITIVE, bool COLLAPSE = true, bool DETERMINISM = false> satisfy_t condition(const literal_t);
        void undo();
        bnc::node::weights* get_weights();
        void add_cache(bnc::node*);
        void remove_cache(bnc::node*);
        bool is_cached();
        bnc::node* get_cached();
        void stats();
        void clear_weight();
        const bnc::node::computed_table& get_computed_table();
        bnc::node* get_trivial();
        bnc::node* get_current(unsigned int,unsigned int context = 0);
    private:
        void add_weight(const bnc::node * const &);
        bnc::node::computed_table &table;
        ite_stack_t stack;

        unsigned int BDDS;
        bnc::node::weights *W;
        bnc::node *cached_node;
        #ifdef DEBUG
        unsigned int cache_hits;
        size_t hashmap_resizes;
        #endif
        bnc::manager_t *manager;
};

typedef ite ite_t;

#endif
