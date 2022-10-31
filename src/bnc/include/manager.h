#ifndef MANAGER_H
#define MANAGER_H

#include "node.h"
#include "ordering.h"
#include "bayesgraph.h"
#include "support.h"
#include <unordered_map>
#include <bn-to-cnf/bayesnet.h>
#include "types.h"
#include "partition.h"

class manager {
    public:
        typedef weighted_node node;

        manager();
        ~manager();

        void encode(bayesnet *bn);
        void load(bn_partitions_t &, const bdd_t);
        void load(const bdd_t);
        bayesnet* get_bayesnet();
        bayesgraph& get_bayesgraph();
        filename_t& get_files() const;
        std::string get_literal_name(const literal_t &);
        std::vector<std::string> get_literal_names(const biordering_t &biordering);

        void reserve_nodes(unsigned int);
        void reserve_weights(unsigned int);

        node* create(const bool initialize  =  true);
        node* create_terminal(const bool);
        node::weights* create_weights();
        void destroy_node(node*&);
        void force_destroy(node*&);
        void destroy(node*&);
        void destroy_support(node*);
        void destroy_weights(node::weights*&);
        void recursive_destroy(node*&);
        void deallocate(node*);

        size_t get_nr_variables();
        const support_t create_variable_support(unsigned int);
        void add_support(node*, unsigned int);
        void add_support(node*, support_t &);
        const support_t create_support(node*,node*);
        const support_t get_support(node*);
        const ordering_t& get_ordering(unsigned int);
        const ordering_t& get_ordering();
        void set_ordering(ordering_t&);
        bn_partition_t& get_partition(unsigned int);
        const ordering_t create_ordering(const unsigned int, node*);
        const ordering_t create_ordering(const unsigned int, node*, node*);
        const ordering_t create_ordering(const unsigned int, const support_t&);
        void create_ordering_support(const unsigned int, const support_t&);
        void generate_ordering(bayesnet*, bn_partition_t&);
        void generate_ordering(bayesnet*, bn_partitions_t&);
        void generate_composition_ordering();
        //domain_closure_t& get_domain_closure();

        void init();
        double compile_time_ms;
        double spanning_tree_time_ms;
        double join_compile_time;
        double cpt_compile_time;
        double total_compile_time;
        double total_compile_time_ms;
        size_t cumulative_nodes;
        size_t cumulative_operators;
        size_t total_or_nodes;
        size_t total_and_nodes;
        size_t total_nodes;
        size_t total_edges;
        size_t total_operators;
        bn_partitions_t partitions;
        ordering_t composition_ordering;
        node::computed_table ctable;

        std::vector<biordering_t> biorderings;
    private:
        //domain_closure_t domain_closure;
        node* allocate();
        node::weights* allocate_weights();
        bayesgraph g;
        bayesnet *bn;
        //struct vstack : public std::vector<node*> {
        //    inline void push(node*n){ push_back(n); };
        //    inline void pop() { pop_back(); };
        //    inline node*& top() { return back(); };
        //};
        //vstack free_nodes;
        node::stack free_nodes;
        node::weights::stack free_weights;
        std::unordered_map<node*, support_t, node::ptr_hash, node::ptr_equal> node_to_support;

        size_t actual_total_node_allocations;
        size_t actual_total_weight_allocations;
        size_t total_node_allocations;
        size_t total_weight_allocations;

};

#endif
