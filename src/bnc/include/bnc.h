#ifndef BNC_H
#define BNC_H

#include "basicbnc.h"
#include "sat.h"
#include <string>
#include "bound.h"

namespace bnc {
    template <bool COLLAPSE = true, bool DETERMINISM = false> node* conjoin(manager_t*, const unsigned int, node*, node*);
    template <bool COLLAPSE = true, bool DETERMINISM = false> node* conjoin(manager_t*, const unsigned int, node*, node*, ordering_t&);
    template <bool COLLAPSE = true> void add_clause_context(manager_t*, node*&, domain_closure_t&, support_t&, ordering_t&);
    template <bool COLLAPSE = true> node* constraints(manager_t*, domain_closure_t&, support_t&, ordering_t&);
    template <bool COLLAPSE = true> node* complete(manager_t*, node*, domain_closure_t&, support_t&, ordering_t&);
    template <bool COLLAPSE = true, bool DETERMINISM = false> node* solve(manager_t*, sat_t&, const ordering_t&);
    template <bool COLLAPSE = true> void bayesnode_to_wpbdds(bnc::manager_t*, std::queue<bnc::node*>&, BayesNode*, bnc::domain_closure_t&, support_t&, ordering_t&);
    void write_dot(manager_t*, node*, std::string aux = "", unsigned int id = 0, node *current = NULL);
    void write_dot(manager_t*, std::vector< node*>&, node *current = NULL);
    void write_pdf(std::string aux);
    void write_dot_pdf(manager_t *manager, node *n, std::string aux, node *current,bool wait = false);
    void write_dot_ps(manager_t *manager, node *n, std::string aux, node *current = NULL);
    void write_dot_ps(manager_t *manager, std::vector<node*>&, std::string aux, node *current = NULL);
    void write_ps(manager_t *manager, std::string aux = "");
    void write_bdd(manager_t*, node*, std::string aux = "");
    void write_bdd(manager_t*, std::vector<node*>&);
    bool equal(node*, node*);
    node* copy(manager_t*, node*);
    node* recursive_merge(manager_t*, node::table&, node*&);
    void print(node*);
    size_t size(node*);
    size_t operators(node*);
    void collapse(manager_t*, node*);

    // simulated annealing functions
    mpz_class get_upper_bound(bayesnet*, ordering_t&);
    mpz_class get_upper_bound(bound_t<false>&, ordering_t&);
    mpz_class get_upper_bound(bayesnet *bn, partitions_t& partitions, std::vector< ordering_t > &o);
    mpz_class get_upper_bound(bayesnet *bn, partition_t & partition, ordering_t &o);
}

#endif
