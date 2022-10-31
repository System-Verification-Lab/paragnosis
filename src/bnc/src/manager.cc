#include "manager.h"
#include "options.h"
#include "exceptions.h"
#include "defines.h"
#include <algorithm>
#include <string.h>
#include <limits>
#include "timer.h"
#include "bnc.h"
#include "partitioner.h"
#include <utility>
#include "architecturebound.h"
#include "composition.h"
#include <cassert>
#include "spanningtree.h"

using namespace bnc;

manager::manager(){
    actual_total_node_allocations = 0;
    actual_total_weight_allocations = 0;
    total_node_allocations = 0;
    total_weight_allocations = 0;
    bn = NULL;

    init();
}

manager::~manager(){
    while(!free_nodes.empty()){
        node* n = free_nodes.top();
        free_nodes.pop();
        free(n);
    }
    while(!free_weights.empty()){
        node::weights *w = free_weights.top();
        free_weights.pop();
        delete w;
    }
    // #ifdef VERBOSE
    // printf("Total node allocations    : %lu (%lu)\n", total_node_allocations, actual_total_node_allocations);
    // printf("Total weights allocations : %lu (%lu)\n", total_weight_allocations, actual_total_weight_allocations);
    // #endif
}

void manager::init(){
    join_compile_time    = -1;
    cpt_compile_time     = -1;
    total_compile_time   = -1;
    cumulative_nodes     = 0;
    cumulative_operators = 0;
    total_or_nodes       = 0;
    total_and_nodes      = 0;
    total_nodes          = 0;
    total_edges          = 0;
    total_operators      = 0;
}

bayesnet* manager::get_bayesnet(){
    return bn;
}

bayesgraph& manager::get_bayesgraph(){
    return g;
}

filename_t& manager::get_files() const {
    return files;
}

// domain_closure_t& manager::get_domain_closure(){
//     return domain_closure;
// }

size_t manager::get_nr_variables(){
    return bn->get_nr_variables();
}

void manager::generate_composition_ordering(){
    Composition comp;
    comp.SetBayesnet(bn);
    comp.SetPartitions(partitions);

    composition_ordering = comp.FindOrdering();

    printf("found:\n");
    printf("    ordering :");
    for(auto it = composition_ordering.begin(); it != composition_ordering.end(); it++)
        printf(" %d", *it);
    printf("\n");
    printf("    score    : %lf\n", comp.ComputeScore(composition_ordering));
}

void manager::generate_ordering(bayesnet *bn, bn_partitions_t& partitions){
    for(unsigned int i = 0; i < partitions.size(); i++){
        bn_partition_t &partition = partitions[i];

        if(partition.ordering.empty() && partition.variable_ordering.empty() && partition.pseudotree.Empty())
            generate_ordering(bn, partition);

        if(partition.variable_ordering.empty()){
            if(!partition.pseudotree.Empty())
                partition.variable_ordering = partition.pseudotree.GetPreOrdering();
            else if(!partition.ordering.empty())
                partition.variable_ordering = partition.ordering.get_literal_to_variable_ordering(g);
            else throw compiler_ordering_exception("Could not determine ordering");
        }

        assert(!partition.variable_ordering.empty());
        if(partition.ordering.empty()){
            const bool kIncludeWeightLiterals = !(OPT_BDD_TYPE == bdd_t::wpbdd || OPT_BDD_TYPE == bdd_t::tdmultigraph || OPT_BDD_TYPE == bdd_t::multigraph);
            partition.ordering = partition.variable_ordering.get_variable_to_literal_ordering(g, partition.partition, kIncludeWeightLiterals);
        }

        if(!partition.ordering.empty()){
            const std::vector<unsigned int>& literal_to_variable = g.get_literal_to_variable();
            partition.support_to_ordering.resize(g.get_nr_variables());
            for(unsigned int i = 0; i < partition.ordering.size(); i++){
                literal_t l = partition.ordering[i];
                if(!g.is_probability(l) && l != 0){
                    unsigned int variable = literal_to_variable[l];
                    partition.support_to_ordering[variable].push_back(std::make_pair(i,partition.ordering[i]));
                }
            }
        }

        if(partition.pseudotree.Empty()){
            if(OPT_BDD_TYPE == bdd_t::tdmultigraph)
                partition.pseudotree.SetTree(bn, partition.partition, partition.variable_ordering);
            else if(OPT_BDD_TYPE == bdd_t::multigraph)
                partition.pseudotree.SetChain(bn, partition.partition, partition.variable_ordering);
        }
    }
}

void manager::generate_ordering(bayesnet *bn, bn_partition_t& partition){
    ordering_t &ordering = partition.variable_ordering;

    if(OPT_ORDER < 0 && (OPT_READ_PSEUDO_ORDERING
        || OPT_READ_ORDERING
        || OPT_READ_VARIABLE_ORDERING
        || OPT_READ_ELIM_ORDERING)){

        std::string filename;
        if(OPT_READ_ORDERING){
            if(partitions.size() == 1)
                filename = files.get_filename(ORDERING);
            else filename = files.get_filename(ORDERING,stringf("%u",partition.partition.id));

            partition.ordering.read(filename);
            ordering = partition.ordering.get_literal_to_variable_ordering(g);
        } else if(OPT_READ_PSEUDO_ORDERING){
            if(partitions.size() == 1)
                filename = files.get_filename(PSEUDO_ORDERING);
            else filename = files.get_filename(PSEUDO_ORDERING,stringf("%u",partition.partition.id));

            partition.pseudotree.Read(filename);
            ordering = partition.pseudotree.GetMinFillOrdering();
        } else if(OPT_READ_ELIM_ORDERING){
            try {
            if(partitions.size() == 1)
                filename = files.get_filename(ELIM_ORDERING);
            else filename = files.get_filename(ELIM_ORDERING,stringf("%u",partition.partition.id));

            ordering.read(filename, bn);
            } catch (...) {}
        } else if(OPT_READ_VARIABLE_ORDERING){
            if(partitions.size() == 1)
                filename = files.get_filename(VARIABLE_ORDERING);
            else filename = files.get_filename(VARIABLE_ORDERING,stringf("%u",partition.partition.id));

            ordering.read(filename);
        }
    }

    if(ordering.empty()){
        if(OPT_ORDER < 0){
            if(OPT_BDD_TYPE == bdd_t::tdmultigraph){
                OPT_ORDER = 10;
            } else
                OPT_ORDER = 7;
        }
        //assert(!OPT_PARTITION && (OPT_ORDER == 7 || OPT_ORDER == 9 || OPT_ORDER == 10))

        printf("Determining ordering with strategy: ");

        Timer t;
        t.Start();
        switch(OPT_ORDER){
            case 0:
                printf("topological sort\n");
                ordering.generate_variable_ordering(bn,partition.partition);
                break;
            case 1:
                printf("breath first\n");
                ordering.generate_variable_ordering_1(bn,partition.partition);
                break;
            case 2:
                printf("reverse topological sort\n");
                ordering.generate_variable_ordering_2(bn,partition.partition);
                break;
            case 3:
                printf("branch and bound\n");
                ordering.generate_variable_ordering_3(bn,partition.partition);
                break;
            case 4:
                printf("brute force\n");
                ordering.generate_variable_ordering_4(bn,partition.partition);
                break;
            case 5:
                printf("topological sort on weakly connected components\n");
                ordering.generate_variable_ordering_5(bn,partition.partition);
                break;
            case 6:
                printf("lookahead branch and bound\n");
                ordering.generate_variable_ordering_6(bn,partition.partition);
                break;
            case 7:
                printf("simulated annealing (for chain ordering)\n");
                ordering.generate_variable_ordering_7(bn,partition.partition);
                break;
            case 8:
                printf("all topological sorts\n");
                ordering.generate_variable_ordering_8(bn,partition.partition);
                break;
            case 9:
                printf("simulated annealing (for tree ordering)\n");
                ordering.generate_variable_ordering_9(bn,partition.partition);
                break;
            case 10:
                printf("minhill\n");
                ordering = PseudoTree::MinFill(bn,partition.partition);
                break;
            default:
                throw compiler_ordering_exception("Order type not in valid range");
        }
        t.Stop();
        printf("Ordering determined in %.3fs\n", t.GetDuration<Timer::Seconds>());
    }

    //const std::vector<unsigned int>& literal_to_variable = g.get_literal_to_variable();
    //partition.support_to_ordering.resize(g.get_nr_variables());
    //for(unsigned int i = 0; i < ordering.size(); i++){
    //    literal_t l = ordering[i];
    //    if(!g.is_probability(l)){
    //        unsigned int variable = literal_to_variable[l];
    //        partition.support_to_ordering[variable].push_back(std::make_pair(i,ordering[i]));
    //    }
    //}
}

void manager::encode(bayesnet *bn){
    this->bn = bn;
    g.encode(bn);
}

void manager::load(const bdd_t OPT_BDD_TYPE){
    assert(bn != NULL);
    // get (initial) variable ordering
    if(OPT_PARTITION){
        if(OPT_READ_PARTITION){
            // read from .part
            partitions.read(files.get_filename_c(PARTITION),bn);
        } else {
            // detect proper partitions
            Timer t;
            ordering_t ordering;
            // generate initial ordering for entire network

            printf("Determining ordering with strategy: ");
            printf("simulated annealing\n");
            t.Start();
            partitioner::generate_partitions_sa(bn, partitions);
            t.Stop();
            printf("Partitions determined in %.3fs\n", t.GetDuration<Timer::Seconds>());

            //partitions.SetOnePartition(bn);
            //t.Start();
            //ordering.generate_variable_ordering_7(bn,partitions[0].partition);
            //t.Stop();
            //printf("Ordering determined in %.3fs\n", t.GetDuration<Timer::Seconds>());

            //printf("Determining partitions\n");
            //t.Start();
            //partitioner::generate_partitions(bn, &ordering, partitions);
            //t.Stop();
            //printf("Partitions determined in %.3fs\n", t.GetDuration<Timer::Seconds>());
        }

        if(OPT_WRITE_COMPOSITION_ORDERING){
            printf("Determining composition ordering\n");
            Timer t;
            t.Start();
            generate_composition_ordering();
            t.Stop();
            printf("Composition ordering determined in %.3fs\n", t.GetDuration<Timer::Seconds>());
        }
    } else partitions.SetOnePartition(bn);

    // set partition ids
    for(unsigned int i = 0; i < partitions.size(); i++)
        partitions[i].partition.id = i;

    partitions.verify(bn);

    if(!OPT_NO_ORDERING){
        generate_ordering(bn, partitions);
        partitions.VerifyOrdering(bn);
    }

    if(OPT_SHOW_SCORE){
        printf("================ SCORE ================\n");
        for(unsigned int i = 0; i < partitions.size(); i++){
            if(partitions.size() > 1)
                printf(" = Partition %lu\n", i);

            bn_partition_t &partition = partitions[i];
            const partition_t &kPartition = partitions[i].partition;
            partition.variable_ordering.stats(bn,kPartition);

            printf("Tree bound: %lu\n", Bound::ComputeTree(bn,kPartition,partition.variable_ordering));

            if(partition.pseudotree.Empty())
                partition.pseudotree.SetTree(bn, partition.partition, partition.variable_ordering);

            SpanningTree spanningtree;
            spanningtree.SetBayesnet(get_bayesnet());
            spanningtree.SetBayesgraph(&get_bayesgraph());
            spanningtree.Init(partition.partition);
            spanningtree.Create(partition.pseudotree);
            printf("\n");
            spanningtree.PrintStats();
        }
        printf("=======================================\n");
    }
}

void manager::set_ordering(ordering_t& o){
    if(partitions.size() != 1)
        partitions.resize(1);

    bn_partition_t &partition = *partitions.begin();
    ordering_t &ordering = partition.ordering;
    ordering = o;

    if(ordering.empty())
        throw compiler_ordering_exception("Could not determine ordering");

    const std::vector<unsigned int>& literal_to_variable = g.get_literal_to_variable();

    partition.support_to_ordering.clear();
    partition.support_to_ordering.resize(ordering.size());
    for(unsigned int i = 0; i < ordering.size(); i++){
        literal_t l = ordering[i];
        if(!g.is_probability(l)){
            unsigned int variable = literal_to_variable[l];
            partition.support_to_ordering[variable].push_back(std::make_pair(i,ordering[i]));
        }
    }
}

inline manager::node::weights* manager::allocate_weights(){
    #ifdef VERBOSE
    actual_total_weight_allocations++;
    #endif

    node::weights *w = new node::weights;
    #ifdef memory_safe
    if(!w)
        throw compiler_debug_exception("out-of-memory");
    #endif

    return w;
}

inline manager::node* manager::allocate(){
    #ifdef VERBOSE
    actual_total_node_allocations++;
    #endif

    node *n = (node*) malloc(sizeof(node));
    #ifdef MEMORY_SAFE
    if(!n)
        throw compiler_debug_exception("out-of-memory");
    #endif

    return n;
}

manager::node* manager::create(const bool initialize){
    #ifdef VERBOSE
    total_node_allocations++;
    #endif

    node *n;
    #ifndef SIMPLE_ALLOCATION
    if(!free_nodes.empty()){
        n = free_nodes.top();
        free_nodes.pop();
    } else
    #endif
        n = allocate();

    if(initialize)
        node::init(n);

    return n;
}

manager::node::weights* manager::create_weights(){
    #ifdef VERBOSE
    total_weight_allocations++;
    #endif

    node::weights *w;
    #ifndef SIMPLE_ALLOCATION
    if(!free_weights.empty()){
        w = free_weights.top();
        free_weights.pop();
        (*w).clear();
    } else
    #endif
        w = allocate_weights();

    return w;
}

void manager::reserve_nodes(unsigned int RESERVE){
    node *nodes = (node*) malloc(RESERVE * sizeof(node));
    if(!nodes)
        throw compiler_debug_exception("could not allocate reserve");

    for(unsigned int i = 0; i < RESERVE; i++)
        free_nodes.push(&(nodes[RESERVE-1-i]));

    //for(unsigned int i = 0; i < RESERVE; i++)
    //    free_nodes.push(nodes[i]);
    //while(free_nodes.size() < RESERVE){
    //    #ifdef VERBOSE
    //    total_node_allocations++;
    //    #endif
    //    free_nodes.push(allocate());
    //}
}

void manager::reserve_weights(unsigned int RESERVE){
    while(free_weights.size() < RESERVE){
        #ifdef VERBOSE
        total_weight_allocations++;
        #endif
        free_weights.push(allocate_weights());
    }
}


const ordering_t& manager::get_ordering(){
    return get_ordering(0);
}

const ordering_t& manager::get_ordering(unsigned int i){
    return partitions[i].ordering;
}

bn_partition_t & manager::get_partition(unsigned int i){
    return partitions[i];
}

const ordering_t manager::create_ordering(const unsigned int partition_id, const support_t& support) {
    bn_partition_t &partition = get_partition(partition_id);

    support_set_t order_set;
    for(auto it = support.begin(); it != support.end(); it++){
        unsigned int variable = *it;
        order_set.insert(partition.support_to_ordering[variable].begin(), partition.support_to_ordering[variable].end());
    }

    ordering_t ordering;
    for(auto it = order_set.begin(); it != order_set.end(); it++)
        ordering.push_back((*it).second);

    return ordering;
}

const ordering_t manager::create_ordering(const unsigned int partition_id, node *m, node *n) {
    const support_t support = create_support(m,n);
    return std::move(create_ordering(partition_id, support));
}

const ordering_t manager::create_ordering(const unsigned int partition_id, node *n) {
    auto hit = node_to_support.find(n);
    if(hit != node_to_support.end())
        return std::move(create_ordering(partition_id, hit->second));
    else return ordering_t();
}

const support_t manager::create_variable_support(unsigned int variable) {
    const std::vector<unsigned int> &literal_to_variable = g.get_literal_to_variable();
    BayesNode *gn = g.get_node(variable);
    unsigned int variables = gn->GetNrVariables();
    support_t support;
    for(unsigned int i = 0; i < variables; i++){
        literal_t l = gn->GetLiteral(i);
        support.insert(literal_to_variable[l]);
    }
    return std::move(support);
}

void manager::add_support(node *n, unsigned int variable) {
    node_to_support[n] = create_variable_support(variable);
}

void manager::add_support(node *n, support_t &support) {
    node_to_support[n] = support;
}

const support_t manager::get_support(node *n) {
    auto hit = node_to_support.find(n);
    if(hit != node_to_support.end())
        return hit->second;
    else return support_t();
}

const support_t manager::create_support(node *m, node *n) {
    auto hitm = node_to_support.find(m);
    auto hitn = node_to_support.find(n);

    if(hitm != node_to_support.end() && hitn != node_to_support.end()){
        support_t &supportm = hitm->second;
        support_t &supportn = hitn->second;

        support_t support;
        support.insert(supportm.begin(), supportm.end());
        support.insert(supportn.begin(), supportn.end());

        return std::move(support);
    } else throw compiler_manager_exception("could not find support set");
}

manager::node* manager::create_terminal(const bool satisfiable){
    node *n = create(true);
    if(satisfiable)
        n->l = -1;
    return n;
}

void manager::destroy_weights(node::weights *&w){
    if(w){
        #ifndef SIMPLE_ALLOCATION
        free_weights.push(w);
        #else
        delete w;
        #endif
        #ifdef DEBUG
        w = NULL;
        #endif
    }
}

void manager::force_destroy(node *&n){
    if(!node::is_terminal(n)){
        node::dereference(n->t);
        node::dereference(n->e);
    }
    destroy_weights(n->W);

    #ifndef SIMPLE_ALLOCATION
    #ifdef DEBUG
    node::init(n);
    #endif
    n->ref = 0;
    free_nodes.push(n);
    #else
    free(n);
    #endif
    #ifdef DEBUG
    n = node::negate(NULL);
    #endif
}

void manager::destroy_node(node *&n){
    destroy_weights(n->W);

    #ifndef SIMPLE_ALLOCATION
    #ifdef DEBUG
    node::init(n);
    #endif
    free_nodes.push(n);
    #else
    free(n);
    #endif
    #ifdef DEBUG
    n = node::negate(NULL);
    #endif
}

void manager::destroy(node *&n){
    if(n->is_dead(n)){
        if(!node::is_terminal(n)){
            node::dereference(n->t);
            node::dereference(n->e);
        }
        destroy_weights(n->W);

        #ifndef SIMPLE_ALLOCATION
        #ifdef DEBUG
        node::init(n);
        #endif
        free_nodes.push(n);
        #else
        free(n);
        #endif
        #ifdef DEBUG
        n = node::negate(NULL);
        #endif
    }
}

void manager::destroy_support(node *n){
    node_to_support.erase(n);
}

void manager::recursive_destroy(node *&n){
    node::stack s;
    if(n && node::is_dead(n))
        s.push(n);

    while(!s.empty()){
        node *n = s.top();
        s.pop();

        if(!node::is_terminal(n)){
            node *&t = n->t;
            node *&e = n->e;

            node::dereference(e);
            if(node::is_dead(e))
               s.push(e);

            node::dereference(t);
            if(node::is_dead(t))
               s.push(t);
        }
        destroy_node(n);
    }
    #ifdef DEBUG
    n = node::negate(NULL);
    #endif
}

void manager::deallocate(node *n){
    node::set done;
    node::stack s;
    s.push(n);
    done.insert(NULL);
    while(!s.empty()){
        node *n = s.top();

        if(!done.contains(n)){
            done.insert(n);
            if(n->e)
                s.push(n->e);
            if(n->t)
                s.push(n->t);

            destroy_weights(n->W);
            free(n);

        } else s.pop();
    }
}

std::string manager::get_literal_name(const literal_t &l) {
    const unsigned int LITERALS = g.get_nr_literals();

    if(l <= (literal_t) LITERALS){
        // This is a variable-value literal
        variable_t variable = g.get_literal_to_variable()[l];
        literal_t base_l = g.get_variable_to_literal()[variable];
        unsigned int value = l - base_l;

        std::string variable_name = get_bayesnet()->get_node_name(variable);
        std::string value_name = get_bayesnet()->get_node_value_name(variable,value);

        return stringf("%u: %s=%s",l , variable_name.c_str(), value_name.c_str());
    } else {
        // this is a probability
        const std::vector<probability_t> &probability =  g.get_weight_to_probability();
        probability_t p = probability[l-(LITERALS+1)];
        return stringf("%u: %.3f", p, l);
    }
}

std::vector<std::string> manager::get_literal_names(const biordering_t &biordering) {
    std::vector<std::string> names(biordering.back.size());

    for(unsigned int i = 0; i < biordering.back.size(); i++){
        const literal_t &l = biordering.back[i];
        names[i] = get_literal_name(l);
    }

    return std::move(names);
}


