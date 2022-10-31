#include "compiler.h"
#include <vector>
#include "compiler.h"
#include "timer.h"
#include "options.h"
#include "defines.h"
#include "threading.h"
#include <thread>
#include "mutex.h"
#include "misc.h"
#include <algorithm>
#include "bnc.h"
#include <exception/stringf.h>

using namespace std;
using namespace bnc;

struct compiler_data_t {
    bnc::manager_t *manager;
    vector<sat_t> sat;
    vector<domain_closure_t> closure;
    vector<support_t> support;
    vector<ordering_t> ordering;
    vector<ite_t> ite;
    int id;
    robust_cond_mutex mutex;
    compiler *c;
    std::vector<unsigned int> counter;
    synchronous_queue<unsigned int> q_init;
    synchronous_queue<unsigned int> q_td;
    synchronous_queue<bnc::node*> q_bu;
    std::vector<bnc::node*> bdds;
    Timer cpt_timer;
    Timer join_timer;
    Timer total_timer;
    size_t CUMULATIVE_CPT_SIZE;
    size_t CUMULATIVE_OPT_SIZE;
};

template <bool COLLAPSE, bool DETERMINISM, compilation_t COMPILATION_TYPE>
void* compiler::compile_async(void *v){
    compiler_data_t &data = *((compiler_data_t*)v);
    compiler &c = *(data.c);
    bnc::manager_t &manager = c.manager;
    const unsigned int partition_id = 0;
    vector<node*> cpts(manager.get_nr_variables());
    bn_partitions_t &partitions = manager.partitions;

    data.total_timer.Start();
    if(COMPILATION_TYPE == compilation_t::topdown_bottomup){
        // ==================  init sat solvers ============================
        bayesgraph &g = manager.get_bayesgraph();
        for(unsigned int variable = 0; variable < data.counter[0]; variable++){

            const bool kIncludeWeightLiterals = c.BDD_TYPE != bdd_t::wpbdd || OPT_PARALLELISM;
            data.ordering[variable].set_ordering(manager.get_ordering(partitions.variable_to_partition_id(variable)),g.get_node(variable),kIncludeWeightLiterals);
            data.sat[variable].set_manager(&manager);
            data.sat[variable].init<DETERMINISM>(g, variable);
        }


        // =========================== compile topdown per cpt ===================
        data.cpt_timer.Start();
        //synchronous_queue<bnc::node*> q_compare;
        for(unsigned int variable = 0; variable < data.counter[1]; variable++){
            #ifdef VERBOSE
            Timer cpt_local_timer;
            cpt_local_timer.Start();
            #endif

            bnc::node *cpt = bnc::solve<COLLAPSE,DETERMINISM>(&manager, data.sat[variable], data.ordering[variable]);
            #ifdef DEBUG
            if(!cpt)
                throw compiler_debug_exception("topdown compiler did not create WPBDD for CPT %u\n", variable);
            #endif
            bnc::add_support(&manager,cpt,variable);

            cpts[variable] = cpt;
            data.q_bu.push_async(cpt);

            #ifdef VERBOSE
            cpt_local_timer.Stop();
            printf("Compiled CPT %u/%u (%.3lfs)\n", variable+1, data.counter[0],cpt_local_timer.GetDuration<Timer::Seconds>());
            #endif

        }
        data.cpt_timer.Stop();

    } else {
        // init closure and orderings
        for(unsigned int variable = 0; variable < data.counter[0]; variable++){
            data.support[variable] = create_variable_support(&manager, variable);
            data.ordering[variable] = create_ordering(&manager, partition_id, data.support[variable]);
        }
        domain_closure_t closure;
        bayesgraph &g = manager.get_bayesgraph();
        closure.init(&g);

        // =========================== compile bottomup per cpt ===================
        data.cpt_timer.Start();
        for(unsigned int variable = 0; variable < data.counter[1]; variable++){
            #ifdef VERBOSE
            Timer cpt_local_timer;
            cpt_local_timer.Start();
            #endif

            // create  wpbdds for each CPT
            std::queue<bnc::node*> cpt_clauses;
            BayesNode *n = g.get_node(variable);
            bnc::bayesnode_to_wpbdds<false>(&c.manager, cpt_clauses, n, closure, data.support[variable], data.ordering[variable]);
            while(cpt_clauses.size() >= 2){
                bnc::node *bdd1 = cpt_clauses.front();
                cpt_clauses.pop();
                bnc::node *bdd2 = cpt_clauses.front();
                cpt_clauses.pop();

                //bnc::node *bdd = bnc::conjoin<false,true>(&c.manager, partition_id, bdd1, bdd2, data.ordering[variable]);
                bnc::node *bdd = bnc::conjoin<COLLAPSE,DETERMINISM>(&c.manager, partition_id, bdd1, bdd2, data.ordering[variable]);
                cpt_clauses.push(bdd);

                bnc::node::dereference(bdd1);
                bnc::recursive_destroy(&manager, bdd1);

                bnc::node::dereference(bdd2);
                bnc::recursive_destroy(&manager, bdd2);
            }
            if(cpt_clauses.size() != 1)
                throw compiler_debug_exception("CPT clauses not generated");

            bnc::node *cpt = cpt_clauses.front();
            cpt_clauses.pop();

            #ifdef DEBUG
            if(!cpt)
                throw compiler_debug_exception("topdown compiler did not create WPBDD for CPT %u\n", variable);
            #endif
            bnc::add_support(&manager,cpt,data.support[variable]);

            if(COLLAPSE)
                bnc::collapse(&manager, cpt);

            cpts[variable] = cpt;
            data.q_bu.push_async(cpt);

            #ifdef VERBOSE
            cpt_local_timer.Stop();
            printf("Compiled CPT %u/%u (%.3lfs)\n", variable+1, data.counter[0],cpt_local_timer.GetDuration<Timer::Seconds>());
            #endif

        }
        data.cpt_timer.Stop();
    }
    data.total_timer.Stop();
    data.total_timer.Add();


    // compute cumulative size
    synchronous_queue<bnc::node*> tmp;
    while(!data.q_bu.empty())
        tmp.push_async(data.q_bu.pop_async());

    while(!tmp.empty()){
        bnc::node *bdd = tmp.pop_async();
        data.CUMULATIVE_CPT_SIZE += bnc::size(bdd);
        data.CUMULATIVE_OPT_SIZE += bnc::operators(bdd);
        data.q_bu.push_async(bdd);
    }
    manager.cpt_compile_time = data.cpt_timer.GetDuration<Timer::Seconds>();
    manager.cumulative_nodes = data.CUMULATIVE_CPT_SIZE;
    manager.cumulative_operators = data.CUMULATIVE_OPT_SIZE;

    printf("\n");
    printf("Compiled CPTs in              : %.3fs\n", manager.cpt_compile_time);
    printf("Cumulative #nodes per CPT     : %u\n", manager.cumulative_nodes);
    printf("Cumulative #operators per CPT : %u\n\n", manager.cumulative_operators);


    // ======================== conjoin bdds bottom up =================================

    data.q_bu.clear_async();
    data.total_timer.Start();
    data.join_timer.Start();
    #ifdef VERBOSE
    unsigned int TOTAL = data.counter[2];
    #endif

    // compile per partition
    for(unsigned int partition_id = 0; partition_id < partitions.size(); partition_id++){
        partition_t &partition = partitions[partition_id].partition;

        // select cpts per partition
        for(auto it = partition.set.begin(); it != partition.set.end(); it++)
            data.q_bu.push_async(cpts[*it]);

        // conjoin partition variables
        while(data.q_bu.size() > 1){
            data.counter[2]--;
            bnc::node* bdd1 = data.q_bu.pop_async();
            bnc::node* bdd2 = data.q_bu.pop_async();

            #ifdef VERBOSE
            Timer t; t.Start();
            #endif

            #ifdef INTERMEDIATE_DOT
            bnc::write_dot_pdf(&manager, bdd1, "bdd1",NULL,true);
            bnc::write_dot_pdf(&manager, bdd2, "bdd2",NULL,true);
            #endif
            bnc::node *bdd = bnc::conjoin<COLLAPSE,DETERMINISM>(&manager, partition_id, bdd1, bdd2);
            #ifdef INTERMEDIATE_DOT
            bnc::write_dot_pdf(&manager, bdd, "bdd",NULL,true);
            #endif

            //if(COLLAPSE){
            //    support_t support = create_support(manager,bdd1, bdd2);
            //    ordering_t ordering = create_ordering(manager, partition_id, support);
            //    bnc::node *tmp = bnc::conjoin<false,DETERMINISM>(&manager, partition_id, bdd1, bdd2);
            //}

            #ifdef VERBOSE
            t.Stop();
            printf("Conjoining with CPT %2u/%u in %6.3fs\n", TOTAL-(data.counter[2]), TOTAL, t.GetDuration<Timer::Seconds>());
            #endif
            data.q_bu.push_async(bdd);

            bnc::node::dereference(bdd1);
            bnc::destroy_support(&manager, bdd1);
            bnc::recursive_destroy(&manager, bdd1);

            bnc::node::dereference(bdd2);
            bnc::destroy_support(&manager, bdd1);
            bnc::recursive_destroy(&manager, bdd2);
        }
        data.bdds.push_back(data.q_bu.pop_async());
    }

    data.join_timer.Stop();
    data.total_timer.Stop();
    data.total_timer.Add();

    return NULL;
}

template <bool COLLAPSE, bool DETERMINISM, compilation_t COMPILATION_TYPE>
void* compiler::compile_sync(void *v){
    return NULL;
}

pthread_function_t compiler::get_compile(bool SYNC, bool COLLAPSE, bool DETERMINISM, compilation_t COMPILATION_TYPE){
    if(SYNC){
        if(COLLAPSE){
            if(COMPILATION_TYPE == compilation_t::topdown_bottomup){
                if(DETERMINISM)
                    return compiler::compile_sync<true,true,compilation_t::topdown_bottomup>;
                else return compiler::compile_sync<true,false,compilation_t::topdown_bottomup>;
            } else if(COMPILATION_TYPE == compilation_t::topdown){
                if(DETERMINISM)
                    return compiler::compile_sync<true,true,compilation_t::topdown>;
                else return compiler::compile_sync<true,false,compilation_t::topdown>;
            } else {
                if(DETERMINISM)
                    return compiler::compile_sync<true,true,compilation_t::bottomup>;
                else return compiler::compile_sync<true,false,compilation_t::bottomup>;
            }
        } else {
            if(COMPILATION_TYPE == compilation_t::topdown_bottomup){
                if(DETERMINISM)
                    return compiler::compile_sync<false,true,compilation_t::topdown_bottomup>;
                else return compiler::compile_sync<false,false,compilation_t::topdown_bottomup>;
            } else if(COMPILATION_TYPE == compilation_t::topdown){
                if(DETERMINISM)
                    return compiler::compile_sync<false,true,compilation_t::topdown>;
                else return compiler::compile_sync<false,false,compilation_t::topdown>;
            } else {
                if(DETERMINISM)
                    return compiler::compile_sync<false,true,compilation_t::bottomup>;
                else return compiler::compile_sync<false,true,compilation_t::bottomup>;
            }
        }
    } else {
        if(COLLAPSE){
            if(COMPILATION_TYPE == compilation_t::topdown_bottomup){
                if(DETERMINISM)
                    return compiler::compile_async<true,true,compilation_t::topdown_bottomup>;
                else return compiler::compile_async<true,false,compilation_t::topdown_bottomup>;
            } else if(COMPILATION_TYPE == compilation_t::topdown){
                if(DETERMINISM)
                    return compiler::compile_async<true,true,compilation_t::topdown>;
                else return compiler::compile_async<true,false,compilation_t::topdown>;
            } else {
                if(DETERMINISM)
                    return compiler::compile_async<true,true,compilation_t::bottomup>;
                else return compiler::compile_async<true,true,compilation_t::bottomup>;
            }
        } else {
            if(COMPILATION_TYPE == compilation_t::topdown_bottomup){
                if(DETERMINISM)
                    return compiler::compile_async<false,true,compilation_t::topdown_bottomup>;
                else return compiler::compile_async<false,false,compilation_t::topdown_bottomup>;
            } else if(COMPILATION_TYPE == compilation_t::topdown){
                if(DETERMINISM)
                    return compiler::compile_async<false,true,compilation_t::topdown>;
                else return compiler::compile_async<false,false,compilation_t::topdown>;
            } else {
                if(DETERMINISM)
                    return compiler::compile_async<false,true,compilation_t::bottomup>;
                else return compiler::compile_async<false,false,compilation_t::bottomup>;
            }
        }
    }
}

int compiler::compile_wpbdd(){
    bayesgraph &g = manager.get_bayesgraph();
    bn_partitions_t &partitions = manager.partitions;

    #ifndef SIMPLE_ALLOCATION
    if(OPT_PARALLELISM)
        throw compiler_exception("Reusing allocated nodes in not thread safe (add '#define SIMPLE_ALLOCATION' to defines.h and recompile)");

    if(OPT_RESERVE > 0)
        manager.reserve_nodes(OPT_RESERVE);

    //const unsigned int RESERVE = 2*(g.get_nr_literals()+g.get_nr_weights());
    //manager.reserve_nodes(RESERVE);
    #endif

    // hashmap settings
    bnc::computed_table_set_load_factor(&manager,OPT_COMPUTED_TABLE_LOAD_FACTOR);
    bnc::computed_table_rehash(&manager,OPT_COMPUTED_TABLE_BUCKETS);

    const unsigned int NODES = g.size();

    compiler_data_t data;
    data.manager = &manager;
    data.CUMULATIVE_CPT_SIZE = 0;
    data.CUMULATIVE_OPT_SIZE = 0;

    data.id = 0;
    data.c = this;
    if(OPT_COMPILATION_TYPE != compilation_t::topdown){
        data.sat.resize(NODES);
        data.closure.resize(NODES);
        data.ordering.resize(NODES);
        data.support.resize(NODES);
        data.counter.resize(3,NODES);
        data.counter[2]--; // number of conjoins between CPT generated bdds
        for(unsigned int i = NODES; i-- > 0;)
            data.q_init.push_async(i);

        pthread_function_t f = get_compile(false, OPT_COLLAPSE, OPT_DETERMINISM, COMPILATION_TYPE);
        f((void*) &data);

        // bdd is located in queue
        wpbdd = data.bdds;

    } else {
        if(partitions.size() > 1)
            throw compiler_exception("bottom-up compilation of partitions not yet supported");

        data.total_timer.Start();
        // =========================  init sat solver ============================
        sat_t sat;
        sat.set_manager(&manager);

        if(OPT_DETERMINISM)
            sat.init<true>(g);
        else sat.init<false>(g);

        // =========================== compile topdown ===========================
        const ordering_t &ordering = manager.get_ordering(0);
        if(OPT_COLLAPSE){
            if(OPT_DETERMINISM)
                wpbdd.push_back(bnc::solve<true,true>(&manager, sat, ordering));
            else wpbdd.push_back(bnc::solve<true,false>(&manager, sat, ordering));
        } else {
           if(OPT_DETERMINISM)
               wpbdd.push_back(bnc::solve<false,true>(&manager, sat, ordering));
           else wpbdd.push_back(bnc::solve<false,false>(&manager, sat, ordering));
        }

        data.total_timer.Stop();
        data.total_timer.Add();
    }

    manager.total_compile_time = data.total_timer.GetTotal<Timer::Seconds>();
    manager.total_compile_time_ms = data.total_timer.GetTotal<Timer::Milliseconds>();
    manager.join_compile_time = data.join_timer.GetDuration<Timer::Seconds>();
    manager.cpt_compile_time = data.cpt_timer.GetDuration<Timer::Seconds>();
    manager.total_nodes = 0;
    manager.total_operators = 0;
    for(auto it = wpbdd.begin(); it != wpbdd.end(); it++){
        manager.total_nodes += bnc::size(*it);
        manager.total_operators += bnc::operators(*it);
    }
    manager.total_edges = manager.total_nodes * 2;

    print_result();

    return 0;
}


template void* compiler::compile_sync <true,  true, compilation_t::topdown>(void*);
template void* compiler::compile_sync <false, true, compilation_t::topdown>(void*);
template void* compiler::compile_sync <true,  true, compilation_t::bottomup>(void*);
template void* compiler::compile_sync <false, true, compilation_t::bottomup>(void*);
template void* compiler::compile_sync <true,  true, compilation_t::topdown_bottomup>(void*);
template void* compiler::compile_sync <false, true, compilation_t::topdown_bottomup>(void*);

template void* compiler::compile_sync <true,  false, compilation_t::topdown>(void*);
template void* compiler::compile_sync <false, false, compilation_t::topdown>(void*);
template void* compiler::compile_sync <true,  false, compilation_t::bottomup>(void*);
template void* compiler::compile_sync <false, false, compilation_t::bottomup>(void*);
template void* compiler::compile_sync <true,  false, compilation_t::topdown_bottomup>(void*);
template void* compiler::compile_sync <false, false, compilation_t::topdown_bottomup>(void*);

template void* compiler::compile_async<true,  true, compilation_t::topdown>(void*);
template void* compiler::compile_async<false, true, compilation_t::topdown>(void*);
template void* compiler::compile_async<true,  true, compilation_t::bottomup>(void*);
template void* compiler::compile_async<false, true, compilation_t::bottomup>(void*);
template void* compiler::compile_async<true,  true, compilation_t::topdown_bottomup>(void*);
template void* compiler::compile_async<false, true, compilation_t::topdown_bottomup>(void*);

template void* compiler::compile_async<true,  false, compilation_t::topdown>(void*);
template void* compiler::compile_async<false, false, compilation_t::topdown>(void*);
template void* compiler::compile_async<true,  false, compilation_t::bottomup>(void*);
template void* compiler::compile_async<false, false, compilation_t::bottomup>(void*);
template void* compiler::compile_async<true,  false, compilation_t::topdown_bottomup>(void*);
template void* compiler::compile_async<false, false, compilation_t::topdown_bottomup>(void*);


