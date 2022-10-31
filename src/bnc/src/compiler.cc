#include <vector>
#include "compiler.h"
#include "domainclosure.h"
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
#include "ordering.h"
#include <limits>
#include "architecturebound.h"

using namespace std;
using namespace bnc;

#ifdef TIMED
const unsigned int BARRIERS = 10;
pthread_barrier_t barrier[BARRIERS];
#endif
robust_cond_mutex mutex;
unsigned int MERGES;

compiler::compiler(){
    BDD_TYPE = bdd_t::none;
    COMPILATION_TYPE = compilation_t::none;
}

compiler::~compiler(){
    if(wpbdd.size() != 0){
        for(auto it = wpbdd.begin(); it != wpbdd.end(); it++){
            bnc::node::dereference(*it);
            bnc::recursive_destroy(&manager, *it);
        }
        wpbdd.clear();
    }
}

void compiler::write(file_t ft){
    switch(ft){
        case ORDERING:
            {
                if(manager.partitions.size() == 1)
                    manager.partitions[0].ordering.write(files.get_filename(ORDERING));
                else {
                    for(unsigned int i = 0; i < manager.partitions.size(); i++)
                        manager.partitions[i].ordering.write(files.get_filename(ORDERING,stringf("%u",i)));
                }
            }
            break;
        case ELIM_ORDERING:
            {
                if(manager.partitions.size() == 1)
                    manager.partitions[0].variable_ordering.write(files.get_filename(ELIM_ORDERING),manager.get_bayesnet());
                else {
                    for(unsigned int i = 0; i < manager.partitions.size(); i++)
                        manager.partitions[i].variable_ordering.write(files.get_filename(ELIM_ORDERING,stringf("%u",i)),manager.get_bayesnet());
                }
            }
            break;
        case VARIABLE_ORDERING:
            {
                if(manager.partitions.size() == 1)
                    manager.partitions[0].variable_ordering.write(files.get_filename(VARIABLE_ORDERING));
                else {
                    for(unsigned int i = 0; i < manager.partitions.size(); i++)
                        manager.partitions[i].variable_ordering.write(files.get_filename(VARIABLE_ORDERING,stringf("%u",i)));
                }
            }
            break;
        case UAI:
            {
                bayesgraph &g = manager.get_bayesgraph();
                g.WriteUAI();
            }
            break;
        case PARTITION:
            {
                std::string filename = files.get_filename_c(PARTITION);
                manager.partitions.write(filename, manager.get_bayesnet());
            }
            break;
        case PSEUDO_ORDERING:
            {
                if(manager.partitions.size() == 1){
                    std::string filename = files.get_filename(PSEUDO_ORDERING);
                    bn_partition_t &bn_partition = manager.partitions[0];
                    PseudoTree &pseudotree = bn_partition.pseudotree;
                    if(!pseudotree.Empty()){
                        pseudotree.Write(filename);
                    } else if(bn_partition.variable_ordering.size() != 0){
                        fprintf(stderr, "Warning: Pseudo tree is empty, building a chain\n");
                        pseudotree.SetChain(manager.get_bayesnet(), bn_partition.partition, bn_partition.variable_ordering);
                        pseudotree.Write(filename);
                    } else fprintf(stderr, "Error: Pseudo tree is empty\n");

                } else {
                    for(unsigned int i = 0; i < manager.partitions.size(); i++){
                        std::string filename = files.get_filename(PSEUDO_ORDERING,stringf("%u",i));
                        bn_partition_t &bn_partition = manager.partitions[i];
                        PseudoTree &pseudotree = bn_partition.pseudotree;
                        if(!pseudotree.Empty()){
                            pseudotree.Write(filename);
                        } else if(bn_partition.variable_ordering.size() != 0){
                            fprintf(stderr, "Warning: Pseudo tree is empty, building a chain\n");
                            pseudotree.SetChain(manager.get_bayesnet(), bn_partition.partition, bn_partition.variable_ordering);
                            pseudotree.Write(filename);
                        } else fprintf(stderr, "Error: Pseudo tree is empty\n");
                    }
                }
            }
            break;

        case DOT:
            if(BDD_TYPE == bdd_t::wpbdd && !wpbdd.empty()){
                bnc::write_dot(&manager, wpbdd);
            } else if(BDD_TYPE == bdd_t::multigraph || BDD_TYPE == bdd_t::tdmultigraph){
                if(OPT_USE_PROBABILITY)
                    MultiGraphProbability::DumpDot(pmultigraphs);
                else
                    MultiGraph::DumpDot(multigraphs);
            } else
                throw compiler_exception("Cannot write dot without bdd");
            break;

        case SPANNING:
            {
                bn_partitions_t &bn_partitions = manager.partitions;
                if(bn_partitions.size() == 0 || bn_partitions[0].variable_ordering.empty())
                    compiler_exception("Missing variable ordering");

                if(OPT_USE_PROBABILITY){
                    if(pmultigraphs.empty())
                        pmultigraphs.resize(bn_partitions.size());
                } else {
                    if(multigraphs.empty())
                        multigraphs.resize(bn_partitions.size());
                }

                for(unsigned int i = 0; i < bn_partitions.size(); i++){
                    bn_partition_t &bn_partition = bn_partitions[i];
                    if(OPT_USE_PROBABILITY){
                        if(!pmultigraphs[i].HasSpanningTree()){ //
                            auto &graph = pmultigraphs[i];
                            if(bn_partition.pseudotree.Empty())
                                graph.InitSpanningChain(&manager, bn_partition);
                            else graph.InitSpanningTree(&manager, bn_partition);
                        }

                        const SpanningTree& kSpanningTree = pmultigraphs[i].GetSpanningTree();
                        kSpanningTree.Write((bn_partitions.size()>1?i:-1));
                    } else {
                        if(!multigraphs[i].HasSpanningTree()){ //
                            auto &graph = multigraphs[i];
                            if(bn_partition.pseudotree.Empty())
                                graph.InitSpanningChain(&manager, bn_partition);
                            else graph.InitSpanningTree(&manager, bn_partition);
                        }

                        const SpanningTree& kSpanningTree = multigraphs[i].GetSpanningTree();
                        kSpanningTree.Write((bn_partitions.size()>1?i:-1));
                    }
                }

            }
            break;
        case AC:
            if (OPT_BDD_TYPE == bdd_t::wpbdd){
                if(BDD_TYPE == bdd_t::wpbdd && !wpbdd.empty())
                    bnc::write_bdd(&manager, wpbdd);
                else throw compiler_exception("Cannot write circuit without bdd");
            } else if (OPT_BDD_TYPE == bdd_t::multigraph || BDD_TYPE == bdd_t::tdmultigraph) {
                if(OPT_USE_PROBABILITY)
                    MultiGraphProbability::DumpDD(pmultigraphs);
                else
                    MultiGraph::DumpDD(multigraphs);
            } else throw compiler_exception("No write function available of diagram type");
            break;
        case MAPPING:
            {
                const bayesgraph &g = manager.get_bayesgraph();
                if(g.defined()){
                    g.write();
                } else throw compiler_exception("Bayesian Network mapping missing");
            }
            break;
        case COMPOSITION_ORDERING:
            {
                const ordering_t &ordering = manager.composition_ordering;
                std::string filename = files.get_filename(COMPOSITION_ORDERING);
                ordering.write(filename);
            }
            break;
        default:
            throw compiler_exception("Undefined write function");
            break;
    }
}

void compiler::set_compilation_type(compilation_t type){
    COMPILATION_TYPE = type;
}

void compiler::set_bdd_type(bdd_t type){
    BDD_TYPE = type;
}

void compiler::encode(bayesnet *bn){
    if(BDD_TYPE == bdd_t::none)
        throw compiler_exception("BDD type is required to be set");

    manager.encode(bn);
}

void compiler::load(){
    if(BDD_TYPE == bdd_t::none)
        throw compiler_exception("BDD type is required to be set");

    manager.load(BDD_TYPE);

    if(OPT_INFO)
        info();
}

void compiler::info(){
    bayesgraph &g = manager.get_bayesgraph();
    if(g.defined()){
        printf("BN Variables  : %u\n", g.get_nr_variables());
        printf("Literals      : %u\n", g.get_nr_literals());
        printf("Probabilities : %u\n", g.get_nr_weights());
    }
}

size_t compiler::get_nr_partitions(){
    return manager.partitions.size();
}


int perm(unsigned int n) {
    std::vector<int> list(n);
    for(unsigned int i = 0; i < n; i++)
        list[i] = i;


    do {
        for(unsigned int i = 0; i < n; i++)
            printf(" %d", list[i]);
        printf("\n");

    } while (next_permutation(list.begin(), list.end()));
	return 0;
}


void comb(int N, int K){
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    // print integers and permute bitmask
    do {
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) printf(" %d", i);
        }
        printf("\n");
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

void compiler::compile(){
    switch(BDD_TYPE){
        case bdd_t::multigraph:
        case bdd_t::tdmultigraph:
            if(OPT_USE_PROBABILITY)
                compile_pmultigraph();
            else
                compile_multigraph();
            break;
        case bdd_t::wpbdd:
            compile_wpbdd();
            break;
        default:
           throw compiler_exception("BDD type is required to be set");
    }
}

void compiler::print_result(){
    printf("\nFINAL RESULT:\n\n");

    if(OPT_BDD_TYPE == bdd_t::multigraph || OPT_BDD_TYPE == bdd_t:: tdmultigraph){
        printf("    Spanning tree     : %.3fms\n", manager.spanning_tree_time_ms);
        printf("    Compilation       : %.3fms\n", manager.compile_time_ms);
        printf("    Total Or #nodes   : %lu\n", manager.total_or_nodes);
        printf("    Total And #nodes  : %lu\n", manager.total_and_nodes);
    } else {
        printf("    Compiled CPTs in  : %.3fs\n", manager.cpt_compile_time);
        printf("    Conjoined CPTs in : %.3fs\n", manager.join_compile_time);
    }
    printf("    Total time        : %.3fs\n", manager.total_compile_time);
    printf("    Total time        : %.3fms\n", manager.total_compile_time_ms);
    printf("\n");
    printf("    Total #nodes      : %lu\n", manager.total_nodes);
    printf("    Total #edges      : %lu\n", manager.total_edges);
    printf("    Total #operators  : %lu\n", manager.total_operators);
    printf("\n");
}

bnc::node* compiler::get_wpbdd(){
    return wpbdd[0];
}

