#include "compiler.h"
#include "multigraph.h"
#include "options.h"
#include "timer.h"

using namespace bnc;

int compiler::compile_multigraph(){
    const bn_partitions_t &kBnPartitions = manager.partitions;

    Timer compile_timer;
    Timer init_timer;
    // compile per partition

    multigraphs.resize(kBnPartitions.size());
    for(unsigned int partition_id = 0; partition_id < kBnPartitions.size(); partition_id++){
        init_timer.Start();
        auto &graph = multigraphs[partition_id];

        if(OPT_BDD_TYPE == bdd_t::multigraph)
            graph.InitSpanningChain(&manager, kBnPartitions[partition_id]);
        else
            graph.InitSpanningTree(&manager, kBnPartitions[partition_id]);

        init_timer.Stop();
        init_timer.Add();

        const SpanningTree& kSpanningTree = graph.GetSpanningTree();
        if(kBnPartitions.size() > 1)
            printf("== Partition %lu ==\n", partition_id);

        printf("Estimated required space: %.3lfMb\n", kSpanningTree.RequiredMb());
        if(OPT_BDD_TYPE == bdd_t::tdmultigraph && !kSpanningTree.IsTree())
            printf("Note: ordering is a chain\n");

        #ifdef DEBUG
        kSpanningTree.PrintStats();
        #endif

        if(!OPT_PARALLELISM){
            compile_timer.Start();
            if(bnc::OPT_ENCODE_STRUCTURE)
                graph.Compile<true>();
            else graph.Compile<false>();
            compile_timer.Stop();
            compile_timer.Add();
        }
    }

    if(OPT_PARALLELISM){
        MultiGraph::ParallelCompile(multigraphs, OPT_WORKERS, compile_timer);
    }

    manager.total_compile_time = compile_timer.GetTotal<Timer::Seconds>()+init_timer.GetTotal<Timer::Seconds>();
    manager.total_compile_time_ms = compile_timer.GetTotal<Timer::Milliseconds>()+init_timer.GetTotal<Timer::Milliseconds>();
    manager.compile_time_ms = compile_timer.GetTotal<Timer::Milliseconds>();
    manager.spanning_tree_time_ms = init_timer.GetTotal<Timer::Milliseconds>();
    manager.total_and_nodes = 0;
    manager.total_or_nodes = 0;
    manager.total_nodes = 0;
    manager.total_operators = 0;
    Timer traverse_timer;
    for(unsigned int partition_id = 0; partition_id < kBnPartitions.size(); partition_id++){
        const auto &kGraph = multigraphs[partition_id];
        const bn_partition_t &kBnPartition = kBnPartitions[partition_id];

        traverse_timer.Start();
        auto size = kGraph.GetSize();
        traverse_timer.Stop();
        traverse_timer.Add();
        if(kBnPartitions.size() > 1)
            printf("[Partition %u] ", partition_id);
        printf("Traversed in: %.3lfms\n", traverse_timer.GetDuration<Timer::Milliseconds>());

        manager.total_nodes += size.nodes;
        manager.total_edges += size.edges;
        manager.total_operators += size.operators;
        manager.total_and_nodes += size.and_nodes;
        manager.total_or_nodes += size.or_nodes;
    }

    printf("Total traversal time: %.3lfms\n", traverse_timer.GetTotal<Timer::Milliseconds>());
    print_result();

    return 0;
}

int compiler::compile_pmultigraph(){
    const bn_partitions_t &kBnPartitions = manager.partitions;

    Timer compile_timer;
    Timer init_timer;
    // compile per partition

    assert(!OPT_PARALLELISM && "parallelism not implemented in combination with use_probabilities");

    pmultigraphs.resize(kBnPartitions.size());
    for(unsigned int partition_id = 0; partition_id < kBnPartitions.size(); partition_id++){
        init_timer.Start();
        auto &graph = pmultigraphs[partition_id];

        if(OPT_BDD_TYPE == bdd_t::multigraph)
            graph.InitSpanningChain(&manager, kBnPartitions[partition_id]);
        else
            graph.InitSpanningTree(&manager, kBnPartitions[partition_id]);

        init_timer.Stop();
        init_timer.Add();

        const SpanningTree& kSpanningTree = graph.GetSpanningTree();
        if(kBnPartitions.size() > 1)
            printf("== Partition %lu ==\n", partition_id);

        printf("Estimated required space: %.3lfMb\n", kSpanningTree.RequiredMb(true));
        if(OPT_BDD_TYPE == bdd_t::tdmultigraph && !kSpanningTree.IsTree())
            printf("Note: ordering is a chain\n");

        #ifdef DEBUG
        kSpanningTree.PrintStats();
        #endif

        compile_timer.Start();
        if(bnc::OPT_ENCODE_STRUCTURE)
            graph.Compile<true>();
        else graph.Compile<false>();
        compile_timer.Stop();
        compile_timer.Add();
        if(kBnPartitions.size() > 1)
            printf("  Partition %lu compilation time: %.3lfms\n", partition_id, compile_timer.GetDuration<Timer::Milliseconds>());
    }

    manager.total_compile_time = compile_timer.GetTotal<Timer::Seconds>()+init_timer.GetTotal<Timer::Seconds>();
    manager.total_compile_time_ms = compile_timer.GetTotal<Timer::Milliseconds>()+init_timer.GetTotal<Timer::Milliseconds>();
    manager.compile_time_ms = compile_timer.GetTotal<Timer::Milliseconds>();
    manager.spanning_tree_time_ms = init_timer.GetTotal<Timer::Milliseconds>();
    manager.total_and_nodes = 0;
    manager.total_or_nodes = 0;
    manager.total_nodes = 0;
    manager.total_operators = 0;
    Timer traverse_timer;
    for(unsigned int partition_id = 0; partition_id < kBnPartitions.size(); partition_id++){
        const auto &kGraph = pmultigraphs[partition_id];
        const bn_partition_t &kBnPartition = kBnPartitions[partition_id];

        traverse_timer.Start();
        auto size = kGraph.GetSize();
        traverse_timer.Stop();
        traverse_timer.Add();
        if(kBnPartitions.size() > 1)
            printf("[Partition %u] ", partition_id);
        printf("Traversed in: %.3lfms\n", traverse_timer.GetDuration<Timer::Milliseconds>());

        manager.total_nodes += size.nodes;
        manager.total_edges += size.edges;
        manager.total_operators += size.operators;
        manager.total_and_nodes += size.and_nodes;
        manager.total_or_nodes += size.or_nodes;
    }
    printf("Total traversal time: %.3lfms\n", traverse_timer.GetTotal<Timer::Milliseconds>());
    print_result();

    return 0;
}

