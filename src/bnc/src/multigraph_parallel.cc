#include "multigraph.h"

#include <pthread.h>
#include <thread>
#include "options.h"
#include "exceptions.h"
#include "timer.h"
#include "bayesgraph.h"
#include <unistd.h>

#undef likely
#undef unlikely

#define OPTION 2

#if OPTION == 1
#include "lockfreequeue.h"
#elif OPTION == 2
#include <boost/lockfree/spsc_queue.hpp>
#elif OPTION == 3
#include <boost/lockfree/queue.hpp>
#elif OPTION == 4
#include "boostqueue.h"
#endif

namespace bnc {

#define MAX_THREADS 64
#define MAX_PARTITIONS 12

MultiGraph::Cache* g_cache[MAX_PARTITIONS];
const SpanningTree* g_kSpanningTree[MAX_PARTITIONS];
const uint32_t *g_kDimensions;
bayesgraph const * g_kBayesgraph;
unsigned int g_kThreadPoolSize;

struct Task {
    enum Type { kCompileOrNode = 0, kCompileOrNodes = 1, kCompilePartition = 2, kFinished = 3 };

    Type type;
    const SpanningTree * kSpanningTree;
    const SpanningTreeNode * kSpanningNode;
    const SpanningTreeNode * kChildSpanningNode;
    MultiGraph::Cache *cache;
    unsigned int id;
};

unsigned int thread_id;

#if OPTION == 1
LockFreeQueue<Task, 1024> task_queue[MAX_THREADS];
#elif OPTION == 2
boost::lockfree::spsc_queue<Task, boost::lockfree::capacity<1024> > task_queue[MAX_THREADS];
#elif OPTION == 3
boost::lockfree::queue<Task> task_queue(4096);
#elif OPTION == 4
concurrent_queue<Task> task_queue;
#endif


pthread_barrier_t start_barrier;

void compile_layer_and(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode, MultiGraph::Cache &cache){
    assert(false && "not implemented");
}

void scheduler_compile_layer_and(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode, MultiGraph::Cache &cache){
    assert(false && "not implemented");
    if(OPT_PARALLEL_LEVEL == 3){
        // create tasks that compute nodes in each layer in parallel
        //XAry::Dimension dimension;
        //dimension.resize(kSpanningNode->spanning.size());
        //XAry::SetDimension(dimension.begin(), kSpanningNode->spanning.begin(), kSpanningNode->spanning.end(), g_kDimensions);

    } else {
        printf("Compile layer %lu/%lu (AND)\n", kSpanningNode->index, kSpanningTree->GetSize());
    }
}

void set_edges(
        MultiGraph::Node * const node,
        const XAry::Bit * const xary_parent_ctr,
        const SpanningTreeNode * const kChildSpanningNode,
        MultiGraph::Cache &cache){

    const XAry::Dim * const kChildDimension = &(kChildSpanningNode->xary.dimension[0]);
    const unsigned int kChildCtrSize = kChildSpanningNode->spanning.size();
    XAry::Bit xary_child_ctr[kChildCtrSize];
    XAry::Set(xary_child_ctr, &(xary_child_ctr[kChildCtrSize]), xary_parent_ctr, &(kChildSpanningNode->xary.map[0]));

    // ============ set edges
    if(kChildSpanningNode->xary.expand){
        const auto kChildId = XAry::GetDecimal(xary_child_ctr, kChildDimension, kChildCtrSize);
        auto * const kChildNode = cache[kChildSpanningNode->index].GetOrNode(kChildId);

        auto *edge = &(node->edges[0]);
        const auto * kEdgeEnd = &(node->edges[node->size]);
        while(edge != kEdgeEnd){
            edge->to = kChildNode;
            ++edge;
        }
    } else {
        auto &child_ctr_parent_value = xary_child_ctr[kChildSpanningNode->xary.pos];

        auto *edge = &(node->edges[0]);
        const auto * kEdgeEnd = &(node->edges[node->size]);
        while(edge != kEdgeEnd){
            const auto kChildId = XAry::GetDecimal(xary_child_ctr, kChildDimension, kChildCtrSize);
            auto * const kChildNode = cache[kChildSpanningNode->index].GetOrNode(kChildId);

            edge->to = kChildNode;
            ++child_ctr_parent_value;
            ++edge;
        }
    }
}

void set_weights(
        MultiGraph::Node * const node,
        const XAry::Bit * const xary_parent_ctr,
        const SpanningTree * const kSpanningTree,
        const SpanningTreeNode * const kSpanningNode,
        MultiGraph::Cache &cache){

    // ============ add weights
    // collect cpt info
    const auto *kMessages = &(kSpanningNode->messages[0]);
    const auto kNrWeights = kSpanningNode->messages.size();

    const unsigned int* kCptDimension[kNrWeights];
    unsigned int kCptNrVariables[kNrWeights];
    const unsigned int* kCptWeights[kNrWeights];

    for(unsigned int i = 0; i < kNrWeights; i++){
        auto *kCpt = g_kBayesgraph->get_node(kMessages[i]);
        kCptDimension[i] = &(kCpt->GetDimensions()[0]);
        kCptNrVariables[i] = kCpt->GetNrVariables();
        kCptWeights[i] = &(kCpt->GetWeights()[0]);
    }

    MultiGraph::Edge *edge = &(node->edges[0]);
    const auto * kEdgeEnd = &(node->edges[node->size]);
    unsigned int value = 0;
    while(edge != kEdgeEnd){
        auto *weights = &(edge->weights[0]);
        for(unsigned int i = 0; i < kNrWeights; i++){
            const unsigned int kCptId = kMessages[i];

            XAry::Bit xary_cpt_ctr[kCptNrVariables[i]];
            XAry::Set(xary_cpt_ctr, &(xary_cpt_ctr[kCptNrVariables[i]]), xary_parent_ctr, &(kSpanningTree->cpt_map_[kCptId][0]));
            xary_cpt_ctr[kSpanningTree->cpt_variable_pos_[kCptId]] = value;

            const unsigned int kWeightIndex = XAry::GetDecimal(xary_cpt_ctr, kCptDimension[i], kCptNrVariables[i]);
            const auto kWeight = kCptWeights[i][kWeightIndex];
            if(kWeight == 0){ // determinism
                edge->size = 1;
                weights[0] = 0;
                edge->to = cache.GetTerminal();
                break;
            } else {
                weights[i] = kWeight;
            }
        }
        ++value;
        ++edge;
    }
}

void create_node_or(
        const SpanningTree * const kSpanningTree,
        const SpanningTreeNode * const kSpanningNode,
        const SpanningTreeNode * const kChildSpanningNode,
        MultiGraph::Cache &cache,
        const unsigned int kParentId){

    //printf("create node %lu:%lu\n", kSpanningNode->variable, kParentId);

    const XAry::Dim * const kParentDimension = &(kSpanningNode->xary.dimension[0]);
    const unsigned int kParentCtrSize = kSpanningNode->spanning.size();
    XAry::Bit xary_parent_ctr[kParentCtrSize];
    XAry::SetDecimal(xary_parent_ctr, kParentDimension, kParentCtrSize, kParentId);

    // ============ create node
    auto *node = cache[kSpanningNode->index].CreateOrNode(kParentId);
    node->variable = kSpanningNode->variable;

    // add edges
    set_edges(node, xary_parent_ctr, kChildSpanningNode, cache);

    // add weights
    set_weights(node, xary_parent_ctr, kSpanningTree, kSpanningNode, cache);

}

void compile_layer_or(const SpanningTree * const kSpanningTree, const SpanningTreeNode * const kSpanningNode, MultiGraph::Cache &cache){
    const SpanningTreeNode * const kChildSpanningNode = &(kSpanningNode->children[0]);
    const auto kNrParents = kSpanningNode->GetOrUpperbound();
    for(unsigned int parent_id = 0; parent_id < kNrParents; parent_id++)
        create_node_or(kSpanningTree, kSpanningNode, kChildSpanningNode, cache, parent_id);
}

void scheduler_compile_layer_or(const SpanningTree * const kSpanningTree, const SpanningTreeNode * const kSpanningNode, MultiGraph::Cache &cache){
    const SpanningTreeNode * const kChildSpanningNode = &(kSpanningNode->children[0]);

    if(OPT_PARALLEL_LEVEL == 3){
        Task task;
        task.type = Task::Type::kCompileOrNode;
        task.kSpanningTree = kSpanningTree;
        task.kSpanningNode = kSpanningNode;
        task.kChildSpanningNode = kChildSpanningNode;
        task.cache = &cache;

        const unsigned int kNrThreads = g_kThreadPoolSize;


        const auto kNrParents = kSpanningNode->GetOrUpperbound();
        for(unsigned int parent_id = 0; parent_id < kNrParents; parent_id++){
            task.id = parent_id;
            #if OPTION == 4 || OPTION == 3
            task_queue.push(task);
            #else
            while(!task_queue[thread_id].push(task)){
                #ifdef DEBUG
                printf("Queue %lu is full\n", thread_id);
                #endif
                thread_id = (thread_id+1) % kNrThreads;
            }
            thread_id = (thread_id+1) % kNrThreads;
            #endif
        }

    } else {
        compile_layer_or(kSpanningTree, kSpanningNode, cache);
    }
}

void compile_partition(const SpanningTree *kSpanningTree, MultiGraph::Cache &cache){
    assert(kSpanningTree->GetSize() == cache.size());

    const auto kSpanningNodes = kSpanningTree->GetSpanningNodes();
    const auto kSpanningNodeEnd = kSpanningNodes.end();
    auto kSpanningNode_it = kSpanningNodes.begin();
    while(kSpanningNode_it != kSpanningNodeEnd){
        const SpanningTreeNode *kSpanningNode = *kSpanningNode_it;
        if(!kSpanningNode->IsRoot() && !kSpanningNode->IsTerminal()){
            if(kSpanningNode->IsAndLayer())
                compile_layer_and(kSpanningTree, kSpanningNode, cache);
            else compile_layer_or(kSpanningTree, kSpanningNode, cache);
        }
        ++kSpanningNode_it;
    }
}


void scheduler_compile_partition(const SpanningTree *kSpanningTree, MultiGraph::Cache &cache){
    if(OPT_PARALLEL_LEVEL == 2){
        // create tasks that compute layers in parallel

        Task task;
        task.type = Task::Type::kCompileOrNodes;
        task.kSpanningTree = kSpanningTree;
        task.cache = &cache;

        const unsigned int kNrThreads = g_kThreadPoolSize;

        const auto kSpanningNodes = kSpanningTree->GetSpanningNodes();
        const auto kSpanningNodeEnd = kSpanningNodes.end();
        auto kSpanningNode_it = kSpanningNodes.begin();
        while(kSpanningNode_it != kSpanningNodeEnd){
            const SpanningTreeNode *kSpanningNode = *kSpanningNode_it;
            if(!kSpanningNode->IsRoot() && !kSpanningNode->IsTerminal()){
                if(kSpanningNode->IsAndLayer()){
                    assert(false && "not implemented");
                } else {
                    task.type = Task::Type::kCompileOrNodes;
                    task.kSpanningNode = kSpanningNode;
                    #if OPTION == 4 || OPTION == 3
                    task_queue.push(task);
                    #else
                    while(!task_queue[thread_id].push(task)){
                        #ifdef DEBUG
                        printf("Queue %lu is full\n", thread_id);
                        #endif
                        thread_id = (thread_id+1) % kNrThreads;
                    }
                    thread_id = (thread_id+1) % kNrThreads;
                    #endif
                }
            }
            ++kSpanningNode_it;
        }

    } else {
        assert(kSpanningTree->GetSize() == cache.size());

        const auto kSpanningNodes = kSpanningTree->GetSpanningNodes();
        const auto kSpanningNodeEnd = kSpanningNodes.end();
        auto kSpanningNode_it = kSpanningNodes.begin();
        while(kSpanningNode_it != kSpanningNodeEnd){
            const SpanningTreeNode *kSpanningNode = *kSpanningNode_it;
            if(!kSpanningNode->IsRoot() && !kSpanningNode->IsTerminal()){
                if(kSpanningNode->IsAndLayer())
                    scheduler_compile_layer_and(kSpanningTree, kSpanningNode, cache);
                else scheduler_compile_layer_or(kSpanningTree, kSpanningNode, cache);
            }
            ++kSpanningNode_it;
        }
    }
}

void* pthread_worker(void *in){
    const size_t kThreadId = (size_t)in;

    pthread_barrier_wait(&start_barrier);

    #if OPTION == 4 || OPTION == 3
    auto &queue = task_queue;
    #else
    auto &queue = task_queue[kThreadId];
    #endif

    Task task;
    while(true){
        while(!queue.pop(task));

        switch(task.type){
            case Task::Type::kFinished:
                return 0;

            case Task::Type::kCompileOrNode:
                create_node_or(task.kSpanningTree, task.kSpanningNode, task.kChildSpanningNode, *task.cache, task.id);
                break;

            case Task::Type::kCompileOrNodes:
                compile_layer_or(task.kSpanningTree, task.kSpanningNode, *task.cache);
                break;

            case Task::Type::kCompilePartition:
                compile_partition(task.kSpanningTree, *task.cache);
                break;
            default:
                assert(false && "Unknown task option");
        }
    }



    return 0;
}

void MultiGraph::ParallelScheduler(std::vector<MultiGraph> &multigraphs, const unsigned int kNrThreads){
    if(OPT_PARALLEL_LEVEL == 1){
        Task task;
        task.type = Task::Type::kCompilePartition;

        const unsigned int kNrThreads = g_kThreadPoolSize;
        for(unsigned int i = 0; i < multigraphs.size(); i++){
            task.kSpanningTree = &(multigraphs[i].spanningtree_);
            task.cache = &(multigraphs[i].cache_);
            #if OPTION == 4 || OPTION == 3
            task_queue.push(task);
            #else
            while(!task_queue[thread_id].push(task)){
                #ifdef DEBUG
                printf("Queue %lu is full\n", thread_id);
                #endif
                thread_id = (thread_id+1) % kNrThreads;
            }
            thread_id = (thread_id+1) % kNrThreads;
            #endif
        }
    } else {
        for(unsigned int i = 0; i < multigraphs.size(); i++){
            scheduler_compile_partition(&(multigraphs[i].spanningtree_), multigraphs[i].cache_);
        }
    }
}

void MultiGraph::ParallelCompile(std::vector<MultiGraph> &multigraphs, const unsigned int kNrThreads, Timer &timer){
    assert(!OPT_ENCODE_STRUCTURE && "Parallel compilation with local structure is not supported");
    assert(multigraphs.size() > 0 && multigraphs.size() <= MAX_PARTITIONS);
    assert(OPT_PARALLEL_LEVEL >= 0 && OPT_PARALLEL_LEVEL <= 3 && "Unknown parallel level");
    for(unsigned int i = 0; i < multigraphs.size(); i++)
        assert(!multigraphs[i].GetSpanningTree().IsTree() && "trees not yet supported");

    const unsigned int kMaxThreads = std::thread::hardware_concurrency();
    const unsigned int kMaxThreadPoolSize = (kMaxThreads-1<MAX_THREADS?kMaxThreads-1:MAX_THREADS);
    const unsigned int kThreadPoolSize = (kNrThreads==0?kMaxThreads-1:kNrThreads);
    thread_id = 0;

    if(kThreadPoolSize > kMaxThreadPoolSize)
        throw compiler_parallel_exception("Trying to start thread pool with %lu threads (max: %lu)",kThreadPoolSize, kMaxThreadPoolSize);

    printf("Threadpool size : %lu\n", kThreadPoolSize);
    printf("Parallel level  : %lu\n", OPT_PARALLEL_LEVEL);
    //printf("Lock free queue : %s\n",(task_queue.is_lock_free()?"yes":"no"));
    printf("Queue           : ");
    #if OPTION == 1
    printf("lockfree queue\n");
    #elif OPTION == 2
    printf("lockfree boost single producer/consumer\n");
    #elif OPTION == 3
    printf("lockfree boost multi producer/consumer\n");
    #elif OPTION == 4
    printf("locking boost multi producer/consumer\n");
    #endif

    //for(unsigned int i = 0; i < MAX_THREADS; i++)
    //    task_queue[i].reserve(1024);

    // set global variables
    g_kThreadPoolSize = kThreadPoolSize;
    g_kBayesgraph = multigraphs[0].g_;
    g_kDimensions = multigraphs[0].bn_->get_states();
    for(unsigned int i = 0; i < multigraphs.size(); i++){
        g_cache[i] = &(multigraphs[i].cache_);
        g_kSpanningTree[i] = &(multigraphs[i].spanningtree_);
    }
    // allocate
    timer.Start();
    for(unsigned int i = 0; i < multigraphs.size(); i++)
        multigraphs[i].AllocateCache(false);
    timer.Stop();
    timer.Add();

    pthread_barrier_init(&start_barrier, NULL, kThreadPoolSize+1);
    std::array<pthread_t,MAX_THREADS> threads;
    for(unsigned int i = 0; i < kThreadPoolSize; i++){
        if(pthread_create(&(threads[i]), NULL, pthread_worker, (void*)i)){
            fprintf(stderr, "Error creating thread %u of %u\n", i+1, kThreadPoolSize);
            exit(1);
        }
    }
    pthread_barrier_wait(&start_barrier);
    usleep(1000000);

    // start thread pool

    timer.Start();
    MultiGraph::ParallelScheduler(multigraphs,kThreadPoolSize);

    Task task;
    task.type = Task::Type::kFinished;
    for(unsigned int i = 0; i < kThreadPoolSize; i++){
        #if OPTION == 1 || OPTION == 2
        while(!task_queue[i].push(task));
        #else
        task_queue.push(task);
        #endif
    }

    // wait until pool is finished
    for(unsigned int i = 0; i < kThreadPoolSize; i++){
        if(pthread_join(threads[i],NULL)){
            fprintf(stderr, "Error canceling thread %u of %u\n", i+1, kThreadPoolSize);
            exit(1);
        }
    }
    timer.Stop();
    timer.Add();

    printf("Threading time: %.3lfms\n", timer.GetDuration<Timer::Milliseconds>());

    // store roots
    for(unsigned int i = 0; i < multigraphs.size(); i++)
        multigraphs[i].root_ = multigraphs[i].cache_.GetRoot();

    pthread_barrier_destroy(&start_barrier);
}

} // namespace bnc

