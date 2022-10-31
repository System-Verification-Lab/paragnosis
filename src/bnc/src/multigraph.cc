#include "multigraph.h"
#include "options.h"
#include "bound.h"
#include "exceptions.h"
#include <unordered_map>
#include <limits>
#include "timer.h"


namespace bnc {

#define MULT_ASSIGN(a, b) __builtin_umull_overflow (a, b, &a)
#define ADD_ASSIGN(a, b)  __builtin_uaddl_overflow (a, b, &a)

template <class T>
inline bool add_assign_or_fail(Bound::BoundSize &a, const T b){
    return ADD_ASSIGN(a, b);
}

template <class T>
inline bool mult_assign_or_fail(Bound::BoundSize &a, const T b){
    return MULT_ASSIGN(a, b);
}


MultiGraph::MultiGraph(){
    manager_ = NULL;
    bn_ = NULL;
    g_ = NULL;
    terminal_ = NULL;
    root_ = NULL;
    initialized_ = false;
}

void MultiGraph::SetManager(bnc::manager_t *manager){
    manager_ = manager;
    bn_ = manager_->get_bayesnet();
    g_ = &manager_->get_bayesgraph();
    spanningtree_.SetBayesnet(bn_);
    spanningtree_.SetBayesgraph(g_);
}

void MultiGraph::RuntimeInit(const partition_t &kPartition){
    const unsigned int kNrBnVariables = bn_->get_nr_variables();

    cpt_ctr_.resize(kNrBnVariables);
    cpt_weights_.resize(kNrBnVariables);
    for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++){
        BayesNode *n = g_->get_node(*it);
        cpt_weights_[*it] = &(n->GetWeights()[0]);
        cpt_ctr_[*it].SetDimension(n->GetDimensions());
    }
}

template <bool STRUCTURE>
void MultiGraph::Compile(){
    assert(initialized_);
    Compile<STRUCTURE>(spanningtree_);
}

bool Expand(const SpanningTreeNode *kSpanningNode, const Variable kVariable, bool *expands){
    auto kFirst = kSpanningNode->children.begin();
    const auto kLast = kSpanningNode->children.end();

    auto expand_it = &(expands[0]);
    bool expand = false;
    while(kFirst != kLast){
        *expand_it = kFirst->spanning.find(kVariable) == kFirst->spanning.end();
        if(*expand_it)
            expand = true;

        ++kFirst;
        ++expand_it;
    }
    return expand;
}

void MultiGraph::AddWeights(MultiGraph::Node *node, const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode, const XAry &xary, const unsigned int kDimension){
    const auto &kMessages = kSpanningNode->messages;

    Edge *edge = &(node->edges[0]);
    const auto kNrWeights = kMessages.size();
    for(unsigned int value = 0; value < kDimension; value++){
        auto *weights = &(edge->weights[0]);
        edge->size = kNrWeights;
        for(unsigned int i = 0; i < kNrWeights; i++){
            const auto kCpt = kMessages[i];

            XAry &ctr = cpt_ctr_[kCpt];
            ctr.Set(xary,kSpanningTree->cpt_map_[kCpt]);
            ctr[kSpanningTree->cpt_variable_pos_[kCpt]] = value;

            const unsigned int kWeightIndex = ctr.GetDecimal();
            const auto kWeight = cpt_weights_[kCpt][kWeightIndex];
            if(kWeight == 0){ // determinism
                edge->size = 1;
                weights[0] = 0;
                edge->to = terminal_;
                break;
            } else {
                weights[i] = kWeight;
            }
        }
        ++edge;
    }
}

void MultiGraph::AddEdges(MultiGraph::Node **to, const XAry &xary, XAry &xary_child, const SpanningTreeNode *kChildSpanningNode, const unsigned int kDimension, const bool kExpand){
    #ifdef DEBUG
    assert((xary_child.Max(kChildSpanningNode->xary.restrict_dimension) == 1 || xary_child.Max(kChildSpanningNode->xary.restrict_dimension) == kDimension) && "Parent has too many children");
    #endif

    xary_child.Set(xary,kChildSpanningNode->xary.map);
    unsigned int value = 0; // kVariable is equal to 'value'
    do {
        auto * const kChildNode = cache_[kChildSpanningNode->index].map[xary_child.GetDecimal()];
        // add edge(s) kParentId -> kChildId

        do {
            to[value] = kChildNode;
            ++value;
        } while (kExpand && value < kDimension);
    } while(xary_child.RestrictedIncrement(kChildSpanningNode->xary.restrict_dimension));
}

void MultiGraph::AddEdges(MultiGraph::Node *node, const XAry &xary, XAry &xary_child, const SpanningTreeNode *kChildSpanningNode, const unsigned int kDimension, const bool kExpand){
    #ifdef DEBUG
    assert((xary_child.Max(kChildSpanningNode->xary.restrict_dimension) == 1 || xary_child.Max(kChildSpanningNode->xary.restrict_dimension) == kDimension) && "Parent has too many children");
    #endif

    xary_child.Set(xary,kChildSpanningNode->xary.map);
    unsigned int value = 0; // kVariable is equal to 'value'
    do {
        auto * const kChildNode = cache_[kChildSpanningNode->index].map[xary_child.GetDecimal()];

        // add edge(s) kParentId -> kChildId
        do {
            node->edges[value].to = kChildNode;
            ++value;
        } while (kExpand && value < kDimension);
    } while(xary_child.RestrictedIncrement(kChildSpanningNode->xary.restrict_dimension));
}

template<bool STRUCTURE>
void MultiGraph::CreateOrNodes(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode){
    const SpanningTreeNode *kChildSpanningNode = &(kSpanningNode->children[0]);

    // prepare counters
    XAry xary;
    xary.SetDimension(bn_,kSpanningNode->spanning);
    xary.SetZero();

    XAry xary_child;
    xary_child.SetDimension(bn_,kChildSpanningNode->spanning);

    // store common variables
    const auto kVariable = kSpanningNode->variable;
    const auto kDimension = kSpanningNode->dimension;
    const auto kIndex = kSpanningNode->index;
    const bool kExpand = kChildSpanningNode->spanning.find(kVariable) == kChildSpanningNode->spanning.end();

    // layer caching
    CacheLayer &cache = cache_[kIndex];
    auto *map = &(cache.map[0]);
    auto *node = cache.CreateOrNode();

    // create level connections
    unsigned int parent_id = 0;
    do {
        // create OR node
        node->variable = kVariable;
        AddEdges(node, xary, xary_child, kChildSpanningNode, kDimension, kExpand);
        AddWeights(node, kSpanningTree, kSpanningNode, xary, kDimension);

        // canonical check and caching
        auto *canonical = Canonical<STRUCTURE>(node);
        if(canonical){
            map[parent_id] = canonical;
        } else {
            cache.StoreOrNode();
            map[parent_id] = node;
            node = cache.CreateOrNode();
        }

        ++parent_id;
    } while(xary.Increment());
}

template<bool STRUCTURE>
void MultiGraph::CreateAndNodes(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode){
    // prepare counters
    XAry xary;
    xary.SetDimension(bn_,kSpanningNode->spanning);
    xary.SetZero();

    const SpanningTreeNode *kChildSpanningNodeBegin = &(kSpanningNode->children[0]);
    const auto kNrChildren = kSpanningNode->children.size();
    XAry xary_child[kNrChildren];
    //std::vector<XAry> xary_child(kNrChildren);
    auto kChildSpanningNode = kChildSpanningNodeBegin;
    unsigned int child = 0;
    while(child < kNrChildren){
        xary_child[child].SetDimension(bn_,kChildSpanningNode->spanning);
        ++kChildSpanningNode;
        ++child;
    }

    // store common variables
    const auto kVariable = kSpanningNode->variable;
    const auto kDimension = bn_->get_states()[kVariable];
    const auto kIndex = kSpanningNode->index;
    bool expand_child[kNrChildren];
    const bool kExpand = Expand(kSpanningNode, kVariable, expand_child);
    assert(kNrChildren * kDimension < 10000 && "requesting a lot of stack memory, adjust limit?");

    // layer caching
    CacheLayer &cache = cache_[kIndex];
    auto *map = &(cache.map[0]);

    // store AND node connections before creating actual node
    Node* and_to[kNrChildren][kDimension];
    Node* or_to[kDimension];

    unsigned int parent_id = 0;
    Node *and_node = NULL;
    Node *node = cache.CreateOrNode();
    do {
        // determine connections for AND nodes
        child = 0;
        kChildSpanningNode = kChildSpanningNodeBegin;
        while(child < kNrChildren){
            AddEdges(and_to[child],xary, xary_child[child], kChildSpanningNode, kDimension, expand_child[child]);
            ++kChildSpanningNode;
            ++child;
        }

        // create #kDimension AND nodes
        unsigned int dim = 0;
        while(dim < kDimension){
            if(and_node == NULL)
                and_node = cache.CreateAndNode();

            and_node->variable = kVariable;
            and_node->SetAnd();

            child = 0;
            while(child < kNrChildren){
                and_node->edges[child].to = and_to[child][dim];
                ++child;
            }

            // canonical check and caching
            auto *canonical = Canonical<STRUCTURE>(and_node);
            if(canonical){
                or_to[dim] = canonical;
            } else {
                cache.StoreAndNode();
                or_to[dim] = and_node;
                and_node = NULL;
            }

            ++dim;
        }

        // create OR node with AND node children
        node->variable = kVariable;

        dim = 0;
        while(dim < kDimension){
            node->edges[dim].to = or_to[dim];
            ++dim;
        }
        AddWeights(node, kSpanningTree, kSpanningNode, xary, kDimension);

        // canonical check and caching
        auto *canonical = Canonical<STRUCTURE>(node);
        if(canonical){
            map[parent_id] = canonical;
        } else {
            cache.StoreOrNode();
            map[parent_id] = node;
            node = cache.CreateOrNode();
        }

        ++parent_id;
    } while(xary.Increment());
}


template<bool STRUCTURE>
void MultiGraph::CreateNodes(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningNode){

    AllocateCacheLayer(kSpanningNode,true);

    if(kSpanningNode->children.size() > 1)
        CreateAndNodes<STRUCTURE>(kSpanningTree, kSpanningNode);
    else CreateOrNodes<STRUCTURE>(kSpanningTree, kSpanningNode);
}

template <bool STRUCTURE>
MultiGraph::Node* MultiGraph::Canonical(MultiGraph::Node *node){
    if(STRUCTURE){
        //#ifdef DEBUG
        //auto hash = NodeHash(node);
        //bool equal = false;
        //for(auto it = table_.begin(); it != table_.end(); it++)
        //    if((*it)->operator==(*node))
        //        equal = true;
        //#endif

        return table_.find_or_insert(node);
    }

    return NULL;
}

template<bool STRUCTURE>
MultiGraph::Node* MultiGraph::Compile(const SpanningTree *kSpanningTree, const SpanningTreeNode *kSpanningRoot, bool traversed[]){
    assert(!kSpanningRoot->IsRoot());

    // compile given spanning tree (based on pseudo tree)
    std::stack<const SpanningTreeNode*> s;
    s.push(kSpanningRoot);
    #ifdef DEBUG
    size_t kNrVariables = kSpanningTree->GetNrVariables();
    size_t nr_variables = 0;
    Timer timer;
    #endif
    while(!s.empty()){
        // traverse depth first
        const auto *kSpanningNode = s.top();
        if(!traversed[kSpanningNode->index]){
            traversed[kSpanningNode->index] = true;
            if(!kSpanningNode->IsLeaf())
                for(auto it = kSpanningNode->children.rbegin(); it != kSpanningNode->children.rend(); it++)
                    s.push(&(*it));

            continue;
        }
        s.pop();

        #ifdef DEBUG
        timer.Start();
        #endif
        CreateNodes<STRUCTURE>(kSpanningTree, kSpanningNode);
        #ifdef DEBUG
        timer.Stop();
        printf("Variable %3lu [%3lu/%3lu] %10lu %.3lfms\n", kSpanningNode->variable, ++nr_variables, kNrVariables, kSpanningNode->cardinality, timer.GetDuration<Timer::Milliseconds>());
        #endif
    }

    assert(cache_[kSpanningRoot->index].map.GetSize() == 2);
    Node *kRoot = cache_[kSpanningRoot->index].map[0];
    return kRoot;
}

template <bool STRUCTURE>
void MultiGraph::Compile(const SpanningTree &kSpanningTree){
    const unsigned int kNrLayers = kSpanningTree.GetSize();

    // initialize graph by preallocating levels
    cache_.Resize(kNrLayers);
    terminal_ = cache_.GetTerminal();

    // init computed table
    table_.reserve(kSpanningTree.GetUpperbound());
    table_.max_load_factor(INFINITY);

    bool traversed[kNrLayers] = {0};

    auto *kSpanningNode = &kSpanningTree.GetRoot();
    assert(kSpanningNode->IsRoot());
    if(kSpanningNode->children.size() > 1){
        assert(false && "disconnected components not yet supported");
    } else if(kSpanningNode->children.size() == 1){
        root_ = Compile<STRUCTURE>(&kSpanningTree, &(kSpanningNode->children[0]), traversed);
    }
}

void MultiGraph::DumpDot(std::vector<MultiGraph>& graphs){
    if(graphs.size() == 1){
        graphs[0].DumpDot();
    } else {
        for(unsigned int i = 0; i < graphs.size(); i++)
            graphs[i].DumpDot(i);
    }
}

void MultiGraph::DumpDot(const int kPartitionId){
    auto& files = manager_->get_files();

    std::string aux;
    if(kPartitionId >= 0)
        aux = stringf("%d",kPartitionId);

    DumpDot(files.get_filename_c(DOT,aux));
}

void MultiGraph::DumpDotOrdering(FILE *file){
    // dump variable ordering
    fprintf(file, "\n    // ==================== ORDER ====================\n\n");
    {
        std::vector< std::vector< Variable > > level_to_variable;
        std::map<Variable,size_t> variable_to_level;
        std::queue<const SpanningTreeNode*> s;
        s.push(&(spanningtree_.GetRoot()));

        // print pseudo tree variables
        fprintf(file, "    { // pseudo nodes\n");
        fprintf(file, "        node [shape = plaintext];\n");
        fprintf(file, "\n");

        size_t level = 0;
        bool add_delim = true;
        while(!s.empty()){
            const SpanningTreeNode *kSpanningNode = s.front();
            s.pop();

            if(!kSpanningNode){
                level++;
                add_delim = true;
            } else {
                if(add_delim){
                    add_delim = false;
                    s.push(NULL);
                }

                if(!kSpanningNode->IsTerminal()){
                    variable_to_level[kSpanningNode->variable] = level;
                    for(auto it =kSpanningNode->children.begin(); it != kSpanningNode->children.end(); it++)
                        s.push(&(*it));

                    // store variables per level
                    level_to_variable.resize(level+1);
                    level_to_variable[level].push_back(kSpanningNode->variable);

                    // print variables
                    if(kSpanningNode->variable == std::numeric_limits<Variable>::max()){
                        fprintf(file, "        variable_%u [ label = \"Ordering\"];\n", kSpanningNode->variable);
                    } else {
                        std::string node_name = bn_->get_node_name(kSpanningNode->variable);
                        fprintf(file, "        variable_%u [ label = \"%u:%s\"];\n", kSpanningNode->variable, kSpanningNode->variable,node_name.c_str());
                    }
                }
            }
        }
        fprintf(file, "        \"TERMINAL\" [style = invis];\n");
        fprintf(file, "    }\n\n");

        // print left-right
        for(unsigned int level = 0; level < level_to_variable.size(); level++){
            if(level_to_variable[level].size() > 1){

                fprintf(file, "    { // pseudo rank\n");
                fprintf(file, "        rank = same;\n");
                fprintf(file, "        rankdir = LR;\n\n");
                for(unsigned int i = 0; i < level_to_variable[level].size(); i++){
                    fprintf(file, "        variable_%u",level_to_variable[level][i]);
                    if(i == level_to_variable[level].size()-1)
                        fprintf(file, " [ style=invis ];\n");
                    else fprintf(file, " ->\n");
                }
                fprintf(file, "    }\n\n");
            }
        }

        // print edges
        fprintf(file, "    { // pseudo edges\n");
        //fprintf(file, "        edge [style = invis];\n");
        s.push(&spanningtree_.GetRoot());
        while(!s.empty()){
            const SpanningTreeNode *kSpanningNode = s.front();
            s.pop();

            if(kSpanningNode->IsLeaf()){
                for(auto it = kSpanningNode->children.begin(); it != kSpanningNode->children.end(); it++)
                    fprintf(file, "        variable_%u -> \"TERMINAL\" [style = invis];\n",kSpanningNode->variable);
            } else if(!kSpanningNode->IsTerminal()){
                for(auto it = kSpanningNode->children.begin(); it != kSpanningNode->children.end(); it++){
                    if(kSpanningNode->variable == std::numeric_limits<Variable>::max())
                        fprintf(file, "        variable_%u -> variable_%u [style = invis];\n", kSpanningNode->variable, it->variable);
                    else
                        fprintf(file, "        variable_%u -> variable_%u;\n", kSpanningNode->variable, it->variable);
                    s.push(&(*it));
                }
            }
        }
        fprintf(file, "    }\n\n");
    }
};

void MultiGraph::DumpDotNodes(FILE *file, std::map< Variable, std::vector<const MultiGraph::Node*>>& and_nodes, std::vector< std::vector<const MultiGraph::Node*> >& nodes_per_level, std::vector< std::set< Variable > >& variables_per_level){

    //fprintf(file, "\n    // ==================== NODES ===================\n\n");

    // dump AND nodes
    for(auto and_it = and_nodes.begin(); and_it != and_nodes.end(); and_it++){
        // node rank
        fprintf(file, "    {   // and nodes\n");
        fprintf(file, "        rank = same;\n");
        fprintf(file, "        node [ shape = circle,label = \"*\" ];\n");
        const auto &kNodes = and_it->second;
        for(auto it = kNodes.begin(); it != kNodes.end(); it++){
            const Node *kNode = *it;
            fprintf(file, "        %u;\n", (size_t) kNode);
        }
        fprintf(file, "    }\n\n");

        if(kNodes.size() > 1){
            fprintf(file, "    {    // and nodes\n");
            fprintf(file, "        rank = same;\n");
            fprintf(file, "        rankdir = LR;\n");
            for(auto it = kNodes.begin(); it != kNodes.end(); it++){
                const Node *kNode = *it;
                if(it != kNodes.begin())
                    fprintf(file, "->");
                fprintf(file, "\n");
                fprintf(file, "        %lu", (size_t) kNode);
            }
            fprintf(file, " [ style=invis ];\n");
            fprintf(file, "    }\n\n");
        }
    }

    // dump OR nodes
    for(unsigned int level = 0; level < nodes_per_level.size(); level++){
        const auto &kNodes = nodes_per_level[level];

        if(kNodes.size() == 0)
            continue;

        // node rank
        fprintf(file, "    {   // nodes for level %u\n",level);
        fprintf(file, "        rank = same;\n");
        for(auto it = variables_per_level[level].begin(); it != variables_per_level[level].end(); it++)
            fprintf(file, "        variable_%u;\n", *it);

        for(auto it = kNodes.begin(); it != kNodes.end(); it++){
            const Node *kNode = *it;
            const Variable kVariable = kNode->GetVariable();
            fprintf(file, "        %lu [ label = \"+\\n%s\" ];\n", (size_t) kNode, bn_->get_node_name(kVariable).c_str());
        }
        fprintf(file, "    }\n\n");

        // node order in level
        if(kNodes.size() > 1){
            fprintf(file, "    {    // node order in level %u\n",level);
            fprintf(file, "        rank = same;\n");
            fprintf(file, "        rankdir = LR;\n");
            for(auto it = kNodes.begin(); it != kNodes.end(); it++){
                const Node *kNode = *it;
                if(it != kNodes.begin())
                    fprintf(file, "->");
                fprintf(file, "\n");
                fprintf(file, "        %lu", (size_t) kNode);
            }
            fprintf(file, " [ style=invis ];\n");
            fprintf(file, "    }\n\n");
        }
    }
    fprintf(file, "    {   // TERMINAL\n");
    fprintf(file, "        rank = same;\n");
    fprintf(file, "        \"TERMINAL\";\n");
    fprintf(file, "        %lu [shape=square,label=\"T\"];\n", (size_t) terminal_);
    fprintf(file, "    }\n\n");
}

void MultiGraph::DumpDotEdges(FILE *file, std::map< Variable, std::vector<const MultiGraph::Node*>>& and_nodes, std::vector< std::vector<const MultiGraph::Node*> >& nodes_per_level){
    // dump AND edges
    fprintf(file, "\n    // ==================== EDGES ===================\n\n");
    bayesgraph& g = manager_->get_bayesgraph();
    for(auto and_it = and_nodes.begin(); and_it != and_nodes.end(); and_it++){
        const auto &kNodes = and_it->second;
        for(auto it = kNodes.begin(); it != kNodes.end(); it++){
            const Node * kNode = *it;
            for(unsigned int value = 0; value < kNode->size; value++){
                fprintf(file,"    %lu -> %lu;\n",(size_t) kNode, (size_t) kNode->edges[value].to);
            }
        }
    }

    // dump OR edges
    for(unsigned int level = 0; level < nodes_per_level.size(); level++){
        const auto &kNodes = nodes_per_level[level];
        for(auto it = kNodes.begin(); it != kNodes.end(); it++){
            const Node * kNode = *it;
            for(unsigned int value = 0; value < kNode->size; value++){
                std::string value_name = bn_->get_node_value_name(kNode->variable,value);
                std::string weights;
                for(unsigned int i = 0; i < kNode->edges[value].size; i++){
                    if(i > 0)
                        weights += ",";

                    auto w = kNode->edges[value].weights[i];
                    if(g.defined()){
                        if(OPT_DETERMINISM && (w == 1 || w == 0))
                            weights += stringf("W%u(%.2lf)", w, (double) w);
                        else
                            weights += stringf("W%u(%.2lf)", w, g.get_weight_to_probability()[w-(g.get_nr_literals()+1)]);
                    } else  weights += stringf("W%u", w);
                }

                fprintf(file, "    %lu -> %lu [fontsize=8,label=\"%s\\n%s\"];\n",(size_t) kNode, (size_t) kNode->edges[value].to,value_name.c_str(),weights.c_str());
            }
        }
    }

}

void MultiGraph::DumpDotPrepare(std::map< Variable, std::vector<const MultiGraph::Node*>>& and_nodes, std::vector< std::vector<const MultiGraph::Node*> >& nodes_per_level, std::vector< std::set< Variable > >& variables_per_level){

    std::queue<const MultiGraph::Node*> s;
    s.push(GetRoot());
    size_t level = 0;
    bool add_delim = true;

    std::unordered_map<size_t,bool> done;
    done[(size_t) terminal_] = true;

    while(!s.empty()){
        const auto *kNode = s.front();
        const size_t kNodeId = (size_t) kNode;
        s.pop();

        if(kNode == NULL){
            level++;
            add_delim = true;
        } else if(!done[kNodeId]){
            done[kNodeId] = true;

            if(add_delim){
                add_delim = false;
                s.push(NULL);
            }

            if(level >= nodes_per_level.size()){
                nodes_per_level.resize(level+1);
                variables_per_level.resize(level+1);
            }

            if(kNode->IsAnd()){
                and_nodes[kNode->variable].push_back(kNode);
            } else {
                variables_per_level[level].insert(kNode->variable);
                nodes_per_level[level].push_back(kNode);
            }

            for(unsigned int i = 0; i < kNode->size; i++){
                assert(kNode->edges[i].to);
                s.push(kNode->edges[i].to);
            }
        }
    }
    nodes_per_level.resize(nodes_per_level.size()+1);
    nodes_per_level[nodes_per_level.size()-1].push_back(terminal_);
    variables_per_level.resize(nodes_per_level.size());
}

void MultiGraph::DumpDot(std::string filename){
    std::map< Variable, std::vector<const MultiGraph::Node*>> and_nodes;
    std::vector< std::vector<const MultiGraph::Node*> > nodes_per_level;
    std::vector< std::set< Variable > > variables_per_level;

    DumpDotPrepare(and_nodes,nodes_per_level,variables_per_level);

    FILE *file = fopen(filename.c_str(), "w");
    if(file){
        fprintf(file, "digraph %s {\n", "MULTIGRAPH");
        fprintf(file, "    center = true;\n");
        fprintf(file, "    rankdir = TB;\n");
        fprintf(file, "    edge [arrowsize=.5,dir = normal];\n");

        DumpDotOrdering(file);
        DumpDotNodes(file, and_nodes, nodes_per_level, variables_per_level);
        DumpDotEdges(file, and_nodes, nodes_per_level);

        fprintf(file, "}\n");
        fclose(file);
    } else throw compiler_write_exception("Could not open dot file '%s'",filename.c_str());
}

void MultiGraph::DumpDD(std::vector<MultiGraph> &graphs){
    if(graphs.size() == 1){
        graphs[0].DumpDD();
    } else {
        for(unsigned int i = 0; i < graphs.size(); i++)
            graphs[i].DumpDD(i);
    }
}

void MultiGraph::DumpDD(const int kPartitionId){
    auto& files = manager_->get_files();

    std::string aux;
    if(kPartitionId >= 0)
        aux = stringf("%d",kPartitionId);

    DumpDD(files.get_filename_c(AC,aux));
}

void MultiGraph::DumpDD(std::string filename){
    FILE *file = fopen(filename.c_str(), "w");
    if(file){

        // build numeric index
        std::unordered_map<size_t,size_t> index;
        index[(size_t)terminal_] = 0;   // terminal

        std::unordered_set<size_t> done;
        done.insert((size_t)terminal_);

        const auto *kRoot = GetRoot();
        std::stack<const MultiGraph::Node*> s;
        s.push(kRoot);
        while(!s.empty()){
            const MultiGraph::Node *kNode = s.top();
            const size_t kNodeId = kNode->GetId();
            s.pop();

            if(done.insert(kNodeId).second){
                for(auto edge = kNode->rEdgeBegin(); edge != kNode->rEdgeEnd(); edge--){
                    assert(edge->to);
                    s.push(edge->to);
                }

                index[kNodeId] = index.size()-1;
            }
        }

        // dump nodes
        Size size = GetSize();
        fprintf(file, "multigraph");
        if(spanningtree_.IsTree())
            fprintf(file, " tree");
        else fprintf(file, " chain");
        fprintf(file, " %lu %lu\n", size.nodes+1, size.edges);

        // terminal
        fprintf(file, "0 0\n");

        // format
        // weights : nr_of weights weight_1 ... weight_n
        // edge    : index weights
        // edges   : nr_of_edges edge_1 ... edge_m
        // type    : < + | * >
        // node    : type variable edges
        // line    : node
        done.clear();
        done.insert((size_t) terminal_);
        size_t lines = 1;
        s.push(kRoot);
        while(!s.empty()){
            const auto *kNode = s.top();
            const size_t kNodeId = kNode->GetId();
            s.pop();

            if(done.insert(kNodeId).second){
                for(auto edge = kNode->rEdgeBegin(); edge != kNode->rEdgeEnd(); edge--){
                    assert(edge->to);
                    s.push(edge->to);
                }

                // print node
                if(kNode->IsAnd()){
                    fprintf(file, "*");
                    fprintf(file, " %lu", kNode->GetVariable());    // variable
                } else {
                    fprintf(file, "+");
                    fprintf(file, " %lu", kNode->variable);         // variable
                }

                fprintf(file, " %lu",kNode->size);                  // nr_of_edges
                for(unsigned int i = 0; i < kNode->size; i++){
                    fprintf(file, " %lu", index[(size_t) kNode->edges[i].to]);         // index
                    fprintf(file, " %lu", kNode->edges[i].size);    // nr_of_weights
                    for(unsigned int w = 0; w < kNode->edges[i].size; w++)
                        fprintf(file, " %lu", kNode->edges[i].weights[w]);             // weight_1 ... weight_n
                }
                fprintf(file,"\n");
                lines++;
            }
        }

        fclose(file);
    } else compiler_write_exception("Could not open Multigraph file '%s'", filename.c_str());
}


MultiGraph::Size MultiGraph::GetSize() const {
    Size size = {0,0,0,0,0,0};
    if(root_ == NULL)
        return size;

    std::unordered_map<size_t,bool> done;
    done[(size_t) terminal_] = true;
    std::stack<const Node*> s;
    s.push(GetRoot());
    while(!s.empty()){
        const auto *kNode = s.top();
        const size_t kNodeId = (size_t) kNode;
        s.pop();

        if(!done[kNodeId]){
            done[kNodeId] = true;
            for(int i = kNode->size-1; i >= 0; i--)
                s.push(kNode->edges[i].to);

            // add sizes
            size.edges += kNode->size;
            if(kNode->IsAnd()){
                size.and_nodes++;
                for(int i = kNode->size-1; i >= 0; i--)
                    assert(kNode->edges[i].size == 0);

                size.operators += kNode->size-1;
            } else {
                size.or_nodes++;
                size_t weights = 0;
                for(int i = kNode->size-1; i >= 0; i--)
                    weights += kNode->edges[i].size;

                size.weights += weights;
                size.operators += kNode->size * 2  + weights;
            }
        }
    }
    size.nodes = size.and_nodes + size.or_nodes;

    return size;
}

bool MultiGraph::HasSpanningTree() const {
    return !spanningtree_.Empty();
}

const SpanningTree& MultiGraph::GetSpanningTree() const {
    return spanningtree_;
}

void MultiGraph::InitSpanningChain(bnc::manager_t* m, const bn_partition_t &kBnPartition){
    SetManager(m);
    assert(!initialized_);
    const partition_t &kPartition = kBnPartition.partition;
    const ordering_t &kOrdering = kBnPartition.variable_ordering;

    RuntimeInit(kPartition);
    spanningtree_.Init(kPartition);

    PseudoTree pseudotree;
    pseudotree.SetChain(m->get_bayesnet(), kPartition, kOrdering);
    spanningtree_.Create(pseudotree);
    initialized_ = true;

    #ifdef DEBUG
    Bound bound;
    bound.SetBayesnet(bn_);
    bound.Init(kBnPartition.partition);
    assert(spanningtree_.GetUpperbound()-1 == bound.Compute<Bound::kNodes>(kOrdering));
    #endif
}

void MultiGraph::AllocateCacheLayer(const SpanningTreeNode* const kSpanningNode, const bool kHasMap){
    assert(kSpanningNode->index < cache_.size());

    auto &local_cache = cache_[kSpanningNode->index];
    local_cache.SetNodeProperties(kSpanningNode->GetNrEdgesPerOr(), kSpanningNode->GetNrEdgesPerAnd(), kSpanningNode->GetNrWeightsPerEdge());
    local_cache.Resize(
        (kHasMap?kSpanningNode->GetMapAllocUpperbound():0),
        kSpanningNode->GetOrAllocUpperbound(),
        kSpanningNode->GetOrEdgeAllocUpperbound(),
        kSpanningNode->GetWeightAllocUpperbound(),
        kSpanningNode->GetAndAllocUpperbound(),
        kSpanningNode->GetAndEdgeAllocUpperbound());
}

void MultiGraph::AllocateCache(const bool kHasMap){
    assert(!spanningtree_.Empty() && "Spanning tree was not initialized");

    const unsigned int kNrLayers = spanningtree_.GetSize();
    cache_.Resize(kNrLayers);
    terminal_ = cache_.GetTerminal();

    const auto &kSpanningNodes = spanningtree_.GetSpanningNodes();
    for(auto it = kSpanningNodes.begin(); it != kSpanningNodes.end(); it++)
        AllocateCacheLayer(*it, kHasMap);
}

void MultiGraph::InitSpanningTree(manager* manager, const bn_partition_t &kBnPartition){
    SetManager(manager);
    assert(!initialized_);
    const partition_t &kPartition = kBnPartition.partition;
    const PseudoTree &kPseudoTree = kBnPartition.pseudotree;

    RuntimeInit(kPartition);
    spanningtree_.Init(kPartition);
    spanningtree_.Create(kPseudoTree);
    initialized_ = true;
}

const MultiGraph::Node* MultiGraph::GetRoot() const {
    return root_;
}

template void MultiGraph::Compile<true>();
template void MultiGraph::Compile<false>();

}

