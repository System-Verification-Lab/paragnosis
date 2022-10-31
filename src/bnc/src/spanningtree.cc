#include "spanningtree.h"
#include "xary.h"
#include <stack>
#include "options.h"
#include "exceptions.h"
#include "multigraphdef.h"
#include "multigraphpdef.h"
namespace bnc {

SpanningTree::SpanningTree() {
    bn_ = NULL;
    g_ = NULL;

    or_upperbound_ = 0;
    and_upperbound_ = 0;
    size_ = 2;
    is_tree_ = false;
    nr_variables_ = 0;
    max_cardinality_ = 0;
    height_ = 0;
    width_ = 0;
    or_edge_upperbound_ = 0;
    and_edge_upperbound_ = 0;
    weights_upperbound_ = 0;
}

void SpanningTree::SetBayesnet(bayesnet *bn){
    bn_ = bn;

    const auto kNrVariables = bn_->get_nr_variables();
    cpt_map_.resize(kNrVariables);
    cpt_variable_pos_.resize(kNrVariables);
}

void SpanningTree::SetBayesgraph(bayesgraph *g){
    g_ = g;
}


void SpanningTree::Create(const PseudoTreeNode &kPseudoNode, SpanningTreeNode &node){
    const auto kDimensions = bn_->get_states();
    // init node
    node.variable = kPseudoNode.variable;
    node.index = size_++;
    node.dimension = (node.IsRoot()?0:kDimensions[node.variable]);

    nodes_.push_back(&node);
    assert(nodes_.size() == size_);

    if(kPseudoNode.children.size() > 1)
        is_tree_ = true;

    // traverse children
    if(kPseudoNode.children.size() > 0){
        node.children.resize(kPseudoNode.children.size());
        auto it = node.children.begin();
        for(auto pseudo_it = kPseudoNode.children.begin(); pseudo_it != kPseudoNode.children.end(); pseudo_it++, it++)
            Create(*pseudo_it, *it);
    } else {
        // create terminal
        node.children.resize(1);
        terminals_.push_back(&(node.children[0]));
    }

    // introduce messages
    if(!node.IsRoot()){
        for(auto constraint_it = variable_to_constraint_[node.variable].begin(); constraint_it != variable_to_constraint_[node.variable].end(); constraint_it++){
            if(!message_added_[*constraint_it]){
                message_added_[*constraint_it] = true;
                node.messages.push_back(*constraint_it);
            }
        }
    }

    // compute spanning set
    for(auto it = node.children.begin(); it != node.children.end(); it++)
        node.spanning.insert(it->spanning.begin(), it->spanning.end());

    for(auto it = node.messages.begin(); it != node.messages.end(); it++)
        node.spanning.insert(constraint_to_variable_[*it].begin(),constraint_to_variable_[*it].end());

    node.spanning.erase(node.variable);

    // compute cardinalities
    if(node.IsRoot())
        node.cardinality = 0;
    else node.cardinality = XAry::Max(bn_, node.spanning);

    if(max_cardinality_ < node.cardinality)
        max_cardinality_ = node.cardinality;

    // update bounds
    weights_upperbound_ += node.GetWeightUpperbound();
    or_upperbound_ += node.GetOrUpperbound();
    and_upperbound_ += node.GetAndUpperbound();
    or_edge_upperbound_ += node.GetOrEdgeUpperbound();
    and_edge_upperbound_ += node.GetAndEdgeUpperbound();

    assert(node.cardinality < std::numeric_limits<unsigned int>::max() && "Node cardinality intractable");
    assert(node.GetNodeUpperbound() < std::numeric_limits<unsigned int>::max() && "Total number of nodes intractable");

    // init xary info and create mapping from parent to child spanning set(s)
    for(auto it = node.children.begin(); it != node.children.end(); it++){
        SpanningTreeNode &child = *it;
        XAry::CreateMap(node.spanning, child.spanning, child.xary.map);
        XAry::CreateRestrictDimension(child.xary.map, child.xary.restrict_dimension);

        child.xary.dimension.resize(child.spanning.size());
        XAry::SetDimension(&(child.xary.dimension[0]), child.spanning.begin(), child.spanning.end(), bn_->get_states());

        auto parent_variable_it = child.spanning.find(node.variable);
        child.xary.expand = parent_variable_it == child.spanning.end();
        child.xary.pos = (!child.xary.expand?std::distance(child.spanning.begin(),parent_variable_it):child.spanning.size());
    }

    // create spanning to message mappings
    for(auto it = node.messages.begin(); it != node.messages.end(); it++){
        const unsigned int kCpt = *it;
        const BayesNode *n = g_->get_node(kCpt);
        const auto &kCptVariables = n->GetVariables();
        XAry::CreateMap(node.spanning, kCptVariables, cpt_map_[kCpt]);
        cpt_variable_pos_[kCpt] = std::distance(kCptVariables.begin(),std::find(kCptVariables.begin(), kCptVariables.end(), node.variable));
    }

    #ifdef DEBUG
    {
        if(!node.IsRoot()){
            for(auto it = node.children.begin(); it != node.children.end(); it++){
                SpanningTreeNode &child = *it;
                XAry xary;
                xary.SetDimension(bn_,child.spanning);

                if(!(xary.Max(child.xary.restrict_dimension) == 1 || xary.Max(child.xary.restrict_dimension) == node.dimension))
                    assert((xary.Max(child.xary.restrict_dimension) == 1 || xary.Max(child.xary.restrict_dimension) == node.dimension) && "Child node is not restricted as expected");

                unsigned int unrestricted = 0;
                for(unsigned int i = 0; i < child.xary.map.size(); i++){
                    if(!child.xary.restrict_dimension[i]){
                        ++unrestricted;
                        auto spanning_variable_it = child.spanning.begin();
                        std::advance(spanning_variable_it, i);
                        assert(*spanning_variable_it == node.variable && "Unrestricted variable must be parent spanning node variable");
                    }
                }
                assert(unrestricted <= 1 && "Too many unrestricted dimensions");
            }
        }
    }
    #endif
}

void SpanningTree::Create(const PseudoTree &kPseudoTree){
    assert(nr_variables_ == kPseudoTree.GetNrVariables() && "Number of variables in ordering and partition differ" );

    nodes_.clear();
    width_ = kPseudoTree.GetWidth();
    height_ = kPseudoTree.GetHeight();
    std::fill(message_added_.begin(), message_added_.end(), false);
    pre_ordering_ = kPseudoTree.GetPreOrdering();
    size_ = 0;
    is_tree_ = false;
    or_upperbound_ = 0;
    and_upperbound_ = 0;
    or_edge_upperbound_ = 0;
    and_edge_upperbound_ = 0;
    max_cardinality_ = 0;
    weights_upperbound_ = 0;

    Create(kPseudoTree.GetPseudoTree(), root_);

    root_.cardinality = 0;
    for(auto it = terminals_.begin(); it != terminals_.end(); it++){
        auto *terminal = (*it);
        if(nodes_.size() == size_)
            nodes_.push_back(terminal);
        terminal->variable = SpanningTreeNode::kRootVariable-1;
        terminal->cardinality = 1;
        terminal->index = size_;
    }
    size_++;


    assert(nodes_.size()==size_);
    assert(root_.IsRoot());
}

void SpanningTree::Init(){
    assert(bn_ != NULL && "Bayesnet must be set first");
    assert(g_ != NULL && "Bayesgraph must be set first");
    const unsigned int kNrVariables = bn_->get_nr_variables();
    if(nr_variables_ == 0)
        nr_variables_ = kNrVariables;

    messages_.clear();
    variable_to_constraint_.clear();
    constraint_to_variable_.clear();


    if(enabled_.empty())
        enabled_.resize(kNrVariables,true);

    message_added_.resize(kNrVariables);
    variable_to_constraint_.resize(kNrVariables);
    constraint_to_variable_.resize(kNrVariables);

    for(unsigned int cpt = 0; cpt < kNrVariables; cpt++){
        if(enabled_[cpt]){
            messages_.push_back(cpt);
            variable_to_constraint_[cpt].push_back(cpt);
            constraint_to_variable_[cpt].push_back(cpt);

            const uint32_t* parent = bn_->get_parent(cpt);
            const unsigned int kNrParents = bn_->get_parent_size(cpt);
            for(unsigned int i = 0; i < kNrParents; i++){
                variable_to_constraint_[parent[i]].push_back(cpt);
                constraint_to_variable_[cpt].push_back(parent[i]);
            }
        }
    }
}

void SpanningTree::Init(const partition_t &kPartition){
    assert(bn_ != NULL && "Bayesnet must be set first");
    assert(g_ != NULL && "Bayesgraph must be set first");
    const unsigned int kNrVariables = bn_->get_nr_variables();
    nr_variables_ = kPartition.set.size() + kPartition.cutset.size();

    enabled_.resize(kNrVariables);
    std::fill(enabled_.begin(),enabled_.end(),false);
    for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++)
        enabled_[*it] = true;

    Init();
}

void SpanningTree::Write(int partition_id) const {
    std::string filename;
    if(partition_id >= 0)
        filename = files.get_filename(SPANNING, stringf("%lu",partition_id));
    else filename = files.get_filename(SPANNING);

    Write(filename);
}

template <class T>
std::string stl_to_string(const T &input){
    std::string output;
    output += "{";
    bool comma = false;
    for(auto it = input.begin(); it != input.end(); it++){
        if(comma)
            output += ",";
        else comma = true;
        output += stringf("%lu", *it);
    }
    output += "}";
    return output;
}

void SpanningTree::Write(std::string filename) const {
    FILE *file = fopen(filename.c_str(), "w");
    if(file){
        fprintf(file, "digraph %s {\n", "SPANNINGTREE");
        fprintf(file, "    forcelabels = true; // omit no labels\n");
        fprintf(file, "    center = true;\n");
        fprintf(file, "    rankdir = TB;\n");
        fprintf(file, "    edge [arrowsize=.5,dir = normal];\n\n");


        const auto *kRoot = &GetRoot();
        assert(kRoot->IsRoot());

        fprintf(file, "\n    // ==================== NODES ===================\n\n");
        fprintf(file, "    %lu [ label = \"order\", shape = plaintext ];\n", (size_t) kRoot);
        std::queue<const SpanningTreeNode*> q;
        q.push(kRoot);
        while(!q.empty()){
            const SpanningTreeNode *kNode = q.front();
            q.pop();

            for(auto it = kNode->children.begin(); it != kNode->children.end(); it++)
                if(!it->IsTerminal())
                    q.push(&(*it));

            if(kNode == kRoot)
                continue;

            std::string label;
            label += stringf("Variable %lu ",kNode->variable);
            #ifdef DEBUG
            label += stl_to_string(variable_to_constraint_[kNode->variable]);
            #endif

            label += "\\nSpanning: ";
            label += stl_to_string(kNode->spanning);
            label += "\\nCardinality: ";
            label += stringf("%lu", kNode->cardinality);

            #ifdef DEBUG
            std::vector <std::string> messages;
            for(auto it = kNode->messages.begin(); it != kNode->messages.end(); it++){
                std::string constraints = stl_to_string(constraint_to_variable_[*it]);
                std::string message = stringf("m %lu : ",*it);
                message += constraints;
                messages.push_back(message);
            }
            for(auto it = messages.begin(); it != messages.end(); it++){
                label += "\\n";
                label += *it;
            }
            #endif

            fprintf(file, "    %lu [ fontsize=8, label = \"%s\", ];\n", (size_t) kNode, label.c_str());

        }

        fprintf(file, "\n    // ==================== EDGES ===================\n\n");
        q.push(kRoot);
        while(!q.empty()){
            const SpanningTreeNode *kNode = q.front();
            q.pop();

            for(auto it = kNode->children.begin(); it != kNode->children.end(); it++){
                if(!it->IsTerminal()){
                    q.push(&(*it));
                    fprintf(file, "    %lu -> %lu;\n", (size_t) kNode, (size_t) &(*it));
                }
            }
        }

        fprintf(file, "}\n");
        fclose(file);
    } else throw compiler_write_exception("Could not open dot file '%s'",filename.c_str());
}

void SpanningTree::PrintStats() const {
    printf("Spanning tree stats\n");
    printf("===================\n");
    printf("    Is tree           : %s\n",  (IsTree()?"true":"false"));
    printf("    Upperbound        : %lu\n", GetUpperbound());
    printf("    Or Upperbound     : %lu\n", GetOrUpperbound());
    printf("    And Upperbound    : %lu\n", GetAndUpperbound());
    printf("    Weight Upperbound : %lu/%lu (ratio %.2lf)\n", GetWeightUpperbound(),
            bn_->get_nr_probabilities(),
            (double)GetWeightUpperbound()/(double)bn_->get_nr_probabilities());
    printf("    Edge Upperbound   : %lu\n", GetEdgeUpperbound());
    printf("    Max cardinality   : %lu\n", GetMaxCardinality());
    printf("    Avg cardinality   : %.2lf\n", GetAvgCardinality());
    printf("    Variables         : %lu\n", GetNrVariables());
    printf("    Width/height      : %lu/%lu\n", GetWidth(), GetHeight());
    printf("    Estimated p space : %.3lfMb\n", RequiredMb(true));
    printf("    Estimated space   : %.3lfMb\n", RequiredMb(false));
}

double SpanningTree::RequiredMb(const bool kProbabilityRep) const {
    size_t size = 0;
    if(!kProbabilityRep){
        size += sizeof(MultiGraphProbabilityDef::Node) * (or_upperbound_ + and_upperbound_);            // nodes
        size += sizeof(MultiGraphProbabilityDef::Edge) * (or_edge_upperbound_ + and_edge_upperbound_);  // edges
        size += sizeof(MultiGraphProbabilityDef::Node*) * (or_upperbound_ + and_upperbound_);           // computed table
        size += sizeof(MultiGraphProbabilityDef::Node*) * (or_upperbound_);                             // node map
    } else {
        size += sizeof(MultiGraphDef::Node) * (or_upperbound_ + and_upperbound_);           // nodes
        size += sizeof(MultiGraphDef::Edge) * (or_edge_upperbound_ + and_edge_upperbound_); // edges
        size += sizeof(MultiGraphDef::Node*) * (or_upperbound_ + and_upperbound_);          // computed table
        size += sizeof(MultiGraphDef::Node*) * (or_upperbound_);                            // node map
        size += sizeof(MultiGraphDef::Weight) * weights_upperbound_;                        // weights
    }

    double Mb = (double) size / (double) 1000000;
    return Mb;
}

size_t SpanningTree::GetAndUpperbound() const {
    return and_upperbound_;
}

size_t SpanningTree::GetOrUpperbound() const {
    return or_upperbound_;
}

size_t SpanningTree::GetUpperbound() const {
    return GetAndUpperbound() + GetOrUpperbound() + 1; // +1 for terminal
}

size_t SpanningTree::GetNrVariables() const {
    return size_-2;
}

size_t SpanningTree::GetSize() const {
    return size_;
}

const SpanningTreeNode& SpanningTree::GetRoot() const {
    return root_;
}

const ordering_t & SpanningTree::GetPreOrdering() const {
    return pre_ordering_;
}

const std::vector<Variable>& SpanningTree::GetMessages() const {
    return messages_;
}

bool SpanningTree::IsTree() const {
    return is_tree_;
}

bool SpanningTree::Empty() const {
    return root_.IsTerminal();
}

size_t SpanningTree::GetMaxCardinality() const {
    return max_cardinality_;
}

size_t SpanningTree::GetHeight() const {
    return height_;
}

size_t SpanningTree::GetWidth() const {
    return width_;
}

double SpanningTree::GetAvgCardinality() const {
    return (double)GetOrUpperbound() / (double)GetNrVariables();
}

size_t SpanningTree::GetOrEdgeUpperbound() const {
    return or_edge_upperbound_;
}

size_t SpanningTree::GetAndEdgeUpperbound() const {
    return and_edge_upperbound_;
}

size_t SpanningTree::GetEdgeUpperbound() const {
    return GetAndEdgeUpperbound() + GetOrEdgeUpperbound();
}

size_t SpanningTree::GetWeightUpperbound() const {
    return weights_upperbound_;
}

}

