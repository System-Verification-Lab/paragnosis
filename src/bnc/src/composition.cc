#include "composition.h"
#include <cassert>
#include <algorithm>
#include <list>
#include <limits>
#include "xary.h"

Composition::Node::Node() : parent(NULL), id(0), partition_id(0) {
}

Composition::Composition(){
    bn_ = NULL;
    composition_ordering_ = NULL;
}

bool Composition::HasPartitions() const {
    return partition_to_shared_variables_.size() > 0;
}

bool Composition::HasOrdering() const {
    return composition_ordering_;
}

size_t Composition::Size() const {
    return nodes_.size();
}

Composition::~Composition(){
    Composition::Destroy(composition_ordering_);
}

void Composition::Destroy(Composition::Node *node){
    if(node != NULL){
        for(auto it = node->children.begin(); it != node->children.end(); it++)
            Destroy(*it);

        delete node;
    }
}

const std::vector<Variable> &Composition::GetPartitionContext(unsigned int partition_id) const {
    #ifdef DEBUG
    assert(partition_id < partition_to_node_.size() && "Out of bounds");
    #endif
    return partition_to_node_[partition_id]->context;
}

const std::vector<variable_t>& Composition::GetPartitionShared(unsigned int partition_id) const {
    #ifdef DEBUG
    assert(partition_id < partition_to_shared_variables_.size() && "Out of bounds");
    #endif

    return partition_to_shared_variables_[partition_id];
}


const std::vector<Variable> &Composition::GetContext(unsigned int i) const {
    #ifdef DEBUG
    assert(i < nodes_.size() && "Out of bounds");
    #endif
    return nodes_[i]->context;
}

void Composition::SetBayesnet(bayesnet *bn){
    this->bn_ = bn;
}

ordering_t Composition::FindOrdering() const{
    assert(bn_ != NULL && "Bayesian network not set for composition ordering");
    assert(variable_to_occurrences_.size() > 0 && "Composition must be initialized first");
    assert(partition_to_shared_variables_.size() <= 12 && "Too many partitions to brute force good ordering");

    // create initial ordering
    const unsigned int kNrPartitions = partition_to_shared_variables_.size();
    ordering_t ordering;
    for(unsigned int i = 0; i < kNrPartitions; i++)
        ordering.push_back(i);

    // init
    ordering_t found_ordering = ordering;
    double best_score = ComputeScore(ordering);
    unsigned int kMaxIterations = 10000000;
    unsigned int iteration = 1;

    // find a componosition ordering with brute force
    while( std::next_permutation(ordering.begin(),ordering.end())){
        double current_score = ComputeScore(ordering);
        if(current_score < best_score){
            best_score = current_score;
            found_ordering = ordering;
        }

        iteration++;
        if(iteration >= kMaxIterations){
            fprintf(stderr, "Warning: brute force limited to ~%u iterations\n", kMaxIterations);
            break;
        }
    }

    return found_ordering;
}


void Composition::SetPartitions(const bn_partitions_t &kPartitions){
    assert(bn_ != NULL && "Bayesian network not set for composition ordering");

    // store variable count from each partition (count > 1 == persistence variable)
    const unsigned int kVariables = bn_->get_nr_variables();
    variable_to_occurrences_.resize(kVariables);
    std::fill(variable_to_occurrences_.begin(), variable_to_occurrences_.end(),0);
    for(unsigned int partition_id = 0; partition_id < kPartitions.size(); partition_id++){
        const partition_t &partition = kPartitions[partition_id].partition;
        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++)
            variable_to_occurrences_[*itv]++;
        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++)
            variable_to_occurrences_[*itv]++;
    }

    // determine shared variables for each partition
    partition_to_shared_variables_.resize(kPartitions.size());
    for(unsigned int partition_id = 0; partition_id < kPartitions.size(); partition_id++){
        const partition_t &partition = kPartitions[partition_id].partition;
        auto &shared = partition_to_shared_variables_[partition_id];

        shared.clear();
        for(auto itv = partition.set.begin(); itv != partition.set.end(); itv++){
            const variable_t &kVariable = *itv;
            if(variable_to_occurrences_[kVariable] > 1)
                shared.push_back(kVariable);
        }
        for(auto itv = partition.cutset.begin(); itv != partition.cutset.end(); itv++){
            const variable_t &kVariable = *itv;
            if(variable_to_occurrences_[kVariable] > 1)
                shared.push_back(kVariable);
        }
        std::sort(shared.begin(), shared.end()); // necessary for 'IsRelated' function
    }
}

bool IsRelated(const std::vector<variable_t>& kShared, const std::set<Variable> &kContext){
    auto first1 = kShared.begin();
    const auto last1 = kShared.end();

    auto first2 = kContext.begin();
    const auto last2 = kContext.end();

    while (first1!=last1 && first2!=last2){
        if(*first1<*first2)
            ++first1;
        else if(*first2<*first1)
            ++first2;
        else {
            return true;
        }
    }
    return false;
}

double Composition::ComputeBestScore(bayesnet *bn, const bn_partitions_t &kPartitions){
    Composition comp;
    comp.SetBayesnet(bn);
    comp.SetPartitions(kPartitions);
    ordering_t ordering = comp.FindOrdering();
    return comp.ComputeScore(ordering);
}

double Composition::ComputeScore(const ordering_t &kOrdering, const bool kCreateChain) const {
    assert(variable_to_occurrences_.size() > 0 && "Composition must be initialized first");

    const auto kDimensions = bn_->get_states();
    double score = 0;

    unsigned int variable_occurrences[variable_to_occurrences_.size()];
    std::copy(variable_to_occurrences_.begin(), variable_to_occurrences_.end(),variable_occurrences);
    std::list< std::set<Variable> > roots;

    const auto kEnd = kOrdering.rend();
    auto it = kOrdering.rbegin();
    while(it != kEnd){
        const auto kPartitionId = *it;
        const auto &kShared = partition_to_shared_variables_[kPartitionId];

        // the new partition node will be the parent of the partition nodes with which it has variables in common
        const auto kRootEnd = roots.end();
        auto itr = roots.begin();
        if(kCreateChain){
            if(itr == kRootEnd){
                roots.emplace_front();
                itr = roots.begin();
            }
        } else {
            unsigned nr_children = 0;
            while(itr != kRootEnd){
                if(IsRelated(kShared,*itr)){
                    // reuse first found child as new node, to eliminate redundant copying of context.
                    ++nr_children;

                    // continue finding other children
                    auto itr2 = itr;
                    ++itr2;
                    while(itr2 != kRootEnd){
                        if(IsRelated(kShared,*itr2)){
                            itr->insert(itr2->begin(), itr2->end());
                            roots.erase(itr2++);
                            ++nr_children;
                        } else
                            ++itr2;
                    }
                    break;
                }
                ++itr;
            }

            // itr points to either nothing or the root of the ordering.
            // thus, create new root if necessary
            if(nr_children == 0){
                roots.emplace_front();
                itr = roots.begin();
            }
        }

        // update context
        for(auto variable_it = kShared.begin(); variable_it != kShared.end(); variable_it++){
            const Variable kVariable = *variable_it;
            assert(variable_occurrences[kVariable] > 0 && "occurance counter goes below 0");
            if(--variable_occurrences[kVariable] != 0)
                itr->insert(kVariable);
            else itr->erase(kVariable);
        }

        // compute score of current node
        double width_score = 1;
        for(auto itc = itr->begin(); itc != itr->end(); itc++)
            width_score *= kDimensions[*itc];

        // compute score
        score += width_score;

        it++;
    }

    return score;
}

size_t Composition::GetTotalNrPartitions() const {
    return partition_to_shared_variables_.size();
}

void Composition::SetNodeIds(Composition::Node * const node, unsigned int &id) const {
    node->id = id;
    for(auto it = node->children.begin(); it != node->children.end(); it++){
        ++id;
        SetNodeIds(*it, id);
    }
}

void Composition::SetNodeIds(Composition::Node * const node) const {
    assert(node && node->parent == NULL);

    unsigned int id = 0;
    SetNodeIds(node,id);
}

Composition::Node* Composition::CreateOrdering(const ordering_t &kOrdering, const bool kCreateChain) const {
    assert(variable_to_occurrences_.size() > 0 && "Composition must be initialized first");

    unsigned int variable_occurrences[variable_to_occurrences_.size()];
    std::copy(variable_to_occurrences_.begin(), variable_to_occurrences_.end(),variable_occurrences);
    std::vector< std::set<Variable> > context(GetTotalNrPartitions());

    std::list< Node* > roots;
    const auto kEnd = kOrdering.rend();
    auto it = kOrdering.rbegin();
    while(it != kEnd){
        const auto kPartitionId = *it;
        const auto &kShared = partition_to_shared_variables_[kPartitionId];

        Node *node = new Node();
        node->partition_id = kPartitionId;

        // the new partition node will be the parent of the partition nodes with which it has variables in common
        if(kCreateChain){
            if(!roots.empty()){
                Node* &child = *roots.begin();
                child->parent = node;
                context[node->partition_id].insert(context[child->partition_id].begin(), context[child->partition_id].end());
                node->children.push_back(child);
                child = node;
            } else roots.push_front(node);
        } else {
            const auto kRootEnd = roots.end();
            auto itr = roots.begin();
            while(itr != kRootEnd){
                Node *child = *itr;
                if(IsRelated(kShared,context[child->partition_id])){
                    (*itr)->parent = node;
                    context[node->partition_id].insert(context[child->partition_id].begin(), context[child->partition_id].end());
                    node->children.push_back(child);
                    roots.erase(itr++);
                } else ++itr;
            }
            roots.push_front(node);
        }

        // update context
        for(auto variable_it = kShared.begin(); variable_it != kShared.end(); variable_it++){
            const Variable kVariable = *variable_it;
            assert(variable_occurrences[kVariable] > 0 && "occurance counter goes below 0");
            if(--variable_occurrences[kVariable] != 0)
                context[node->partition_id].insert(kVariable);
            else context[node->partition_id].erase(kVariable);
        }
        std::copy(context[node->partition_id].begin(), context[node->partition_id].end(),  std::back_inserter(node->context));

        ++it;
    }

    // create common root
    Node *root = new Node();
    root->partition_id = std::numeric_limits<decltype(Node::partition_id)>::max();
    for(auto it = roots.begin(); it != roots.end(); it++){
        (*it)->parent = root;
        root->children.push_back(*it);
    }

    SetNodeIds(root);

    return root;
}

void Composition::BuildNodeMap(const Composition::Node * const node){
    if(node){
        if(!IsRoot(node))
            partition_to_node_[node->partition_id] = node;

        assert(node->id == nodes_.size() && "Node id should have been assigned in pre-order");
        nodes_.push_back(node);

        for(auto it = node->children.begin(); it != node->children.end(); it++)
            BuildNodeMap(*it);
    }
}

void Composition::BuildOrdering(const ordering_t &kOrdering, const bool kCreateChain) {
    ordering_ = kOrdering;
    Composition::Destroy(composition_ordering_);
    composition_ordering_ = CreateOrdering(kOrdering,kCreateChain);

    partition_to_node_.resize(GetTotalNrPartitions());
    std::fill(partition_to_node_.begin(),partition_to_node_.end(), nullptr);
    nodes_.clear();
    BuildNodeMap(composition_ordering_);
}

const ordering_t &Composition::GetOrdering() const{
    return ordering_;
}

const Composition::Node* Composition::GetCompositionOrdering() const{
    return composition_ordering_;
}

const std::vector<const Composition::Node*>& Composition::GetPartitionNodes() const{
    return partition_to_node_;
}

const std::vector<const Composition::Node*>& Composition::GetNodes() const{
    return nodes_;
}

const std::vector< std::vector<variable_t> > &Composition::GetShared() const{
    return partition_to_shared_variables_;
}

#define  BOX_LOWERLEFTCORNER  "\u2514" // └  U+2514
#define  BOX_VERTICALT        "\u251C" // ├  U+251C
#define  BOX_HORIZONTAL       "\u2500" // ─  U+2500
#define  BOX_VERTICAL         "\u2502" // │  U+2502
#define  BOX_HORIZONTALT      "\u252C" // ┬  U+252C

void Print( const Composition::Node * const kNode, const Composition::Node * const kParentNode, std::string &output, bayesnet *bn){
    if(Composition::IsRoot(kNode))
        printf("(*)\n");
    else {
        assert(kParentNode);
        assert(!kParentNode->children.empty());

        std::string corner;
        if(kParentNode->children.back() == kNode)
            corner = BOX_LOWERLEFTCORNER;
        else
            corner = BOX_VERTICALT;

        auto node_str = stringf("%s %s%s%s(%u)",output.substr(0,output.length()-4).c_str(), corner.c_str(),BOX_HORIZONTAL,BOX_HORIZONTAL, kNode->partition_id);
        printf("%-60s", node_str.c_str() );


        std::string context = "{";
        for(auto itc = kNode->context.begin(); itc != kNode->context.end(); itc++){
            if(itc != kNode->context.begin())
                context += stringf(",");
            context += stringf("%u",*itc);
        }
        context += "}";
        printf("%s : %lu\n", context.c_str(),bnc::XAry::Max(bn, kNode->context));


    }

    for(auto it = kNode->children.begin(); it != kNode->children.end(); it++){
        const Composition::Node * const kChild = *it;
        output += stringf(" %s  ", (it + 1 != kNode->children.end()) ? BOX_VERTICAL : " " );
        Print(kChild, kNode, output, bn);
        output = output.substr(0,output.length()-4);
    }
}

void Composition::PrintAscii(bayesnet *bn) const {
    std::string output;
    ::Print( GetCompositionOrdering(), NULL, output, bn);




}


