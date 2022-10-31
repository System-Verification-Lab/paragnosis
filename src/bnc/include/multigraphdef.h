#ifndef MULTGRAPHDEF_H
#define MULTGRAPHDEF_H

#include "dynamicarray.h"

namespace bnc {

class MultiGraphDef {
    public:

        struct Size {
            size_t operators;
            size_t weights;
            size_t nodes;
            size_t and_nodes;
            size_t or_nodes;
            size_t edges;
        };


        struct Node;

        typedef unsigned int Weight;

        struct __attribute__((__packed__)) Edge {
            Node *to;
            unsigned int size;
            Weight *weights;

            Weight* WeightEnd() {
                return &(weights[size]);
            }

            const Weight* WeightEnd() const {
                return &(weights[size]);
            }

            Weight* WeightBegin() {
                return &(weights[0]);
            }

            const Weight* WeightBegin() const {
                return &(weights[0]);
            }

            Weight* rWeightEnd() {
                return &(weights[-1]);
            }

            const Weight* rWeightEnd() const {
                return &(weights[-1]);
            }

            Weight* rWeightBegin() {
                return &(weights[size-1]);
            }

            const Weight* rWeightBegin() const {
                return &(weights[size-1]);
            }
        };

        struct __attribute__((__packed__)) Node {
            Variable variable;
            uint16_t size;
            Edge *edges;

            static const Variable kTypeMask = 1 << (sizeof(Variable)*8-1);

            Edge* EdgeEnd() {
                return &(edges[size]);
            }

            const Edge* EdgeEnd() const {
                return &(edges[size]);
            }

            Edge* EdgeBegin() {
                return &(edges[0]);
            }

            const Edge* EdgeBegin() const {
                return &(edges[0]);
            }

            Edge* rEdgeEnd() {
                return &(edges[-1]);
            }

            const Edge* rEdgeEnd() const {
                return &(edges[-1]);
            }

            Edge* rEdgeBegin() {
                return &(edges[size-1]);
            }

            const Edge* rEdgeBegin() const {
                return &(edges[size-1]);
            }

            inline size_t GetId() const {
                return (size_t) this;
            }

            inline void UnsetAnd(){
                variable &= ~kTypeMask;
            }

            inline void SetAnd() {
                variable |= kTypeMask;
            }

            inline bool IsAnd() const {
                return variable & kTypeMask;
            }

            inline Variable GetVariable() const {
                return variable & ~kTypeMask;
            }

            inline void SetVariable(Variable v) {
                variable = (variable & kTypeMask) | v;
            }

            inline bool operator==(const Node &other) const {
                if(variable != other.variable || size != other.size)
                    return false;

                const auto kEdgeEnd = &(edges[size]);
                auto *edge = &(edges[0]);
                auto *otheredge = &(other.edges[0]);
                while(edge != kEdgeEnd){
                    if(edge->to != otheredge->to || edge->size != otheredge->size)
                        return false;

                    const auto kWeightEnd = &(edge->weights[edge->size]);
                    auto *weight = &(edge->weights[0]);
                    auto *otherweight = &(otheredge->weights[0]);
                    while(weight != kWeightEnd){
                        if(*weight != *otherweight)
                            return false;

                        ++otherweight;
                        ++weight;
                    }

                    ++otheredge;
                    ++edge;
                }
                return true;
            }
        };

        class CacheLayer {
            public:
                DynamicArray<Node*> map;

            private:
                DynamicArray<Node> nodes;
                DynamicArray<Edge> edges;
                DynamicArray<unsigned int> weights;
                DynamicArray<Node> and_nodes;
                DynamicArray<Edge> and_edges;

                size_t node_top;
                size_t and_node_top;
                size_t edge_top;
                size_t weight_top;

                size_t edges_per_or;
                size_t edges_per_and;
                size_t weights_per_edge;

            public:
                CacheLayer() : node_top(0), and_node_top(0), edge_top(0), weight_top(0), edges_per_or(0), edges_per_and(0), weights_per_edge(0) {};

                inline void SetNodeProperties(const size_t kOrEdges, const size_t kAndEdges, const size_t kWeightsPerEdge){
                    edges_per_or = kOrEdges;
                    edges_per_and = kAndEdges;
                    weights_per_edge = kWeightsPerEdge;
                }

                inline void Resize(const size_t kNrMapNodes, const size_t kOrUpperbound, const size_t kOrEdgeUpperbound, const size_t kWeightUpperbound, const size_t kAndUpperbound, const size_t kAndEdgeUpperbound){
                    map.Resize(kNrMapNodes);
                    nodes.Resize(kOrUpperbound);
                    edges.Resize(kOrEdgeUpperbound);
                    weights.Resize(kWeightUpperbound);
                    and_nodes.Resize(kAndUpperbound);
                    and_edges.Resize(kAndEdgeUpperbound);
                }

                inline Node* CreateNode(const size_t kNrEdges, const size_t kNrWeights){
                    #ifdef DEBUG
                    assert(node_top+1 <= nodes.GetSize() && "failed to allocate node");
                    assert(edge_top+kNrEdges <= edges.GetSize() && "failed to allocate edges");
                    assert(weight_top+kNrEdges * kNrWeights <= weights.GetSize() && "failed to allocate weights");
                    #endif

                    auto *node = &(nodes[node_top]);
                    node->size = kNrEdges;
                    node->edges = &(edges[edge_top]);

                    auto edge_weights = &(weights[weight_top]);
                    const auto kEnd = &(node->edges[kNrEdges]);
                    auto edge = node->edges;
                    while(edge != kEnd){
                        edge->size = kNrWeights;
                        edge->weights = edge_weights;
                        edge_weights += kNrWeights;
                        ++edge;
                    }
                    return node;
                }

                inline Node* CreateOrNode(const size_t kIndex){
                    #ifdef DEBUG
                    assert(!(edges_per_or == 0 && edges_per_and == 0 && weights_per_edge == 0) && "node properties not initialized");
                    assert(kIndex < nodes.GetSize() && "failed to allocate OR node");
                    assert((kIndex*edges_per_or)+edges_per_or <= edges.GetSize() && "failed to allocate OR edges");
                    assert((weights_per_edge == 0 || ((kIndex*edges_per_or*weights_per_edge) + (edges_per_or*weights_per_edge) <= weights.GetSize())) && "failed to allocate weights");
                    #endif

                    auto *node = &(nodes[kIndex]);
                    node->size = edges_per_or;
                    node->edges = &(edges[kIndex*edges_per_or]);

                    auto * edge_weights = &(weights[kIndex*edges_per_or*weights_per_edge]);
                    const auto kEnd = &(node->edges[edges_per_or]);
                    auto edge = node->edges;
                    while(edge != kEnd){
                        edge->size = weights_per_edge;
                        edge->weights = edge_weights;
                        edge_weights += weights_per_edge;
                        ++edge;
                    }
                    return node;
                }

                inline Node* GetOrNode(const size_t kIndex){
                    return &(nodes[kIndex]);
                }

                inline Node* CreateAndNode(const size_t kIndex){
                    #ifdef DEBUG
                    assert(!(edges_per_or == 0 && edges_per_and == 0 && weights_per_edge == 0) && "node properties not initialized");
                    assert(kIndex < and_nodes.GetSize() && "failed to allocate AND node");
                    assert((kIndex*edges_per_and)+edges_per_and <= and_edges.GetSize() && "failed to allocate AND edges");
                    #endif

                    auto *node = &(and_nodes[kIndex]);
                    node->size = edges_per_and;
                    node->edges = &(and_edges[kIndex*edges_per_and]);

                    const auto kEnd = &(node->edges[edges_per_and]);
                    auto edge = node->edges;
                    while(edge != kEnd){
                        edge->size = 0;
                        edge->weights = NULL;
                        ++edge;
                    }
                    return node;
                }

                inline Node* CreateOrNode(){
                    return CreateOrNode(node_top);
                }

                inline Node* CreateAndNode(){
                    return CreateAndNode(and_node_top);
                }

                inline void StoreOrNode(){
                    ++node_top;
                }

                inline void StoreAndNode(){
                    ++and_node_top;
                }

                inline void StoreNode(const size_t kNrEdges, const size_t kNrWeights){
                    ++node_top;
                    edge_top   += kNrEdges;
                    weight_top += kNrEdges * kNrWeights;
                }

                inline Node* CreateTerminal(const size_t kIndex){
                    auto *node = &(nodes[kIndex]);
                    node->variable = 0;
                    node->size = 0;
                    node->edges = NULL;

                    return node;
                }

                inline Node* CreateTerminal(){
                    return CreateTerminal(node_top);
                }

                inline void StoreTerminal(){
                    StoreNode(0,0);
                }

                inline Node* CreateTerminalLayer(){
                    Resize(1,1,0,0,0,0);
                    auto *node = CreateTerminal();
                    StoreTerminal();
                    map[0] = node;

                    return node;
                }

                inline size_t GetMaxNrOrNodes() const {
                    return nodes.GetSize();
                }
        };

        class Cache : public std::vector<CacheLayer> {
            public:
                Cache(){
                    terminal_ = NULL;
                }

                void Resize(size_t size){
                    assert(size > 2 && "Expecting at least an auxiliary root and terminal node");
                    resize(size);
                    auto &terminal_cache = back();
                    terminal_cache.CreateTerminalLayer();
                    terminal_ = terminal_cache.GetOrNode(0);
                }

                Node* GetRoot(){
                    auto it = begin();
                    if(it == end())
                        return NULL;

                    assert(it->GetMaxNrOrNodes() == 0 && "Root is an And node. not implemented.");
                    if(it->GetMaxNrOrNodes() == 0)
                        ++it;

                    if(it == end())
                        return NULL;
                    else return it->GetOrNode(0);
                }

                Node* GetTerminal(){
                    return terminal_;

                }
            private:
                Node *terminal_;
        };
};

}

#endif

