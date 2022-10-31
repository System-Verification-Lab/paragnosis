#ifndef MULTGRAPHPDEF_H
#define MULTGRAPHPDEF_H

#include "dynamicarray.h"
#include "types.h"

namespace bnc {

class MultiGraphProbabilityDef {
    public:
        struct Size {
            size_t operators;
            size_t probabilities;
            size_t nodes;
            size_t and_nodes;
            size_t or_nodes;
            size_t edges;
        };

        struct Node;


        struct __attribute__((__packed__)) Edge {
            Node *to;
            Probability probability;
        };

        static_assert(sizeof(uint16_t) == sizeof(Variable), "Node struct containt missaligned members");

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

                const auto kEdgeEnd = EdgeEnd();
                auto *edge = EdgeBegin();
                auto *otheredge = other.EdgeBegin();
                while(edge != kEdgeEnd){
                    if(edge->to != otheredge->to || edge->probability != otheredge->probability)
                        return false;

                    ++otheredge;
                    ++edge;
                }
                return true;
            }
        };

        class CacheLayer {
            public:
                DynamicArray<Node*> map;

            protected:
                DynamicArray<Node> nodes;
                DynamicArray<Edge> edges;
                DynamicArray<Node> and_nodes;
                DynamicArray<Edge> and_edges;

                size_t node_top;
                size_t and_node_top;
                size_t edge_top;

                size_t edges_per_or;
                size_t edges_per_and;

            public:
                CacheLayer() : node_top(0), and_node_top(0), edge_top(0), edges_per_or(0), edges_per_and(0) {};

                inline bool Empty() const {
                    return !(nodes.GetSize() != 0 || edges.GetSize() != 0 || and_nodes.GetSize() != 0 || and_edges.GetSize() != 0);
                }

                inline size_t GetSize() const {
                    return node_top;
                }

                inline void SetNodeProperties(const size_t kOrEdges, const size_t kAndEdges){
                    edges_per_or = kOrEdges;
                    edges_per_and = kAndEdges;
                }

                inline bool Resize(const size_t kNrMapNodes, const size_t kOrUpperbound, const size_t kOrEdgeUpperbound, const size_t kAndUpperbound, const size_t kAndEdgeUpperbound){
                    if(!map.Resize(kNrMapNodes))
                        return false;
                    if(!nodes.Resize(kOrUpperbound))
                        return false;
                    if(!edges.Resize(kOrEdgeUpperbound))
                        return false;
                    if(!and_nodes.Resize(kAndUpperbound))
                        return false;
                    if(!and_edges.Resize(kAndEdgeUpperbound))
                        return false;

                    return true;
                }

                inline Node* CreateNode(const size_t kNrEdges){
                    #ifdef DEBUG
                    assert(node_top+1 <= nodes.GetSize() && "failed to allocate node");
                    assert(edge_top+kNrEdges <= edges.GetSize() && "failed to allocate edges");
                    #endif

                    auto *node = &(nodes[node_top]);
                    node->size = kNrEdges;
                    node->edges = &(edges[edge_top]);

                    return node;
                }

                inline Node* CreateOrNode(const size_t kIndex){
                    #ifdef DEBUG
                    assert(!(edges_per_or == 0 && edges_per_and == 0) && "node properties not initialized");
                    assert(kIndex < nodes.GetSize() && "failed to allocate OR node");
                    assert((kIndex*edges_per_or)+edges_per_or <= edges.GetSize() && "failed to allocate OR edges");
                    #endif

                    auto *node = &(nodes[kIndex]);
                    node->size = edges_per_or;
                    node->edges = &(edges[kIndex*edges_per_or]);

                    return node;
                }

                inline Node* GetNode(const size_t kIndex){
                    return &(nodes[kIndex]);
                }

                inline const Node* GetNode(const size_t kIndex) const {
                    return &(nodes[kIndex]);
                }

                inline Node* GetOrNode(const size_t kIndex){
                    return &(nodes[kIndex]);
                }

                inline Node* CreateAndNode(const size_t kIndex){
                    #ifdef DEBUG
                    assert(!(edges_per_or == 0 && edges_per_and == 0) && "node properties not initialized");
                    assert(kIndex < and_nodes.GetSize() && "failed to allocate AND node");
                    assert((kIndex*edges_per_and)+edges_per_and <= and_edges.GetSize() && "failed to allocate AND edges");
                    #endif

                    auto *node = &(and_nodes[kIndex]);
                    node->size = edges_per_and;
                    node->edges = &(and_edges[kIndex*edges_per_and]);

                    const auto kEdgeEnd = node->EdgeEnd();
                    auto edge = node->EdgeBegin();
                    while(edge != kEdgeEnd){
                        edge->probability = 1;
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

                inline void StoreNode(const size_t kNrEdges){
                    ++node_top;
                    edge_top += kNrEdges;
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
                    StoreNode(0);
                }

                inline Node* CreateTerminalLayer(){
                    Resize(1,1,0,0,0);
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
            protected:
                Node *terminal_;
        };

};

}

#endif

