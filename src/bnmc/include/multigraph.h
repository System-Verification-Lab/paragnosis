#ifndef BNMC_MULTIGRAPH_H
#define BNMC_MULTIGRAPH_H

#include <bnc/multigraphpdef.h>

class MultiGraph : public bnc::MultiGraphProbabilityDef {
    public:
        class Circuit : public CacheLayer {
            public:
                inline bool Resize(const size_t kNrNodes, const size_t kNrEdges){
                    return CacheLayer::Resize(0,kNrNodes,kNrEdges,0,0);
                }
                inline Node* GetRoot(){ return GetNode(root_index_); }
                inline const Node* GetRoot() const { return GetNode(root_index_); }
                inline size_t GetIndex(const Node * const n) const { return n- &(nodes[0]); }
                inline size_t GetNrTerminals() { return nr_terminals_; }
                inline size_t SetNrTerminals(const size_t kNrTerminals) { return nr_terminals_ = kNrTerminals; }
                inline void SetRootIndex(const size_t kRootIndex){ root_index_ = kRootIndex; }
                inline size_t GetRootIndex() const { return root_index_; }
            private:
                size_t nr_terminals_;
                size_t root_index_;
        };

};

/*

#include <cstdint>
#include <cassert>
#include <limits>
#include <vector>
#include <bnc/types.h>

struct MultiNode;

struct __attribute__((__packed__)) MultiEdge {
    MultiNode* to;
    double w;

    inline void SetIndex(const size_t kTo){
        to = (MultiNode*) kTo;
    }
    inline size_t GetIndex() const {
        return (size_t) to;
    }
};

static_assert(sizeof(size_t) == sizeof(MultiNode*), "Invalid casting between size_t to MultiNode*");
static_assert(sizeof(uint16_t) == sizeof(Variable), "Invalid casting between size_t to MultiNode*");

struct __attribute__((__packed__)) MultiNode {
    uint32_t index;
    uint16_t variable;
    uint16_t size;
    MultiEdge edges[];

    static const uint16_t kTypeMask = 1 << 15;

    inline void UnsetAnd(){
        variable &= ~kTypeMask;
    }

    inline void SetAnd() {
        variable |= kTypeMask;
    }

    inline bool IsAnd() const {
        return variable & kTypeMask;
    }

    inline MultiEdge* rBegin(){
        assert(size > 0);
        return &(edges[size-1]);
    }

    inline const MultiEdge* rBegin() const{
        assert(size > 0);
        return &(edges[size-1]);
    }

    inline MultiEdge* rEnd(){
        return (&(edges[0]))-1;
    }

    inline const MultiEdge* rEnd() const {
        return (&(edges[0]))-1;
    }

    inline MultiEdge* Begin(){
        return &(edges[0]);
    }

    inline const MultiEdge* Begin() const{
        return &(edges[0]);
    }

    inline MultiEdge* End(){
        return &(edges[size]);
    }

    inline const MultiEdge* End() const {
        return &(edges[size]);
    }

    inline Variable GetVariable() const {
        return variable & ~kTypeMask;
    }

    inline uint16_t Size() const {
        return size;
    }

    template <class T>
    inline void SetIndex(const T kIndex){
        assert(kIndex <= std::numeric_limits<uint32_t>::max() && "index too large");
        index = kIndex;
    }

    inline uint32_t GetIndex() const {
        return index;
    }

    template<class T>
    inline void Resize(const T kSize){
        assert(kSize <= std::numeric_limits<uint16_t>::max() && "cannot store this many edges");
        size = kSize;
    }

    inline uint16_t GetVariable(){
        return variable;
    }

    template<class T>
    inline void SetVariable(const T kVariable){
        assert(kVariable <= std::numeric_limits<uint16_t>::max() && "variable too large");
        variable = kVariable;
    }

    template<class T>
    static MultiNode &Cast(T &e){
        return *(MultiNode*) &(e);
    }

    template<class T>
    static MultiNode &Cast(T* &e){
        return *(MultiNode*) e;
    }

    static inline size_t Bytes(const unsigned int kNrEdges){
        return sizeof(MultiNode) + (kNrEdges * sizeof(MultiEdge));
    }

    static inline size_t NodeSize(const unsigned int kNrEdges){
        assert(Bytes(kNrEdges) % sizeof(MultiNode) == 0 && "Simple MultiNode pointer arithmetic will fail");
        return Bytes(kNrEdges)/sizeof(MultiNode);
    }

    static inline void Increment(MultiNode *& ptr){
        //ptr += 1 + (ptr->size * (sizeof(MultiEdge)/sizeof(MultiNode)));
        ptr = (MultiNode*) (((uint8_t*)ptr) + Bytes(ptr->size));
    }
};


class MultiGraphCircuit {
    private:
        typedef std::vector< MultiNode > MultiGraph;
        MultiGraph graph;
        size_t nr_nodes_;
        size_t nr_edges_;

    public:
        inline size_t GetNrNodes() const {
            return nr_nodes_;
        }
        inline size_t Size() const {
            return GetNrNodes();
        }
        inline void SetNrNodes(const size_t kNrNodes){
            nr_nodes_ = kNrNodes;
        }
        inline size_t GetNrEdges() const {
            return nr_edges_;
        }
        inline void SetNrEdges(const size_t kNrEdges){
            nr_edges_ = kNrEdges;
        }
        inline std::vector< MultiNode >& GetCircuit(){
            return graph;
        }

        inline MultiNode *GetRoot(){
            assert(graph.size() > kRootIndex && "MultiGraph is missing");
            return &(graph[kRootIndex]);
        }
        inline const MultiNode *GetRoot() const {
            assert(graph.size() > kRootIndex && "MultiGraph is missing");
            return &(graph[kRootIndex]);
        }
        static const size_t kTerminalIndex = 0;
        static const size_t kRootIndex = 1;
};


static_assert((sizeof(MultiEdge)+sizeof(MultiNode)) % sizeof(MultiNode) == 0, "Pointer arithmetic will fail for MultiNodes");
//static_assert(sizeof(MultiNode) % sizeof(MultiNode) == 0, "Pointer arithmetic will fail for MultiNodes");

//typedef std::vector< MultiNode > MultiGraphCircuit;

*/
#endif
