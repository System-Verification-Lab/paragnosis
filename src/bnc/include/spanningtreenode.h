#ifndef BNC_SPANNINGTREENODE_H
#define BNC_SPANNINGTREENODE_H

#include "types.h"
#include <vector>
#include <set>
#include "xary.h"
#include <cstdint>
#include <cassert>

namespace bnc {

typedef std::set<Variable> SpanningSet;

struct XAryInfo {
    unsigned int pos;
    bool expand;
    XAry::Dimension dimension;
    XAry::Mapping map;
    XAry::RestrictDimension restrict_dimension;
};

struct SpanningTreeNode {
    Variable variable;
    uint16_t index;
    uint16_t dimension;
    size_t cardinality;
    SpanningSet spanning;
    std::vector <Variable> messages;
    std::vector<SpanningTreeNode> children;
    XAryInfo xary;

    static const Variable kRootVariable = std::numeric_limits<Variable>::max();

    inline bool IsRoot() const { return variable == SpanningTreeNode::kRootVariable; }
    inline bool IsTerminal() const { return children.size() == 0; }
    inline bool IsLeaf() const { return children.size() == 1 && children[0].IsTerminal(); }
    inline bool IsAndLayer() const { return children.size() > 1; }

    inline size_t GetNrEdgesPerOr() const {
        return dimension;
    }

    inline size_t GetNrEdgesPerAnd() const {
        return children.size();
    }

    inline size_t GetNrWeightsPerEdge() const {
        return messages.size();
    }

    // upperbounds
    inline size_t GetOrUpperbound() const {
        return cardinality;
    }

    inline size_t GetAndUpperbound() const {
        return (GetNrEdgesPerAnd()>1?GetOrUpperbound() * GetNrEdgesPerOr():0);
    }

    inline size_t GetNodeUpperbound() const {
        return GetOrUpperbound() + GetAndUpperbound();
    }

    inline size_t GetOrEdgeUpperbound() const {
        return GetOrUpperbound() * GetNrEdgesPerOr();
    }

    inline size_t GetAndEdgeUpperbound() const {
        return GetAndUpperbound() * GetNrEdgesPerAnd();
    }

    inline size_t GetEdgeUpperbound() const {
        return GetOrEdgeUpperbound() + GetAndEdgeUpperbound();
    }

    inline size_t GetWeightUpperbound() const {
        return GetOrEdgeUpperbound() * GetNrWeightsPerEdge();
    }

    inline size_t GetMapUpperbound() const {
        return GetOrUpperbound();
    }

    // allocations
    inline size_t GetOrAllocUpperbound() const {
        auto size = GetOrUpperbound();
        return (size==0?0:size+1);
    }

    inline size_t GetAndAllocUpperbound() const {
        auto size = GetAndUpperbound();
        return (size==0?0:size+1);
    }

    inline size_t GetNodeAllocUpperbound() const {
        return GetOrAllocUpperbound() + GetAndAllocUpperbound();
    }

    inline size_t GetOrEdgeAllocUpperbound() const {
        auto size = GetOrEdgeUpperbound();
        return (size==0?0:size+GetNrEdgesPerOr());
    }

    inline size_t GetAndEdgeAllocUpperbound() const {
        auto size = GetAndEdgeUpperbound();
        return (size==0?0:size+GetNrEdgesPerAnd());
    }

    inline size_t GetEdgeAllocUpperbound() const {
        return GetOrEdgeAllocUpperbound() + GetAndEdgeAllocUpperbound();
    }

    inline size_t GetWeightAllocUpperbound() const {
        auto size = GetWeightUpperbound();
        return (size==0?0:size+(GetNrWeightsPerEdge()*GetNrEdgesPerOr()));
    }

    inline size_t GetMapAllocUpperbound() const {
        return GetOrAllocUpperbound();
    }

};

typedef std::vector<SpanningTreeNode*> SpanningTreeNodeList;

class SpanningNodes : public SpanningTreeNodeList {
    public:
        const const_iterator last() const {
            #ifdef DEBUG
            assert((begin() == end() || (*(end() - 1))->IsTerminal()) && "No terminal found");
            #endif


            auto it = end();
            if(it == begin())
                return it;

            do {
                --it;
            } while(it != begin() && (*it)->IsTerminal());

            return it;
        }
        const const_iterator first() const {
            #ifdef DEBUG
            assert((begin() == end() || !(*begin())->IsAndLayer()) && "Disconnected components detected!");
            #endif

            auto it = begin();
            if(it != end())
                return ++it;
            else return it;
        }
};


}

#endif

