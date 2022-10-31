#ifndef MULTI_GRAPH_P_H
#define MULTI_GRAPH_P_H

#include <bn-to-cnf/bayesnet.h>
#include "partition.h"
#include <vector>
#include "bnc.h"
#include "spanningtree.h"
#include "xary.h"
#include "multigraphptable.h"
#include <map>
#include "dynamicarray.h"
#include "timer.h"

namespace bnc {

class MultiGraphProbability : public MultiGraphProbabilityDef {
    public:

        MultiGraphProbability();

        void InitSpanningChain(manager*, const bn_partition_t &);
        void InitSpanningTree(manager*, const bn_partition_t &);

        void SetManager(manager*);

        template<bool STRUCTURE>
        void Compile();

        void DumpDot(const int kPartitionId = -1);
        void DumpDot(std::string filename);
        static void DumpDot(std::vector<MultiGraphProbability>&);

        void DumpDD(const int kPartitionId = -1);
        template <bool kBinary = true>
        void DumpDD(std::string filename);
        static void DumpDD(std::vector<MultiGraphProbability>&);

        Size GetSize() const; // O(n)
        bool HasSpanningTree() const;
        const SpanningTree& GetSpanningTree() const;

        const Node *GetRoot() const;

        static void ParallelCompile(std::vector<MultiGraphProbability> &, const unsigned int kNrThreads, Timer &);

        Cache &GetCache(){ return cache_; };

    private:
        static void ParallelScheduler(std::vector<MultiGraphProbability> &, const unsigned int);

        void AllocateCache(const bool kHasMap);
        void AllocateCacheLayer(const SpanningTreeNode* const kSpanningNode, const bool kHasMap);

        void DumpDotPrepare(
            std::map< Variable, std::vector<const MultiGraphProbability::Node*>>&,
            std::vector< std::vector<const MultiGraphProbability::Node*> >&,
            std::vector< std::set< Variable > >&);

        void DumpDotOrdering(FILE *);
        void DumpDotNodes(FILE *,
            std::map< Variable, std::vector<const MultiGraphProbability::Node*>>&,
            std::vector< std::vector<const MultiGraphProbability::Node*> >&,
            std::vector< std::set< Variable > >&);;
        void DumpDotEdges(FILE *,
            std::map< Variable, std::vector<const MultiGraphProbability::Node*>>&,
            std::vector< std::vector<const MultiGraphProbability::Node*> >&);

        void RuntimeInit(const partition_t &kPartition);

        template<bool STRUCTURE>
        void CreateNodes(
            const SpanningTree      *kSpanningTree,
            const SpanningTreeNode  *kSpanningNode);

        template<bool STRUCTURE>
        void CreateOrNodes(
            const SpanningTree      *kSpanningTree,
            const SpanningTreeNode  *kSpanningNode);

        template<bool STRUCTURE>
        void CreateAndNodes(
            const SpanningTree      *kSpanningTree,
            const SpanningTreeNode  *kSpanningNode);

        void AddProbabilities(
            MultiGraphProbability::Node *node,
            const SpanningTree          *kSpanningTree,
            const SpanningTreeNode      *kSpanningNode,
            const XAry                  &xary,
            const unsigned int          kDimension);

        void AddEdges(
            Node                    *node,
            const XAry              &xary,
            XAry                    &xary_child,
            const SpanningTreeNode  *kChildSpanningNode,
            const unsigned int      kDimension,
            const bool              kExpand);

        void AddEdges(
            Node                    **to,
            const XAry              &xary,
            XAry                    &xary_child,
            const SpanningTreeNode  *kChildSpanningNode,
            const unsigned int      kDimension,
            const bool              kExpand);

        template <bool STRUCTURE>
        MultiGraphProbability::Node* Canonical(MultiGraphProbability::Node *node);

        template<bool STRUCTURE>
        Node* Compile(
            const SpanningTree      *kSpanningTree,
            const SpanningTreeNode  *kSpanningRoot,
            bool traversed[]);

        template<bool STRUCTURE>
        void Compile(const bnc::SpanningTree&);

        bool initialized_;

        bnc::manager_t *manager_;
        bayesnet *bn_;
        bayesgraph *g_;
        ordering_t ordering_;

        bnc::SpanningTree spanningtree_;

        void ComputeSpanning(const ordering_t&);

        Cache cache_;
        Node *terminal_;
        Node *root_;
        // compilation variables
        std::vector < DynamicArray<size_t> > node_map_;
        std::vector<const Probability*> cpt_probabilities_;
        std::vector<XAry> cpt_ctr_;
        MultiGraphProbabilityComputedTable table_;
};


}

#endif
