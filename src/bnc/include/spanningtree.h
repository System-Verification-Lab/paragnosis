#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H

#include "pseudotree.h"
#include "spanningtreenode.h"
#include <bn-to-cnf/bayesnet.h>
#include "manager.h"
#include "bayesgraph.h"
#include <vector>
#include <stdint.h>
#include "basicbnc.h"
#include <limits>
#include "xary.h"
#include "types.h"

namespace bnc {


class SpanningTree {
    public:
        SpanningTree();

        bool Empty() const ;
        void SetBayesnet(bayesnet*);
        void SetBayesgraph(bayesgraph*);
        void Init();
        void Init(const partition_t &kPartition);

        void Create(const PseudoTree &);

        size_t GetOrUpperbound() const;
        size_t GetAndUpperbound() const;
        size_t GetUpperbound() const;
        size_t GetOrEdgeUpperbound() const;
        size_t GetAndEdgeUpperbound() const;
        size_t GetEdgeUpperbound() const;
        size_t GetWeightUpperbound() const;
        size_t GetMaxCardinality() const;
        double GetAvgCardinality() const;
        size_t GetHeight() const;
        size_t GetWidth() const;

        size_t GetSize() const;
        size_t GetNrVariables() const;

        const SpanningTreeNode &GetRoot() const;
        const ordering_t& GetPreOrdering() const;
        const std::vector<Variable> & GetMessages() const;

        std::vector<XAry::Mapping> cpt_map_;
        std::vector<XAry::Map> cpt_variable_pos_;

        bool IsTree() const ;

        void Write(int partition_id = -1) const;
        void Write(std::string filename) const;
        void PrintStats() const;

        double RequiredMb(const bool kProbabilityRep = false) const;

        inline const SpanningTreeNode* GetSpanningNode(unsigned int i) const  { return nodes_[i]; }
        inline const SpanningNodes& GetSpanningNodes() const { return nodes_; }

    private:
        void Create(const PseudoTreeNode &kPseudoNode, SpanningTreeNode &node);

        bayesnet *bn_;
        bayesgraph *g_;

        bool is_tree_;
        size_t nr_variables_;
        std::vector<SpanningTreeNode*> terminals_;
        std::vector<Variable> messages_;
        SpanningTreeNode root_;
        size_t or_upperbound_;
        size_t and_upperbound_;
        size_t or_edge_upperbound_;
        size_t and_edge_upperbound_;
        size_t weights_upperbound_;
        size_t size_;
        size_t max_cardinality_;
        size_t width_;
        size_t height_;
        ordering_t pre_ordering_;

        std::vector< bool > enabled_;
        std::vector< bool > message_added_;
        std::vector< std::vector<unsigned int> > variable_to_constraint_;
        std::vector< std::vector<Variable> > constraint_to_variable_;
        SpanningNodes nodes_;
};

}

#endif
