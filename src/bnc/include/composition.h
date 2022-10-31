#ifndef BNC_COMPOSITION_H
#define BNC_COMPOSITION_H

#include <exception/exception.h>
#include <bn-to-cnf/bayesnet.h>
#include "ordering.h"
#include "partition.h"
#include "types.h"
#include <set>
#include <vector>

create_exception(CompositionException);

class Composition {
    public:
        struct Node {
            Node();
            Node* parent;
            unsigned int id;
            unsigned int partition_id;
            std::vector<Variable> context;
            std::vector<Node*> children;
        };

        Composition();
        ~Composition();

        static void Destroy(Node*);
        inline static bool IsRoot(const Node * const kNode) {
            return kNode->parent == NULL;
        }
        inline static bool IsLeaf(const Node * const kNode) {
            return kNode->children.size() == 0;
        }

        bool HasOrdering() const;
        bool HasPartitions() const;

        static double ComputeBestScore(bayesnet*, const bn_partitions_t &);

        void SetBayesnet(bayesnet*);
        void SetPartitions(const bn_partitions_t &);
        size_t GetTotalNrPartitions() const;
        size_t Size() const;

        ordering_t FindOrdering() const;
        double ComputeScore(const ordering_t &, const bool kCreateChain = false) const;
        Node* CreateOrdering(const ordering_t &, const bool kCreateChain = false) const;
        void BuildOrdering(const ordering_t &, const bool kCreateChain = false);

        const Node* GetCompositionOrdering() const;
        const ordering_t& GetOrdering() const;

        const std::vector< std::vector<variable_t> >& GetShared() const;
        const std::vector<const Node*>& GetPartitionNodes() const;
        const std::vector<const Node*>& GetNodes() const;

        const std::vector<Variable>& GetContext(unsigned int) const;
        const std::vector<Variable>& GetPartitionContext(unsigned int) const;
        const std::vector<variable_t>& GetPartitionShared(unsigned int) const;
        void PrintAscii(bayesnet *) const;

    private:
        void BuildNodeMap(const Node * const);
        void SetNodeIds(Node * const, unsigned int &id) const;
        void SetNodeIds(Node * const) const;
        Node* composition_ordering_;
        ordering_t ordering_;

        bayesnet *bn_;
        std::vector<unsigned int>              variable_to_occurrences_;
        std::vector< std::vector<variable_t> > partition_to_shared_variables_;
        std::vector<const Node*>               partition_to_node_;
        std::vector<const Node*>               nodes_;
};

#endif

