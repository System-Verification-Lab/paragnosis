#ifndef PSEUDOTREE_H
#define PSEUDOTREE_H

#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include "ordering.h"
#include "basicpartition.h"
#include "types.h"
#include <sstream>

namespace bnc {

struct PseudoTreeNode {
    Variable variable;
    std::vector<PseudoTreeNode> children;
};

class PseudoTree {
    public:
        PseudoTree();
        static ordering_t MinFill(bayesnet *bn, const partition_t &);
        void SetTree(bayesnet *bn, const partition_t &, const ordering_t &);
        void SetChain(bayesnet *bn, const partition_t &, const ordering_t &);

        void GenerateTree(bayesnet *bn, const partition_t &);

        const PseudoTreeNode& GetPseudoTree() const;
        const ordering_t& GetPreOrdering() const ;
        const ordering_t& GetMinFillOrdering() const ;
        const size_t GetSize() const ;
        const size_t GetNrVariables() const;
        void Write() const;
        void Write(std::string filename) const;
        void Read(std::string filename);
        void Read();
        bool Empty() const;
        void Clear();

        size_t GetWidth() const;
        size_t GetHeight() const;
    private:
        void CreateTreeString(const PseudoTreeNode* node, std::ostringstream& oss) const;


        size_t width_;
        size_t height_;
        PseudoTreeNode root_;
        ordering_t dfs_ordering_;     // pre-order ordering (dfs traversal of speudotree)
        ordering_t minfill_ordering_; // ordering used to create speudotree
};

}

#endif

