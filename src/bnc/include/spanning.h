#ifndef BNC_SPANNINGTREE_H
#define BNC_SPANNINGTREE_H

#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include <utility>

namespace bnc {

class Spanning {
    public:
        typedef std::vector< std::pair<unsigned int, unsigned int> > EdgeList;

        Spanning();
        void SetBayesnet(bayesnet*);
        void CreateSpanningTree();
        size_t ComputeNrComponents();
        int ComputePartitionMap(const bool * const, unsigned int* const, unsigned int* const, const unsigned int) const;
        const EdgeList& GetEdges() const;
        unsigned int GetSize() const;

    private:
        bayesnet *bn_;
        EdgeList edges_;
};

}

#endif
