#ifndef BRANCHNBOUND_H
#define BRANCHNBOUND_H


#include "bound.h"

#include <vector>
#include <set>
#include "ordering.h"
#include <bn-to-cnf/bayesnet.h>
#include <queue>

struct bbnode {
    unsigned int variable;
    unsigned int children;
    bbnode *parent;
    std::vector<bbnode*> child;
    size_t upper_bound;
};

struct bbnode_ascending {
    bool operator()(const bbnode* x, const bbnode* y) {
        return x->upper_bound > y->upper_bound;
    };
};

typedef std::priority_queue< bbnode*, std::vector<bbnode*>, bbnode_ascending > bbqueue;

class branchnbound {
    public:
        branchnbound();
        ~branchnbound();
        unsigned int get_ancestor_variable(const unsigned int K, bbnode *node);
        bbnode* get_ancestor(const unsigned int K, bbnode *node);
        bool is_ancestor(const unsigned int K, const unsigned int variable, bbnode *node);
        size_t get_treewidth(const unsigned int, const unsigned int, const unsigned int);
        size_t get_treewidth(const unsigned int, const unsigned int);
        void init(bayesnet*,const partition_t &kPartition);
        void bestfirst();

        void depthfirst(std::queue<bbnode*> &queue, std::set<unsigned int> variables, bbnode *root, bbstate &state, const unsigned int K);
        void depthfirst(std::queue<bbnode*> &queue, const unsigned int VARIABLES, const unsigned int K);
        void lookahead(const unsigned int);
        void lookahead(std::queue<bbnode*> &queue, const unsigned int VARIABLES, const unsigned int K);
        void bestguess();
        void destroy(bbnode*);
        void print();
        ordering_t get_ordering();
    private:
        unsigned int counter;
        bbqueue queue;
        bound_t<false> bound;
        bayesnet *bn;
        const partition_t *kPartition;
        std::vector<ordering_t> orders;
        bbnode *root;
        size_t MAX_UPPER_BOUND;
        std::vector<bbnode*> solutions;
};


#endif

