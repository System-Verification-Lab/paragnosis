#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <set>
#include <string>
#include "ordering.h"
#include "support.h"
#include "basicpartition.h"
#include "pseudotree.h"

typedef std::vector< partition_t > partitions_t;

struct bn_partition_t {
    bnc::PseudoTree pseudotree;
    partition_t partition;
    ordering_t ordering;
    ordering_t variable_ordering;
    std::vector< support_order_t > support_to_ordering;
};

class bn_partitions_t : public std::vector<bn_partition_t> {
    public:
        void SetOnePartition(bayesnet *bn);
        void read(std::string, bayesnet* bn = NULL);
        void write(std::string, bayesnet* bn = NULL) const ;
        void verify(bayesnet *bn) const;
        void VerifyOrdering(bayesnet *bn) const;
        bool monolithic();
        unsigned int variable_to_partition_id(variable_t);
};

#endif
