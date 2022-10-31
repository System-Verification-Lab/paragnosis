#ifndef ARCHITECTUREBOUND_H
#define ARCHITECTUREBOUND_H

#include "partition.h"
#include "ordering.h"
#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include <exception/exception.h>
#include <unordered_map>
#include <set>

create_exception(ArchBoundException);

class ArchitectureBound {
    public:
        ArchitectureBound();

        void SetBayesnet(bayesnet*);
        void SetPartitions(const bn_partitions_t &,const std::set<unsigned int>&);
        void SetPartitions(const bn_partitions_t &);
        size_t ComputeBound(const ordering_t &);
    private:
        bayesnet *bn_;
        std::vector< std::vector<variable_t> > component_variables_;
        std::vector<unsigned int> total_variable_occurrences_;
        std::vector<unsigned int> copy_;
};

#endif

