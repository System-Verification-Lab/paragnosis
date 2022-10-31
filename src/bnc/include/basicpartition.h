#ifndef BASICPARTITION_H
#define BASICPARTITION_H

#include <bn-to-cnf/cnf.h>
#include <bn-to-cnf/bayesnet.h>

typedef std::set<variable_t> partition_set_t;

struct partition_t {
    unsigned int id;
    partition_set_t set;
    partition_set_t cutset;

    void determine_cutset(bayesnet*);
    size_t size() const {
        return set.size() + cutset.size();
    }
};

#endif
