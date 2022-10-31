#ifndef MANAGER_H
#define MANAGER_H

#include <bn-to-cnf/exceptions.h>
#include <bn-to-cnf/cnf.h>
#include <bnc/files.h>
#include "evidence.h"
#include <map>
#include <vector>
#include <string>

namespace bnmc {

struct QueryProbabilities {
    probability_t p;
    probability_t pq;
};

struct Manager {
        Manager();
        ~Manager();

        void Read(file_t);
        filename_t files;
        std::string logfile;

        bool have_tdmultigraph;
        bool have_multigraph;
        bool have_pwpbdd;
        bool have_partition;
        bool have_wpbdd;
        bool have_ace;
        bool have_ordering;
        bool have_mapping;
        bool have_bn;
        bool have_filename;
        bool have_evidence;
        size_t total_query_count;

        // parameters for parallel execution
        unsigned int buffer;
        std::vector<unsigned int> workers;
        bool OPT_PARALLEL_PARENTS;
        bool OPT_PARALLEL_JOINTS;
        bool OPT_PARALLEL_COMPONENTS;
        probability_t probability;

        Evidence evidence;
        bayesnet *bn;
        literal_mapping_t mapping;
        std::map<std::string, unsigned int> traversals;
        std::map<std::string, QueryProbabilities> probabilities;
};

}

#endif
