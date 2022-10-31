#ifndef ORDERING_H
#define ORDERING_H

#include <bn-to-cnf/cnf.h>
#include <bn-to-cnf/bayesnet.h>
#include <string>
#include "basicpartition.h"
#include "bayesgraph.h"
#include "files.h"
#include "types.h"

typedef std::vector<std::string> string_ordering_t;

class ordering_t : public std::vector<literal_t> {
    public:
        ordering_t();
        ordering_t& operator=(const ordering_t&);
        ordering_t& operator=(ordering_t&);
        const char* get_filename();
        void write(std::string, bayesnet *bn = NULL) const;
        void read(std::string, bayesnet *bn = NULL);
        void generate_variable_ordering(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_1(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_2(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_3(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_4(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_5(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_6(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_7(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_8(bayesnet *bn,const partition_t &kPartition);
        void generate_variable_ordering_9(bayesnet *bn,const partition_t &kPartition);

        void generate_variable_ordering_simulated_annealing(bayesnet *bn, const partition_t&);

        void literal_to_variable_ordering(bayesgraph &);
        void variable_to_literal_ordering(bayesgraph &, const partition_t &, const bool kIncludeWeightLiterals);
        void set_ordering(const ordering_t&, const BayesNode*, const bool kIncludeWeightLiterals);
        void set_ordering(ordering_t&,ordering_t&,ordering_t&);
        void add_weights(bayesgraph &g, const partition_t &);
        ordering_t get_literal_to_variable_ordering(bayesgraph &);
        ordering_t get_variable_to_literal_ordering(bayesgraph &, const partition_t&, const bool kIncludeWeightLiterals);
        string_ordering_t get_variable_string_ordering(bayesnet *) const ;
        void stats(bayesnet *bn, const partition_t &kPartition);
    private:
        int anneal(bayesnet *bn, const partition_t&, ordering_t &, unsigned int score_function = 0);
};


class biordering_t {
    public:
        ordering_t to;
        ordering_t back;
        void set(const ordering_t&);
};

#endif

