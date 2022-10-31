#ifndef BAYESGRAPH_H
#define BAYESGRAPH_H

#include "bayesnode.h"
#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include <unordered_map>
#include <string>

class bayesgraph;

typedef std::unordered_map< std::string, std::unordered_map<std::string, unsigned int> > variable_literal_map_t;

class literal_mapping_t {
    friend class bayesgraph;
    public:
        literal_mapping_t();
        void write() const;
        void read();
        void read(const char*);

        bayesnet *get_bayesnet() const;
        const std::vector<probability_t>& get_weight_to_probability() const;
        const std::vector<unsigned int>& get_variable_to_literal() const;
        const std::vector<unsigned int>& get_literal_to_variable() const;
        const std::vector<unsigned int>& get_dimension() const;
        const variable_literal_map_t& get_variable_value_name_to_literal() const;
        unsigned int get_literal(std::string&, std::string&) const;

        unsigned int get_nr_literals() const;
        unsigned int get_nr_variables() const;
        unsigned int get_nr_weights() const;
        bool is_probability(unsigned int) const;
        probability_t get_probability(unsigned int) const;
        bool defined() const;
    private:
        bayesnet *bn;
        std::vector<probability_t> weight_to_probability;
        std::vector<unsigned int> dimension;
        std::vector<unsigned int> variable_to_literal;
        std::vector<unsigned int> literal_to_variable;

        variable_literal_map_t variable_value_name_to_literal;

        bool DEFINED;
};

class bayesgraph : public literal_mapping_t {
    public:
        int encode(bayesnet*);


        const BayesNode* get_node(unsigned int) const;
        BayesNode* get_node(unsigned int);
        unsigned int size();
        void WriteUAI(std::string = "") ;

    private:
        std::vector<BayesNode> nodes;
};

#endif

