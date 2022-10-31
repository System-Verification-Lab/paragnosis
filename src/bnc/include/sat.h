#ifndef SAT_H
#define SAT_H

#include <bn-to-cnf/cnf.h>
#include <map>
#include <stack>
#include <vector>
#include <set>
#include <queue>
#include "bayesgraph.h"
#include "basicbnc.h"
#include "satisfy.h"

typedef std::pair<literal_t,unsigned int> history_p;

class sat {
    public:
        sat();
        ~sat();
        template <bool DETERMINISM = false> void init(bayesgraph&, int variable = -1);
        void undo();
        bool weighted();
        template <bool DETERMINISM = false> satisfy_t condition(literal_t);
        bnc::node::weights* get_weights();
        void print(bool hist = false);

        void set_manager(manager*);
    private:
        template <bool DETERMINISM = false> void init(BayesNode*);
        void undo(literal_t, unsigned int);
        void release_constraint(literal_t);
        satisfy_t impose_constraint(literal_t);
        satisfy_t condition(literal_t, std::map<literal_t,unsigned int>&, unsigned int unit_clause = 0);

        // limit variables
        unsigned int CLAUSES;
        unsigned int LITERALS;
        unsigned int VARIABLES;
        unsigned int WEIGHTS;

        // information storage variables
        std::vector< unsigned int > dimension;
        std::vector< unsigned int > literals_left;
        std::vector< unsigned int > literal_to_variable;
        std::vector< unsigned int > variable_to_literal;
        std::map<literal_t, literal_t> literal_map;
        std::vector <weight_t> clause_to_weight;
        std::vector< std::vector<literal_t> > literal_to_clause;

        // run time variables
        std::vector< bool > literal_constraint;
        std::vector< std::vector<literal_t> > history;
        std::vector<bool> sat_variable;
        std::vector<uint8_t> unsat_literals;
        std::vector<literal_t> sat_literal_clause;
        unsigned int total_unsat_clauses;
        unsigned int total_unsat_variables;

        std::vector< std::multiset<literal_t> > clauses;
        std::map< literal_t, unsigned int> unit_clauses;

        bool mapped;

        bnc::node::weights *W;
        bnc::manager_t *manager;
};

typedef sat sat_t;

#endif

