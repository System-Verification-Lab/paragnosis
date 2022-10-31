#ifndef DOMAINCLOSURE_H
#define DOMAINCLOSURE_H

#include <bn-to-cnf/cnf.h>
#include "bayesgraph.h"
#include "satisfy.h"
#include <vector>

class domain_closure {
    public:
        void init(bayesgraph *g);
        satisfy_t condition(literal_t);
        void undo(literal_t);
        domain_closure& operator=(domain_closure&);
        bool is_redundant(literal_t);
    private:
        bayesgraph *g;
        std::vector<unsigned int> ALO;
        std::vector<bool> AMO;
        int total_unsat_variables;
};

typedef domain_closure domain_closure_t;

#endif
