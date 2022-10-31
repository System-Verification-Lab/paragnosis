#ifndef OLD_BOUND_H
#define OLD_BOUND_H


#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include <stack>
#include <set>
#include "basicbnc.h"
#include "partition.h"
#include <gmpxx.h>
#include "ordering.h"

struct bbstate {
    void allocate(const unsigned int);
    std::vector<unsigned int> constraint;
    std::vector<unsigned int> lock;
    std::stack<unsigned int> history;
    std::stack<mpz_class> value;
    std::set<unsigned int> span_variables;
    bbstate& operator=(bbstate&);
    bool operator!=(bbstate&);
};

typedef std::vector< std::vector<unsigned int> > map_backup_t;

template < bool EXT = false>
class bound {
    public:
        bound();
        void init(bayesnet*);
        void init(bbstate&, bayesnet*);
        void init(bbstate&);
        void set_partition(const partition_t&, map_backup_t&);
        void unset_partition(const partition_t&, map_backup_t&);

        void reset();
        void reset(bbstate&);
        void condition(unsigned int);
        void condition(bbstate&, unsigned int);
        void undo();
        void undo(bbstate&);
        void undo(unsigned int);
        void undo(bbstate&, unsigned int);
        mpz_class get_value();
        mpz_class get_value(bbstate&);
        void set_value(size_t);
        void set_value(bbstate&, size_t);
    private:
        std::vector< std::vector<unsigned int> > variable_to_constraint;
        std::vector< std::vector<unsigned int> > constraint_to_variable;
        bayesnet *bn;
        bbstate state;
        unsigned int VARIABLES;
};

template< bool EXT = false> using bound_t = bound<EXT>;

#endif

