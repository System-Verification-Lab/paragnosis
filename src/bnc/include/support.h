#ifndef SUPPORT_H
#define SUPPORT_H

#include <set>
#include <vector>
#include <bn-to-cnf/cnf.h>

typedef std::set<unsigned int> support_set;
typedef std::pair<unsigned int, literal_t> support_pair_t;
typedef std::vector< support_pair_t > support_order_t;
typedef support_set support_t;

struct support_pair_compare {
    bool operator() (const support_pair_t& lhs, const support_pair_t& rhs) const{
        return lhs.first < rhs.first;
    }
};

typedef std::set< support_pair_t,support_pair_compare > support_set_t;

#endif
